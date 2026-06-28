"""
ph_runner.py — single-shot runners for the QE phonon workflow codes.

Public API
----------
Parsers
    parse_ph_frequencies_cm(stdout)     — phonon frequencies in cm⁻¹ from ph.x output
    parse_dynmat_modes(stdout)          — mode table (freq, IR, optional Raman) from dynmat.x
    parse_dynmat_polarizability(stdout) — 3×3 polarizability tensor (Å³) from dynmat.x
    parse_dynmat_eps_electronic(stdout) — 3×3 ε∞ tensor (requires lperm=.true. in &INPUT)
    parse_dynmat_eps_static(stdout)     — 3×3 ε₀ = ε∞ + Δε_ionic (same requirement)

Runner classes
    PhRunner(ph_cmd)
        .run_one(tag, inp, run_dir, force_rerun=False)
        Returns: {tag, wall_s, frequencies_cm}

    Q2rRunner(q2r_cmd)
        .run_one(tag, inp, run_dir, force_rerun=False)
        Returns: {tag, wall_s}

    MatdynRunner(matdyn_cmd)
        .run_one(tag, inp, run_dir, force_rerun=False)
        Returns: {tag, wall_s}

    DynmatRunner(dynmat_cmd)
        .run_one(tag, inp, run_dir, force_rerun=False)
        Returns: {tag, wall_s, freq_cm, freq_thz, ir, raman, polarizability,
                  eps_electronic, eps_static, eps_ionic}

All runners accept any input object with a .to_string() method (PhInputph,
Q2rInput, MatdynInput, or any namelist object).  They write <tag>.in and
<tag>.out to run_dir and skip the calculation if a completed output already
exists (unless force_rerun=True).
"""

from __future__ import annotations

import re
import subprocess
import time
from pathlib import Path

import numpy as np

from convergence_runner import _print, extract_qe_error
from matdyn_results import MatdynResults


# ---------------------------------------------------------------------------
# Output parsers
# ---------------------------------------------------------------------------

def parse_ph_frequencies_cm(stdout: str):
    """Parse phonon frequencies in cm⁻¹ from ph.x stdout.

    Matches lines of the form:
        freq (    N) =   X.XXX [THz] =   Y.YYY [cm-1]

    Negative values (imaginary modes) are preserved as-is.

    Returns
    -------
    ndarray of float, or None if no frequency lines are found.
    """
    pattern = re.compile(
        r'freq\s*\(\s*\d+\s*\)\s*=\s*[-\d.]+\s*\[THz\]\s*=\s*([-\d.]+)\s*\[cm-1\]'
    )
    matches = pattern.findall(stdout)
    return np.array([float(x) for x in matches]) if matches else None


def parse_dynmat_modes(stdout: str) -> dict | None:
    """Parse the mode table from dynmat.x stdout.

    The table header is either:
        ``# mode   [cm-1]    [THz]      IR``
    or (when lraman was active in ph.x):
        ``#  mode   [cm-1]     [THz]       IR      Raman``

    Returns
    -------
    dict with ndarray values, or None if no table found:
        freq_cm  : ndarray[float] — frequencies in cm⁻¹ (negative = imaginary)
        freq_thz : ndarray[float] — frequencies in THz
        ir       : ndarray[float] — IR activity in (D/Å)²/amu
        raman    : ndarray[float] | None — Raman activity in Å⁴/amu, or None
    """
    header = re.search(r'#\s+mode\s+\[cm-1\]\s+\[THz\]\s+IR(\s+Raman)?', stdout)
    if header is None:
        return None
    has_raman = header.group(1) is not None

    # Data rows follow immediately after the header line
    tail = stdout[header.end():]
    # Each row: mode_idx  freq_cm  freq_thz  ir  [raman]
    if has_raman:
        row_pat = re.compile(
            r'^\s+\d+\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$',
            re.MULTILINE,
        )
    else:
        row_pat = re.compile(
            r'^\s+\d+\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$',
            re.MULTILINE,
        )

    rows = row_pat.findall(tail)
    if not rows:
        return None

    arr = np.array([[float(x) for x in r] for r in rows])
    result = {
        'freq_cm':  arr[:, 0],
        'freq_thz': arr[:, 1],
        'ir':       arr[:, 2],
        'raman':    arr[:, 3] if has_raman else None,
    }
    return result


def _parse_3x3_block(stdout: str, marker_re: str):
    """Find the first line matching marker_re, then return the next 3×3 float block.

    Skips non-numeric lines (e.g. the Clausius-Mossotti note that appears before
    the polarizability matrix).  Returns ndarray(3,3) or None.
    """
    m = re.search(marker_re, stdout)
    if m is None:
        return None
    float_row = re.compile(r'^\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$')
    rows = []
    for line in stdout[m.end():].splitlines():
        fm = float_row.match(line)
        if fm:
            rows.append([float(fm.group(1)), float(fm.group(2)), float(fm.group(3))])
            if len(rows) == 3:
                break
    return np.array(rows) if len(rows) == 3 else None


def parse_dynmat_polarizability(stdout: str):
    """Parse the polarizability tensor from dynmat.x stdout.

    Reads the 3×3 matrix (in Å³) printed after ``Polarizability (A^3 units)``.
    The optional Clausius-Mossotti line between the header and the matrix is skipped.

    Returns ndarray(3,3) or None.
    """
    return _parse_3x3_block(stdout, r'Polarizability \(A\^3 units\)')


def parse_dynmat_eps_electronic(stdout: str):
    """Parse the electronic (high-frequency) dielectric tensor ε∞ from dynmat.x stdout.

    Only present when ``lperm = .true.`` is set in the dynmat.x &INPUT namelist.

    Returns ndarray(3,3) or None if the section is absent.
    """
    return _parse_3x3_block(stdout, r'Electronic dielectric permittivity tensor')


def parse_dynmat_eps_static(stdout: str):
    """Parse the static dielectric tensor ε₀ = ε∞ + Δε_ionic from dynmat.x stdout.

    Printed immediately after the electronic tensor when ``lperm = .true.``.

    Returns ndarray(3,3) or None if the section is absent.
    """
    return _parse_3x3_block(
        stdout, r'\.\.\. with zone-center polar mode contributions'
    )


# ---------------------------------------------------------------------------
# Internal helper — shared run logic
# ---------------------------------------------------------------------------

def _run_code(cmd, tag, inp, run_dir, force_rerun, code_name):
    """Write input, run the code, cache on JOB DONE., return (stdout, wall_s)."""
    in_file  = run_dir / f'{tag}.in'
    out_file = run_dir / f'{tag}.out'

    in_file.write_text(inp.to_string())

    cached_text = None
    if out_file.is_file() and not force_rerun:
        cached_text = out_file.read_text()
    if cached_text is not None and 'JOB DONE.' in cached_text:
        return cached_text, np.nan

    t0 = time.perf_counter()
    result = subprocess.run(
        cmd + ['-input', str(in_file)],
        capture_output=True, text=True,
    )
    wall_s = time.perf_counter() - t0
    out_file.write_text(result.stdout)

    if result.returncode != 0:
        print(result.stdout, flush=True)
        if result.stderr.strip():
            print(result.stderr, flush=True)
        qe_msg = extract_qe_error(result.stdout) or extract_qe_error(result.stderr)
        if qe_msg:
            detail = f'QE error:\n{qe_msg}'
        elif result.stderr.strip():
            detail = f'Last stderr:\n{result.stderr[-1200:]}'
        else:
            detail = f'Last stdout:\n{result.stdout[-1200:]}'
        raise RuntimeError(
            f'{code_name} failed for {in_file.name} (return code {result.returncode})\n'
            + detail
        )
    return result.stdout, wall_s


# ---------------------------------------------------------------------------
# PhRunner
# ---------------------------------------------------------------------------

class PhRunner:
    """Run a single ph.x calculation and parse its output.

    Parameters
    ----------
    ph_cmd : list[str]
        Command prefix to invoke ph.x, e.g.
        ``['conda', 'run', '-n', 'qe_env', 'ph.x']``.

    run_one() returns a dict with keys
    ------------------------------------
    tag            : str
    wall_s         : float  — wall time in s (nan if result was cached)
    frequencies_cm : ndarray | None
        Phonon frequencies in cm⁻¹ printed at the end of the run.
        Available for single-q and qplot calculations; None for ldisp
        runs (where frequencies are written to dynamical matrix files).
    """

    def __init__(self, ph_cmd):
        self.ph_cmd = list(ph_cmd)

    def run_one(self, tag, inp, run_dir, force_rerun=False):
        """Run a single ph.x calculation.

        Parameters
        ----------
        tag : str
            Names the .in / .out files written to run_dir.
        inp : PhInputph
            Complete ph.x input object.
        run_dir : str | Path
            Directory for input / output files.
        force_rerun : bool
            If True, overwrite an existing output and rerun.

        Returns
        -------
        dict — see class docstring for keys.
        """
        _print(f'  {tag}: running ph.x …')
        try:
            stdout, wall_s = _run_code(
                self.ph_cmd, tag, inp, Path(run_dir), force_rerun, 'ph.x'
            )
        except Exception:
            _print(f'  {tag}: FAILED')
            raise
        status = 'cached' if np.isnan(wall_s) else f'{wall_s:.1f}s'
        _print(f'  {tag}: {status}')
        return {
            'tag':            tag,
            'wall_s':         wall_s,
            'frequencies_cm': parse_ph_frequencies_cm(stdout),
        }


# ---------------------------------------------------------------------------
# Q2rRunner
# ---------------------------------------------------------------------------

class Q2rRunner:
    """Run a single q2r.x calculation.

    Parameters
    ----------
    q2r_cmd : list[str]
        Command prefix to invoke q2r.x, e.g.
        ``['conda', 'run', '-n', 'qe_env', 'q2r.x']``.

    run_one() returns a dict with keys
    ------------------------------------
    tag    : str
    wall_s : float  — wall time in s (nan if cached)
    """

    def __init__(self, q2r_cmd):
        self.q2r_cmd = list(q2r_cmd)

    def run_one(self, tag, inp, run_dir, force_rerun=False):
        """Run a single q2r.x calculation.

        Parameters
        ----------
        tag : str
        inp : Q2rInput
        run_dir : str | Path
        force_rerun : bool
        """
        _print(f'  {tag}: running q2r.x …')
        try:
            _, wall_s = _run_code(
                self.q2r_cmd, tag, inp, Path(run_dir), force_rerun, 'q2r.x'
            )
        except Exception:
            _print(f'  {tag}: FAILED')
            raise
        status = 'cached' if np.isnan(wall_s) else f'{wall_s:.1f}s'
        _print(f'  {tag}: {status}')
        return {'tag': tag, 'wall_s': wall_s}


# ---------------------------------------------------------------------------
# MatdynRunner
# ---------------------------------------------------------------------------

class MatdynRunner:
    """Run a single matdyn.x calculation.

    Parameters
    ----------
    matdyn_cmd : list[str]
        Command prefix to invoke matdyn.x, e.g.
        ``['conda', 'run', '-n', 'qe_env', 'matdyn.x']``.

    run_one() returns a dict with keys
    ------------------------------------
    tag     : str
    wall_s  : float          — wall time in s (nan if cached)
    results : MatdynResults  — lazy handle for reading and plotting output files

    The ``results`` object gives access to dispersion, DOS, and mode data
    without reading any files until the relevant property is first accessed.
    """

    def __init__(self, matdyn_cmd):
        self.matdyn_cmd = list(matdyn_cmd)

    def run_one(self, tag, inp, run_dir, force_rerun=False):
        """Run a single matdyn.x calculation.

        Parameters
        ----------
        tag : str
        inp : MatdynInput
        run_dir : str | Path
        force_rerun : bool
        """
        _print(f'  {tag}: running matdyn.x …')
        run_dir = Path(run_dir)
        try:
            _, wall_s = _run_code(
                self.matdyn_cmd, tag, inp, run_dir, force_rerun, 'matdyn.x'
            )
        except Exception:
            _print(f'  {tag}: FAILED')
            raise
        status = 'cached' if np.isnan(wall_s) else f'{wall_s:.1f}s'
        _print(f'  {tag}: {status}')
        return {
            'tag':     tag,
            'wall_s':  wall_s,
            'results': MatdynResults.from_input(run_dir, inp),
        }


# ---------------------------------------------------------------------------
# DynmatRunner
# ---------------------------------------------------------------------------

class DynmatRunner:
    """Run a single dynmat.x calculation and parse its output.

    Parameters
    ----------
    dynmat_cmd : list[str]
        Command prefix to invoke dynmat.x, e.g.
        ``['conda', 'run', '-n', 'qe_env', 'dynmat.x']``.

    run_one() returns a dict with keys
    ------------------------------------
    tag            : str
    wall_s         : float         — wall time in s (nan if cached)
    freq_cm        : ndarray       — phonon frequencies in cm⁻¹
    freq_thz       : ndarray       — phonon frequencies in THz
    ir             : ndarray       — IR activities in (D/Å)²/amu
    raman          : ndarray|None  — Raman activities in Å⁴/amu (None if not computed)
    polarizability : ndarray|None  — 3×3 polarizability tensor in Å³
                                     (None if ε∞ / Z* were not on the dynmat file)
    eps_electronic : ndarray|None  — 3×3 ε∞ (electronic) dielectric tensor
                                     (None unless lperm = .true. in &INPUT)
    eps_static     : ndarray|None  — 3×3 ε₀ = ε∞ + Δε_ionic (static dielectric tensor)
                                     (None unless lperm = .true. in &INPUT)
    eps_ionic      : ndarray|None  — 3×3 ionic (lattice) contribution Δε = ε₀ − ε∞
                                     (None unless lperm = .true. in &INPUT)
    """

    def __init__(self, dynmat_cmd):
        self.dynmat_cmd = list(dynmat_cmd)

    def run_one(self, tag, inp, run_dir, force_rerun=False):
        """Run a single dynmat.x calculation.

        Parameters
        ----------
        tag : str
        inp : any object with a .to_string() method
        run_dir : str | Path
        force_rerun : bool
        """
        _print(f'  {tag}: running dynmat.x …')
        try:
            stdout, wall_s = _run_code(
                self.dynmat_cmd, tag, inp, Path(run_dir), force_rerun, 'dynmat.x'
            )
        except Exception:
            _print(f'  {tag}: FAILED')
            raise
        status = 'cached' if np.isnan(wall_s) else f'{wall_s:.1f}s'
        _print(f'  {tag}: {status}')

        modes   = parse_dynmat_modes(stdout) or {}
        eps_el  = parse_dynmat_eps_electronic(stdout)
        eps_st  = parse_dynmat_eps_static(stdout)
        eps_ion = (eps_st - eps_el) if (eps_el is not None and eps_st is not None) else None
        return {
            'tag':            tag,
            'wall_s':         wall_s,
            'freq_cm':        modes.get('freq_cm'),
            'freq_thz':       modes.get('freq_thz'),
            'ir':             modes.get('ir'),
            'raman':          modes.get('raman'),
            'polarizability': parse_dynmat_polarizability(stdout),
            'eps_electronic': eps_el,
            'eps_static':     eps_st,
            'eps_ionic':      eps_ion,
        }
