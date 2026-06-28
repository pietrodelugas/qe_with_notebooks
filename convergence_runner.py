"""convergence_runner.py — parsing, running, and sweep helpers for QE convergence notebooks.

Public API
----------
Constants
    RY_TO_EV, BOHR_TO_ANG, RY_PER_BOHR_TO_EV_PER_ANG

Parsers (pure functions — only need the pw.x stdout string)
    extract_qe_error(text)
    parse_total_energy_ry(stdout)
    parse_irreducible_kpoints(stdout)
    parse_force_z_ev_ang(stdout, atom_index_1based)
    parse_total_scf_correction_ev_ang(stdout)
    parse_stress_zz_kbar(stdout)

Convergence criterion (pure function)
    first_globally_converged_index(values, threshold)

Runner class
    QERunner(pw_cmd)
        .run_sweep(cases, run_dir, force_rerun=False,
                   collect_force_stress=False, atom_index_1based=1,
                   collect_scf_correction=False)

        cases is a list of (tag, PWInput) pairs built in the notebook.
"""

from __future__ import annotations

import re
import subprocess
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Unit conversion constants
# ---------------------------------------------------------------------------

RY_TO_EV = 13.605693009          # 1 Ry in eV
BOHR_TO_ANG = 0.529177210903     # 1 bohr in Angstrom
RY_PER_BOHR_TO_EV_PER_ANG = RY_TO_EV / BOHR_TO_ANG

# ---------------------------------------------------------------------------
# Output parsers — pure functions; only need the pw.x stdout string
# ---------------------------------------------------------------------------

def _is_qe_fence(line: str) -> bool:
    s = line.strip()
    return len(s) >= 60 and s == '%' * len(s)


def extract_qe_error(text: str) -> str | None:
    """Return the QE error block from pw.x output, or None if not present.

    QE wraps fatal errors between two lines of '%' characters (the exact
    count varies by QE version — typically 72 or 78).  Any all-'%' line
    of at least 60 characters is treated as a fence.
    """
    lines = text.splitlines()
    fences = [i for i, ln in enumerate(lines) if _is_qe_fence(ln)]
    if len(fences) < 2:
        return None
    block = lines[fences[0] + 1 : fences[1]]
    return '\n'.join(block).strip() or None


def parse_total_energy_ry(stdout: str) -> float:
    """Return the final total energy in Ry from pw.x stdout.

    For SCF there is one ``! total energy`` line; for relax there are several
    (one per ionic step).  Always returns the last match so the function is
    correct for both calculation types.
    """
    hits = re.findall(r'!\s+total energy\s+=\s+([-0-9.]+)\s+Ry', stdout)
    if not hits:
        raise ValueError('Could not parse total energy from pw.x output.')
    return float(hits[-1])


def parse_irreducible_kpoints(stdout: str):
    """Return the number of irreducible k-points, or None if not found."""
    m = re.search(r'number of k points=\s*(\d+)', stdout)
    return int(m.group(1)) if m else None


def parse_force_z_ev_ang(stdout: str, atom_index_1based: int) -> float:
    """Return the z-component of the force on one atom in eV/Ang.

    pw.x reports forces in Ry/bohr; this function converts automatically.

    Parameters
    ----------
    stdout : str
        Full pw.x standard output.
    atom_index_1based : int
        Atom number as printed in the pw.x output (starts from 1).
    """
    pattern = (
        rf'atom\s+{atom_index_1based}\s+type\s+\d+\s+force\s*=\s*'
        r'([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)'
    )
    m = re.search(pattern, stdout)
    if not m:
        raise ValueError(f'Could not parse force for atom {atom_index_1based}.')
    fz_ry_bohr = float(m.group(3))
    return fz_ry_bohr * RY_PER_BOHR_TO_EV_PER_ANG


def parse_hydrostatic_pressure_kbar(stdout: str) -> float:
    """Return the hydrostatic pressure in kbar from the pw.x stress block.

    pw.x prints: ``total   stress  (Ry/bohr**3)    (kbar)     P= <value>``
    The last occurrence is used so the function works for vc-relax outputs too.
    """
    hits = re.findall(r'P=\s*([-+]?\d*\.?\d+)', stdout)
    if not hits:
        raise ValueError('Could not parse hydrostatic pressure from pw.x output.')
    return float(hits[-1])


def parse_mean_abs_scf_correction_ev_ang(stdout: str) -> float:
    """Return the mean absolute per-component SCF force correction in eV/Ang.

    The "Total SCF correction" printed on a single line can vanish by symmetry
    (equal-and-opposite corrections on symmetry-related atoms sum to zero).
    This function instead parses the per-atom block printed when
    verbosity = 'medium' or 'high' and computes::

        sum(|dF_{i,alpha}|) / (3 * nat)

    which is always non-zero when corrections exist.  Raises ValueError if the
    block is absent (set verbosity = 'medium' in CONTROL).
    """
    m = re.search(
        r'The SCF correction term to forces\s*\n((?:[ \t]*atom\s+\d+.*\n)+)',
        stdout,
    )
    if not m:
        raise ValueError(
            'SCF correction block not found — set verbosity = "medium" in CONTROL.'
        )
    components = []
    for line in m.group(1).splitlines():
        hit = re.search(
            r'force\s*=\s*([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)', line
        )
        if hit:
            components.extend(float(hit.group(k)) for k in (1, 2, 3))
    if not components:
        raise ValueError('SCF correction block found but contained no atom lines.')
    mean_abs_ry_bohr = sum(abs(c) for c in components) / len(components)
    return mean_abs_ry_bohr * RY_PER_BOHR_TO_EV_PER_ANG


def parse_stress_zz_kbar(stdout: str) -> float:
    """Return the zz component of the stress tensor in kbar.

    pw.x prints the stress block as:
        total   stress  (Ry/bohr**3)    (kbar)     P= ...
         sxx  sxy  sxz   sxx(kbar) ...
         syx  syy  syz   syx(kbar) ...
         szx  szy  szz   szx(kbar) ... szz(kbar)

    This function extracts the last number on the third row (szz in kbar).
    """
    block = re.search(
        r'total\s+stress.*?\(kbar\)\s*\n\s*(.*?)\n\s*(.*?)\n\s*(.*?)\n',
        stdout,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if not block:
        raise ValueError('Could not parse stress block from pw.x output.')

    nums = re.findall(r'[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?', block.group(3))
    if len(nums) < 6:
        raise ValueError('Stress block found but zz(kbar) could not be parsed.')
    return float(nums[-1])


# ---------------------------------------------------------------------------
# Convergence criterion — pure function
# ---------------------------------------------------------------------------

def first_globally_converged_index(values, threshold) -> int | None:
    """Return the first index i such that all values[i:] are <= threshold.

    This is the strict criterion used in the exercise: a setting is converged
    only when every denser/larger value also stays within the threshold — not
    just the two nearest neighbours.

    Returns None if no such index exists in the array.
    """
    arr = np.asarray(values, dtype=float)
    # tail_max[i] = max(arr[i], arr[i+1], ..., arr[-1])
    tail_max = np.maximum.accumulate(arr[::-1])[::-1]
    idx = np.where(tail_max <= threshold)[0]
    return int(idx[0]) if len(idx) else None


# ---------------------------------------------------------------------------
# QERunner — holds run-time context, exposes run_sweep
# ---------------------------------------------------------------------------

class QERunner:
    """Run pre-built pw.x inputs and parse their outputs.

    Parameters
    ----------
    pw_cmd : list[str]
        Command prefix to invoke pw.x, e.g. ['/path/to/pw.x'].

    Usage
    -----
    >>> runner = QERunner(pw_cmd)
    >>> results = runner.run_sweep(cases, run_dir, force_rerun=False)
    """

    def __init__(self, pw_cmd):
        self.pw_cmd = list(pw_cmd)

    def run_sweep(
        self,
        cases,
        run_dir,
        force_rerun=False,
        collect_force_stress=False,
        atom_index_1based=1,
        collect_pressure=False,
        collect_scf_correction=False,
    ):
        """Run all cases and return a list of result dicts.

        Parameters
        ----------
        cases : list of (tag, PWInput)
            Each tag names the .in/.out files; each PWInput is a complete input.
        run_dir : Path
            Directory for input/output files.
        force_rerun : bool
            If True, overwrite existing output files and rerun.
        collect_force_stress : bool
            If True, also parse Fz from each output.
        atom_index_1based : int
            Atom number (1-based) for force parsing.
        collect_pressure : bool
            If True, parse the hydrostatic pressure P (kbar) from the stress block.
            Requires ``tstress = .true.`` in the pw.x input.
        collect_scf_correction : bool
            If True, parse the mean absolute per-component SCF force correction
            (eV/Ang) from the verbose force block.  Requires ``verbosity = 'medium'``
            or ``'high'`` in the CONTROL namelist.

        Returns
        -------
        list of dicts — one per case, with keys:
            tag, energy_ry, wall_s, nk_irr,
            and optionally force_z_ev_ang, pressure_kbar, scf_correction_ev_ang.
        """
        results = []
        for i, (tag, inp) in enumerate(cases, 1):
            prefix = f'  [{i}/{len(cases)}] {tag}'
            print(f'{prefix}: running …', flush=True)
            try:
                data = self._run_case(tag, inp, Path(run_dir), force_rerun,
                                      collect_force_stress, atom_index_1based,
                                      collect_pressure, collect_scf_correction)
            except Exception:
                print(f'{prefix}: FAILED', flush=True)
                raise
            if np.isnan(data['wall_s']):
                status = 'cached'
            else:
                status = f'{data["wall_s"]:.1f}s'
            print(f'{prefix}: {status}    ')
            results.append(data)
        return results

    def _run_case(self, tag, inp, run_dir, force_rerun, collect_force_stress,
                  atom_index_1based, collect_pressure=False, collect_scf_correction=False):
        in_file  = run_dir / f'{tag}.in'
        out_file = run_dir / f'{tag}.out'

        in_file.write_text(inp.to_string())

        if out_file.is_file() and not force_rerun:
            _text = out_file.read_text()
        else:
            _text = None
        if _text is not None and 'JOB DONE.' in _text:
            stdout = _text
            wall_s = np.nan
        else:
            t0 = time.perf_counter()
            result = subprocess.run(
                self.pw_cmd + ['-input', str(in_file)],
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
                    f'pw.x failed for {in_file.name} (return code {result.returncode})\n'
                    + detail
                )
            stdout = result.stdout

        data = {
            'tag': tag,
            'energy_ry': parse_total_energy_ry(stdout),
            'wall_s': wall_s,
            'nk_irr': parse_irreducible_kpoints(stdout),
        }
        if collect_force_stress:
            data['force_z_ev_ang'] = parse_force_z_ev_ang(stdout, atom_index_1based)
        if collect_pressure:
            data['pressure_kbar'] = parse_hydrostatic_pressure_kbar(stdout)
        if collect_scf_correction:
            data['scf_correction_ev_ang'] = parse_mean_abs_scf_correction_ev_ang(stdout)
        return data
