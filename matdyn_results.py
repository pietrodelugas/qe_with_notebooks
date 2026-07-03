"""
matdyn_results.py — read and analyse the output files produced by matdyn.x.

File readers (standalone functions)
------------------------------------
    read_freq(path)     — .freq file  → (q_points, frequencies)
    read_freq_gp(path)  — .freq.gp file → (path_lengths, frequencies)
    read_dos(path)      — .dos file   → (energies, dos, pdos)
    read_modes(path)    — .modes file → list of per-q dicts

Results class
-------------
    MatdynResults(run_dir, flfrq, fldos, flvec)
        .from_input(run_dir, inp)     — construct from a MatdynInput

    Lazy-loaded data
        .q_points    ndarray (nq, 3)        q-vectors in reciprocal coords
        .frequencies ndarray (nq, nbnd)     phonon frequencies in cm⁻¹
        .path        ndarray (nq,)          cumulative path lengths (from .freq.gp)
        .energies    ndarray (ndos,)        DOS frequency axis in cm⁻¹
        .dos         ndarray (ndos,)        total phonon DOS in states/cm⁻¹
        .pdos        ndarray (nat, ndos)    per-atom PDOS, or None
        .modes       list of dicts          eigenvectors per q-point (from .modes)

    Inspection
        .dos_data()                  → (energies, dos, pdos)  — x/y arrays for the DOS plot
        .dispersion_data()           → (path, frequencies)    — x/y arrays for the dispersion plot
        .frequencies_at(n)           → ndarray (nbnd,) in cm⁻¹
        .mode_at(n, band)            → ndarray (nat, 3) complex  eigenvector

    Plotting
        .plot_dispersion(ax, xticks, ...)
        .plot_dos(ax, with_pdos, ...)
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np


# ============================================================================
# File readers
# ============================================================================

def read_freq(path) -> tuple[np.ndarray, np.ndarray]:
    """Read a matdyn.x ``.freq`` file.

    Format::

        &plot nbnd=N, nks=M /
                qx  qy  qz
           freq1  freq2  …  freqN       ← 6 per line, wraps for nbnd>6
                qx  qy  qz
           …

    Parameters
    ----------
    path : str | Path

    Returns
    -------
    q_points    : ndarray (nq, 3)       q-vectors in reciprocal coords
    frequencies : ndarray (nq, nbnd)    frequencies in cm⁻¹ (negative = imaginary)
    """
    with open(path) as f:
        header = f.readline()
        nbnd = int(re.search(r'nbnd\s*=\s*(\d+)', header).group(1))
        nq   = int(re.search(r'nks\s*=\s*(\d+)', header).group(1))

        q_points    = np.zeros((nq, 3))
        frequencies = np.zeros((nq, nbnd))

        for n in range(nq):
            q_points[n] = [float(x) for x in f.readline().split()[:3]]
            freqs: list[float] = []
            while len(freqs) < nbnd:
                freqs.extend(float(x) for x in f.readline().split())
            frequencies[n] = freqs[:nbnd]

    return q_points, frequencies


def read_freq_gp(path) -> tuple[np.ndarray, np.ndarray]:
    """Read a matdyn.x ``.freq.gp`` file (gnuplot-ready dispersion).

    Each line contains ``path_length  freq1  freq2  …  freqN``.

    Returns
    -------
    path_lengths : ndarray (nq,)
    frequencies  : ndarray (nq, nbnd)
    """
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1:]


def read_dos(path) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """Read a matdyn.x ``.dos`` file.

    Format::

        # Frequency[cm^-1] DOS PDOS    ← optional header (QE 7.x)
        E   total_DOS   PDOS_atom1   PDOS_atom2  …

    Lines starting with ``#`` are skipped.

    Returns
    -------
    energies : ndarray (ndos,)          frequency axis in cm⁻¹
    dos      : ndarray (ndos,)          total DOS in states/cm⁻¹
    pdos     : ndarray (nat, ndos) or None   per-atom PDOS (None if absent)
    """
    data = np.loadtxt(path, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)

    energies = data[:, 0]
    dos      = data[:, 1]
    pdos     = data[:, 2:].T if data.shape[1] > 2 else None  # (nat, ndos)
    return energies, dos, pdos


def read_modes(path) -> list[dict]:
    """Read a matdyn.x ``.modes`` file.

    Each q-point block is enclosed between ``***…***`` fence lines and contains
    one ``freq (N) = X [THz] = Y [cm-1]`` line per branch followed by *nat*
    eigenvector lines of the form ``( re im  re im  re im )``.

    Returns
    -------
    list of dicts, one per q-point, with keys:

    q        : ndarray (3,)                    q-vector in reciprocal coords
    freq_cm  : ndarray (nbnd,)                 frequencies in cm⁻¹
    freq_thz : ndarray (nbnd,)                 frequencies in THz
    eigvec   : ndarray (nbnd, nat, 3) complex  polarisation vectors
               eigvec[branch, atom, xyz] = displacement amplitude
    """
    q_re    = re.compile(r'q\s*=\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)')
    freq_re = re.compile(
        r'freq\s*\(\s*\d+\s*\)\s*=\s*([-\d.]+)\s*\[THz\]\s*=\s*([-\d.]+)\s*\[cm-1\]'
    )
    ev_re   = re.compile(
        r'\(\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\)'
    )
    fence_re = re.compile(r'^\s*\*{60,}\s*$')

    modes: list[dict] = []
    current: dict | None = None
    cur_evs: list | None = None  # eigenvector rows for the current branch

    with open(path) as f:
        for line in f:
            # ── new q-point ────────────────────────────────────────────────
            q_m = q_re.search(line)
            if q_m:
                current = {
                    'q':        np.array([float(q_m.group(i)) for i in (1, 2, 3)]),
                    'freq_cm':  [],
                    'freq_thz': [],
                    'eigvec':   [],
                }
                cur_evs = None
                continue

            if current is None:
                continue

            # ── new branch ─────────────────────────────────────────────────
            freq_m = freq_re.search(line)
            if freq_m:
                if cur_evs:                          # save previous branch
                    current['eigvec'].append(np.array(cur_evs))
                current['freq_thz'].append(float(freq_m.group(1)))
                current['freq_cm'].append(float(freq_m.group(2)))
                cur_evs = []
                continue

            # ── eigenvector line for one atom ──────────────────────────────
            ev_m = ev_re.search(line)
            if ev_m is not None and cur_evs is not None:
                r = [float(ev_m.group(i)) for i in range(1, 7)]
                cur_evs.append([
                    complex(r[0], r[1]),
                    complex(r[2], r[3]),
                    complex(r[4], r[5]),
                ])
                continue

            # ── fence line (opening or closing) ────────────────────────────
            if fence_re.match(line):
                if cur_evs:                          # save last branch
                    current['eigvec'].append(np.array(cur_evs))
                    cur_evs = None
                if current['freq_cm']:               # closing fence → store block
                    modes.append({
                        'q':        current['q'],
                        'freq_cm':  np.array(current['freq_cm']),
                        'freq_thz': np.array(current['freq_thz']),
                        'eigvec':   np.array(current['eigvec'], dtype=complex),
                    })
                    current = None
                # else: opening fence — keep going

    return modes


# ============================================================================
# MatdynResults
# ============================================================================

class MatdynResults:
    """Read and analyse output files produced by matdyn.x.

    Construct from the run directory and file names::

        res = MatdynResults.from_input(run_dir, inp)  # recommended
        res = MatdynResults(run_dir, flfrq='si.freq', fldos='si.dos')

    All data is loaded **lazily** on first access.

    Dispersion attributes
    ---------------------
    q_points    : ndarray (nq, 3)       q-vectors in reciprocal coords
    frequencies : ndarray (nq, nbnd)    frequencies in cm⁻¹ (negative = imaginary)
    path        : ndarray (nq,)         cumulative path length (from .freq.gp or computed)

    DOS attributes
    --------------
    energies : ndarray (ndos,)          frequency axis in cm⁻¹
    dos      : ndarray (ndos,)          total DOS in states/cm⁻¹
    pdos     : ndarray (nat, ndos)|None per-atom PDOS (None if not in file)

    Mode attributes
    ---------------
    modes : list of dicts   (see read_modes() for dict keys)
    """

    def __init__(self,
                 run_dir,
                 flfrq: str = 'matdyn.freq',
                 fldos: str = 'matdyn.dos',
                 flvec: str = 'matdyn.modes'):
        self._dir  = Path(run_dir)
        self._flfrq = flfrq
        self._fldos = fldos
        self._flvec = flvec

        self._q_points:    np.ndarray | None = None
        self._frequencies: np.ndarray | None = None
        self._path:        np.ndarray | None = None
        self._energies:    np.ndarray | None = None
        self._dos:         np.ndarray | None = None
        self._pdos:        np.ndarray | None = None
        self._modes:       list | None       = None

    @classmethod
    def from_input(cls, run_dir, inp) -> 'MatdynResults':
        """Construct from a :class:`MatdynInput` object.

        File names are read from ``inp._params``; missing keys fall back
        to matdyn.x defaults.
        """
        p = inp._params
        return cls(
            run_dir=run_dir,
            flfrq=Path(p.get('flfrq', 'matdyn.freq')).name,
            fldos=Path(p.get('fldos',  'matdyn.dos')).name,
            flvec=Path(p.get('flvec',  'matdyn.modes')).name,
        )

    # =========================================================================
    # Internal loaders
    # =========================================================================

    def _load_dispersion(self) -> None:
        freq_path = self._dir / self._flfrq
        if not freq_path.exists():
            raise FileNotFoundError(f'Frequency file not found: {freq_path}')
        self._q_points, self._frequencies = read_freq(freq_path)

        gp_path = self._dir / (self._flfrq + '.gp')
        if gp_path.exists():
            self._path = read_freq_gp(gp_path)[0]
        else:
            diffs = np.diff(self._q_points, axis=0)
            self._path = np.concatenate(
                [[0.0], np.cumsum(np.linalg.norm(diffs, axis=1))]
            )

    def _load_dos(self) -> None:
        dos_path = self._dir / self._fldos
        if not dos_path.exists():
            raise FileNotFoundError(f'DOS file not found: {dos_path}')
        self._energies, self._dos, self._pdos = read_dos(dos_path)

    def _load_modes(self) -> None:
        modes_path = self._dir / self._flvec
        if not modes_path.exists():
            raise FileNotFoundError(f'Modes file not found: {modes_path}')
        self._modes = read_modes(modes_path)

    # =========================================================================
    # Lazy properties
    # =========================================================================

    @property
    def q_points(self) -> np.ndarray:
        if self._q_points is None:
            self._load_dispersion()
        return self._q_points

    @property
    def frequencies(self) -> np.ndarray:
        if self._frequencies is None:
            self._load_dispersion()
        return self._frequencies

    @property
    def path(self) -> np.ndarray:
        if self._path is None:
            self._load_dispersion()
        return self._path

    @property
    def energies(self) -> np.ndarray:
        if self._energies is None:
            self._load_dos()
        return self._energies

    @property
    def dos(self) -> np.ndarray:
        if self._dos is None:
            self._load_dos()
        return self._dos

    @property
    def pdos(self) -> np.ndarray | None:
        if self._dos is None:
            self._load_dos()
        return self._pdos

    @property
    def modes(self) -> list[dict]:
        if self._modes is None:
            self._load_modes()
        return self._modes

    # =========================================================================
    # Inspection
    # =========================================================================

    def dos_data(self) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
        """Return the DOS data as a plain tuple, loading the file if needed.

        Returns
        -------
        energies : ndarray (ndos,)          frequency axis in cm⁻¹
        dos      : ndarray (ndos,)          total DOS in states/cm⁻¹
        pdos     : ndarray (nat, ndos)|None per-atom PDOS, or None if absent
        """
        return self.energies, self.dos, self.pdos

    def dispersion_data(self) -> tuple[np.ndarray, np.ndarray]:
        """Return the dispersion data as a plain tuple, loading the file if needed.

        Returns
        -------
        path        : ndarray (nq,)        cumulative path length along the BZ path
        frequencies : ndarray (nq, nbnd)   phonon frequencies in cm⁻¹
        """
        return self.path, self.frequencies

    def frequencies_at(self, n: int) -> np.ndarray:
        """Return phonon frequencies at q-point index *n* in cm⁻¹.

        Parameters
        ----------
        n : int   index into the q-point list (0-based)

        Returns
        -------
        ndarray (nbnd,)
        """
        return self.frequencies[n]

    def mode_at(self, n: int, band: int) -> np.ndarray:
        """Return the complex polarisation vector for branch *band* at q-point *n*.

        Parameters
        ----------
        n    : int   q-point index (0-based)
        band : int   branch index (0-based)

        Returns
        -------
        ndarray (nat, 3) complex
            Cartesian displacement amplitudes for each atom.
        """
        return self.modes[n]['eigvec'][band]

    # =========================================================================
    # Plotting
    # =========================================================================

    def plot_dispersion(self,
                        ax=None,
                        xticks: list[tuple] | None = None,
                        ylabel: str = r'Frequency (cm$^{-1}$)',
                        color: str = 'steelblue',
                        lw: float = 1.0,
                        **kw) -> 'Axes':
        """Plot the phonon dispersion.

        Parameters
        ----------
        ax : matplotlib Axes | None
            If None a new figure is created.
        xticks : list of (path_coord, label) | None
            High-symmetry point positions.  Each ``path_coord`` is a value
            on the x-axis (same units as ``self.path``).  Vertical lines are
            drawn and the x-axis is labelled accordingly.
        ylabel : str
        color, lw, **kw : passed to ``ax.plot`` for every branch.

        Returns
        -------
        ax : matplotlib Axes
        """
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots()

        x    = self.path
        freq = self.frequencies
        for b in range(freq.shape[1]):
            ax.plot(x, freq[:, b], color=color, lw=lw, **kw)

        ax.axhline(0, color='k', lw=0.5, ls='--')
        ax.set_xlim(x[0], x[-1])
        ax.set_ylabel(ylabel)

        if xticks:
            positions = [p for p, _ in xticks]
            labels    = [lbl for _, lbl in xticks]
            ax.set_xticks(positions)
            ax.set_xticklabels(labels)
            for p in positions:
                ax.axvline(p, color='k', lw=0.5)
        else:
            ax.set_xlabel('Wave vector')

        return ax

    def plot_dos(self,
                 ax=None,
                 with_pdos: bool = True,
                 atom_labels: list[str] | None = None,
                 xlabel: str = r'Frequency (cm$^{-1}$)',
                 ylabel: str = r'DOS (states cm)',
                 colors: list | None = None,
                 **kw) -> 'Axes':
        """Plot the phonon DOS and optionally the per-atom PDOS.

        Parameters
        ----------
        ax : matplotlib Axes | None
        with_pdos : bool
            Overlay per-atom PDOS as filled areas (only if PDOS columns
            are present in the ``.dos`` file).
        atom_labels : list of str | None
            Legend labels for each atom's PDOS (e.g. ``['Si', 'Si']``).
            Defaults to ``'Atom 1'``, ``'Atom 2'``, …
        xlabel, ylabel : str
        colors : list of colour specs for PDOS atoms.
        **kw : passed to ``ax.plot`` for the total DOS line.

        Returns
        -------
        ax : matplotlib Axes
        """
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots()

        ax.plot(self.energies, self.dos, color='k', lw=1.2, label='Total', **kw)

        if with_pdos and self.pdos is not None:
            nat   = self.pdos.shape[0]
            clrs  = colors or plt.rcParams['axes.prop_cycle'].by_key()['color']
            for i in range(nat):
                lbl = (atom_labels[i] if atom_labels and i < len(atom_labels)
                       else f'Atom {i + 1}')
                ax.fill_between(self.energies, self.pdos[i],
                                alpha=0.4, color=clrs[i % len(clrs)], label=lbl)
            ax.legend(fontsize='small')

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(self.energies[0], self.energies[-1])
        ax.set_ylim(bottom=0)
        return ax

    # =========================================================================
    # Repr
    # =========================================================================

    def __repr__(self) -> str:
        loaded = []
        if self._frequencies is not None:
            nq, nb = self._frequencies.shape
            loaded.append(f'dispersion({nq}q×{nb}bands)')
        if self._dos is not None:
            loaded.append(f'dos({len(self._energies)}pts)')
        if self._modes is not None:
            loaded.append(f'modes({len(self._modes)}q)')
        status = ', '.join(loaded) if loaded else 'not loaded'
        return (f'MatdynResults(dir={str(self._dir)!r}, '
                f'flfrq={self._flfrq!r}, [{status}])')
