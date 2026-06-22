"""bandstructure_tools.py — helpers for band structure, DOS, and fat-bands workflows.

Public API
----------
Parsers
    parse_fermi_energy(stdout)            -> float  [eV]
    parse_hs_positions(bandsx_stdout)     -> list[float]  (x-coords of high-symmetry points)
    parse_gamma_symmetries(bandsx_stdout) -> list[dict]   (IRREP labels at Γ, requires lsym=.true.)

Band structure
    read_bands_gnu(filepath)              -> (k_coords, bands_ev)
    build_kpath_str(points, npt)          -> str  (K_POINTS crystal_b block)
    print_bands_summary(k_coords, bands_ev, E_F, hs_x, path, n_occ)
    plot_band_structure(k_coords, bands_ev, E_F, hs_x, hs_labels, gap_ev)

Fat bands
    read_projwfc_weights(filproj_path)    -> dict  keys: weights, natomwfc, nkpts, nbnd
"""

import re
import numpy as np


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_fermi_energy(stdout):
    """Return the Fermi / highest-occupied-level energy in eV from pw.x stdout.

    Tries 'highest occupied level' (insulator) first, then 'the Fermi energy is'
    (metal).  Always takes the last match so the function is safe for relax runs
    that print multiple SCF cycles.
    """
    hits = re.findall(r'highest occupied level \(ev\):\s+([-\d.]+)', stdout)
    if hits:
        return float(hits[-1])
    hits = re.findall(r'the Fermi energy is\s+([-\d.]+)', stdout)
    if hits:
        return float(hits[-1])
    raise ValueError('Fermi / highest-occupied-level energy not found in stdout')


def parse_hs_positions(bandsx_stdout):
    """Parse high-symmetry point x-coordinates from bands.x stdout.

    bands.x prints one line per path point:
        high-symmetry point:  kx  ky  kz   x=  X.XXXX

    Returns
    -------
    list of float — cumulative k-path distances for each high-symmetry point,
    in the same units as the x-column of the *.dat.gnu file.
    """
    return [float(x) for x in re.findall(r'high-symmetry point:.*?x(?:\s*=|\s+coordinate)\s+([\d.]+)', bandsx_stdout)]


def parse_gamma_symmetries(bandsx_stdout):
    """Parse band symmetry labels at Γ from bands.x stdout (requires lsym=.true.).

    When lsym=.true., bands.x prints at each high-symmetry k-point a block like:

        Band symmetry, Oh (m-3m) point group:

        e(  1 -  1) =  -73.82 eV   -->  G_1+  (  1)
        e(  3 -  5) =   -4.95 eV   -->  G_15  (  3)
        ...

    This function returns the block for the first Γ = (0,0,0) point found.

    Returns
    -------
    list of dict with keys:
        bands      : (n_start, n_end) — 1-based band indices
        energy_ev  : float
        label      : str  — IRREP label (group prefix stripped)
        deg        : int  — degeneracy
    """
    results = []
    in_gamma = False
    in_symm  = False

    for line in bandsx_stdout.splitlines():
        if re.search(r'xk=\(\s*0\.0+\s*,\s*0\.0+\s*,\s*0\.0+', line):
            in_gamma = True
            in_symm  = False
            results  = []
            continue
        if not in_gamma:
            continue
        if 'Band symmetry' in line:
            in_symm = True
            continue
        if not in_symm:
            continue
        # format: e(  1 -  1) =    -64.55  eV     1   --> A_1g G_1   G_1+
        m = re.match(
            r'\s*e\(\s*(\d+)\s*-\s*(\d+)\s*\)\s*=\s*([-\d.]+)\s*eV\s+(\d+)\s*-->\s*(.*)', line)
        if m:
            n1, n2  = int(m.group(1)), int(m.group(2))
            e_ev    = float(m.group(3))
            deg     = int(m.group(4))
            label   = m.group(5).strip()
            results.append({'bands': (n1, n2), 'energy_ev': e_ev, 'label': label, 'deg': deg})
        elif results and not line.strip():
            break   # blank line ends the symmetry block

    return results


# ---------------------------------------------------------------------------
# Band structure
# ---------------------------------------------------------------------------

def print_bands_summary(k_coords, bands_ev, E_F, hs_x, path, n_occ):
    """Print a concise summary of the computed band structure."""
    gap_ev = bands_ev[n_occ].min() - bands_ev[n_occ - 1].max()
    print(f'k-points on path : {len(k_coords)}')
    print(f'bands computed   : {len(bands_ev)}')
    print(f'E_F (highest occ): {E_F:.4f} eV')
    print(f'Valence band max : {bands_ev[n_occ - 1].max() - E_F:+.4f} eV  (= E_F by definition)')
    print(f'Cond. band min   : {bands_ev[n_occ].min() - E_F:+.4f} eV')
    print(f'DFT-PBE gap      : {gap_ev:.4f} eV  (experiment: ~7.8 eV)')
    print()
    print('High-symmetry points and arc-lengths:')
    labels = ['G' if p[0] == 'G' else p[0] for p in path]
    for label, x in zip(labels, hs_x):
        print(f'  {label:2s}  x = {x:.4f}')


def read_bands_gnu(filepath):
    """Parse a bands.x *.dat.gnu file.

    The file contains one block per band, blocks separated by blank lines.
    Each line within a block is:  k_linear_coord   energy_eV

    Returns
    -------
    k_coords : ndarray, shape (nk,)      — linear k-path coordinate
    bands_ev : ndarray, shape (nbnd, nk) — band energies in eV
    """
    blocks = []
    current = []
    with open(filepath) as fh:
        for line in fh:
            stripped = line.strip()
            if stripped:
                current.append(list(map(float, stripped.split())))
            elif current:
                blocks.append(current)
                current = []
    if current:
        blocks.append(current)

    k_coords = np.array([row[0] for row in blocks[0]])
    bands_ev  = np.array([[row[1] for row in blk] for blk in blocks])
    return k_coords, bands_ev


def build_kpath_str(points, npt):
    """Build a K_POINTS crystal_b input block for pw.x.

    Parameters
    ----------
    points : list of (label, kx, ky, kz)
        High-symmetry points in crystal coordinates.  The last point receives
        npt=1 automatically (it is only an endpoint, not a segment start).
    npt : int or list of int
        Number of k-points per segment.  Pass a single int to use the same
        count for every segment, or a list of len(points)-1 ints to set each
        segment independently (useful for proportional sampling where longer
        segments get more points).

    Returns
    -------
    str — full K_POINTS crystal_b block, ready to paste into a pw.x input.

    Example
    -------
    >>> FCC_PATH = [
    ...     ('L', 0.5,   0.5,   0.5  ),
    ...     ('G', 0.0,   0.0,   0.0  ),
    ...     ('X', 0.5,   0.0,   0.5  ),
    ...     ('W', 0.5,   0.25,  0.75 ),
    ...     ('K', 0.375, 0.375, 0.75 ),
    ...     ('G', 0.0,   0.0,   0.0  ),
    ... ]
    >>> print(build_kpath_str(FCC_PATH, npt=30))
    >>> print(build_kpath_str(FCC_PATH, npt=[50, 60, 30, 20, 64]))
    """
    n = len(points)
    if isinstance(npt, int):
        counts = [npt] * (n - 1) + [1]
    else:
        counts = list(npt) + [1]
    lines = ['K_POINTS crystal_b', str(n)]
    for i, (label, kx, ky, kz) in enumerate(points):
        lines.append(f'  {kx:.6f}  {ky:.6f}  {kz:.6f}  {counts[i]:3d}  ! {label}')
    return '\n'.join(lines)


def plot_band_structure(k_coords, bands_ev, E_F, hs_x, hs_labels, gap_ev):
    """Two-panel band structure plot: full range (semicores visible) and valence/gap zoom."""
    import matplotlib.pyplot as plt

    fig, (ax_full, ax_gap) = plt.subplots(1, 2, figsize=(11, 5))

    for ax in (ax_full, ax_gap):
        for band in bands_ev:
            ax.plot(k_coords, band - E_F, color='steelblue', lw=0.9)
        ax.axhline(0, color='gray', lw=0.7, ls='--')
        for x in hs_x[1:-1]:
            ax.axvline(x, color='k', lw=0.6)
        ax.set_xticks(hs_x)
        ax.set_xticklabels(hs_labels)
        ax.set_xlim(k_coords[0], k_coords[-1])
        ax.set_ylabel('$E - E_F$ (eV)')

    ax_full.set_ylim(-80, 15)
    ax_full.set_title('MgO bands — full range')

    ax_gap.set_ylim(-25, 15)
    ax_gap.set_title('MgO bands — valence & gap')
    ax_gap.annotate(f'gap = {gap_ev:.2f} eV',
                    xy=(hs_x[0] * 0.05 + hs_x[1] * 0.95, gap_ev / 2),
                    fontsize=9, color='firebrick')

    fig.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Fat bands
# ---------------------------------------------------------------------------

def read_projwfc_weights(filproj_path):
    """Parse a projwfc.x *.projwfc_up file (written by write_proj_file in write_proj.f90).

    File structure
    --------------
    Header (formatted text, written by write_io_header):
        Line 0  : title string
        Line 1  : nr1x nr2x nr3x nr1 nr2 nr3  nat  ntyp    (8 integers)
        Line 2  : ibrav  celldm(1..6)
        Lines 3-5: lattice vectors (only if ibrav == 0)
        Line    : gcutm  dual  ecutwfc  9
        ntyp lines : nt  atm  zv
        nat  lines : na  tau1  tau2  tau3  ityp
        Line    : natomwfc  nkstot  nbnd
        Line    : noncolin  lspinorb

    Projection data (one block per atomic wavefunction):
        nwfc  na  atom_label  orbital_label  n  l  m
        ik  ibnd  |<phi_nwfc|psi_{ibnd,ik}>|^2
        …  (nkstot × nbnd lines)

    Returns
    -------
    dict with keys:
        weights  : ndarray shape (nkpts, nbnd, natomwfc)
                   weights[ik, ibnd, i] = |<phi_i|psi_{ibnd,ik}>|^2
        orbitals : list of dicts, one per atomic wfc
                   keys: nwfc (1-based), na, atom (str), els (str), n, l, m
        natomwfc, nkpts, nbnd : int
    """
    with open(filproj_path) as fh:
        lines = fh.readlines()

    # --- parse header ---
    vals   = lines[1].split()
    nat    = int(vals[6])
    ntyp   = int(vals[7])
    ibrav  = int(lines[2].split()[0])
    extra  = 3 if ibrav == 0 else 0          # optional lattice-vector lines

    # line index of "natomwfc  nkstot  nbnd"
    idx = 3 + extra + 1 + ntyp + nat
    natomwfc, nkstot, nbnd = map(int, lines[idx].split())

    # projection data starts after the noncolin/lspinorb line
    i = idx + 2

    weights  = np.zeros((nkstot, nbnd, natomwfc))
    orbitals = []

    for nwfc_idx in range(natomwfc):
        parts = lines[i].split()          # nwfc na atom els n l m
        orbitals.append({
            'nwfc': int(parts[0]),
            'na':   int(parts[1]),
            'atom': parts[2].strip(),
            'els':  parts[3].strip(),
            'n':    int(parts[4]),
            'l':    int(parts[5]),
            'm':    int(parts[6]),
        })
        i += 1
        for ik in range(nkstot):
            for ibnd in range(nbnd):
                tok = lines[i].split()
                weights[ik, ibnd, nwfc_idx] = float(tok[2])
                i += 1

    return {
        'weights':  weights,
        'orbitals': orbitals,
        'natomwfc': natomwfc,
        'nkpts':    nkstot,
        'nbnd':     nbnd,
    }
