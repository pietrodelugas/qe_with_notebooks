"""elastic_tools.py — input builders and analysis helpers for MgO elastic constants.

Public API
----------
Constants
    EV_ANG3_TO_GPA

Input builders  (return a PWInput with calculation='relax')
    build_ortho_input(a0_ang, x, ecutwfc, nk, prefix, pseudo_dir, outdir, pseudos)
    build_mono_input (a0_ang, x, ecutwfc, nk, prefix, pseudo_dir, outdir, pseudos)

Fitting / analysis  (pure functions)
    fit_elastic_quad(x_values, dE_values)          -> coeff (eV)
    extract_c11_c12(coeff_ev, V0_conv_ang3, B_gpa) -> (C11, C12)  GPa
    extract_c44    (coeff_ev, V0_conv_ang3)         -> C44          GPa
"""

import numpy as np
from ase.build import bulk

from pw_input import (
    ControlNamelist, SystemNamelist, ElectronsNamelist, IonsNamelist,
    AtomicSpeciesCard, AtomicPositionsCard, KPointsAutoCard, PWInput,
)

_ANG_TO_BOHR   = 1.0 / 0.529177210903
EV_ANG3_TO_GPA = 160.2176   # 1 eV/Å³ in GPa

# ---------------------------------------------------------------------------
# Input builders
# ---------------------------------------------------------------------------

def build_ortho_input(a0_ang, x, ecutwfc, nk, prefix, pseudo_dir, outdir, pseudos):
    """PWInput for orthorhombic strain at amplitude x (ibrav=8, relax).

    Strain tensor (volume-conserving):
        ε = diag(x,  −x,  x²/(1−x²))

    Resulting celldm for ibrav=8:
        celldm(1) = a₀(1+x)       [bohr]
        celldm(2) = (1−x)/(1+x)
        celldm(3) = 1 / ((1−x²)(1+x))

    Parameters
    ----------
    a0_ang : float  — equilibrium lattice parameter in Å (from BM EOS fit)
    x      : float  — strain amplitude (dimensionless, typically ±0.01 … ±0.06)
    nk     : int    — k-grid size along each direction (nk×nk×nk)
    """
    atoms = bulk('MgO', 'rocksalt', a=a0_ang, cubic=True)

    control   = ControlNamelist(calculation='relax', prefix=prefix,
                                pseudo_dir=str(pseudo_dir), outdir=str(outdir),
                                tprnfor=True, tstress=True)
    system    = SystemNamelist(ibrav=8, nat=8, ntyp=2, ecutwfc=ecutwfc,
                               celldm_1=a0_ang * (1 + x) * _ANG_TO_BOHR,
                               celldm_2=(1 - x) / (1 + x),
                               celldm_3=1 / ((1 - x**2) * (1 + x)))
    electrons = ElectronsNamelist(conv_thr=1.e-9)
    ions      = IonsNamelist()
    species   = AtomicSpeciesCard.from_atoms(atoms, pseudos)
    positions = AtomicPositionsCard.from_atoms(atoms, units='crystal')
    kpoints   = KPointsAutoCard(8, nk1=nk, nk2=nk, nk3=nk)
    return PWInput(control, system, electrons, species, positions, kpoints, ions=ions)


def build_mono_input(a0_ang, x, ecutwfc, nk, prefix, pseudo_dir, outdir, pseudos):
    """PWInput for monoclinic shear strain at amplitude x (ibrav=12, relax).

    Strain tensor (volume-conserving):
        ε = [[0, x/2, 0], [x/2, 0, 0], [0, 0, x²/(4−x²)]]

    Resulting celldm for ibrav=12 (unique axis c):
        celldm(1) = a₀ √(1+x²/4)   [bohr]
        celldm(2) = 1
        celldm(3) = 4 / ((4−x²) √(1+x²/4))
        celldm(4) = cos γ = x / (1+x²/4)

    Parameters
    ----------
    a0_ang : float  — equilibrium lattice parameter in Å
    x      : float  — strain amplitude (dimensionless)
    nk     : int    — k-grid size along each direction (nk×nk×nk)
    """
    atoms = bulk('MgO', 'rocksalt', a=a0_ang, cubic=True)

    f = 1 + x**2 / 4   # shorthand that appears in all four celldm formulas
    control   = ControlNamelist(calculation='relax', prefix=prefix,
                                pseudo_dir=str(pseudo_dir), outdir=str(outdir),
                                tprnfor=True, tstress=True)
    system    = SystemNamelist(ibrav=12, nat=8, ntyp=2, ecutwfc=ecutwfc,
                               celldm_1=a0_ang * np.sqrt(f) * _ANG_TO_BOHR,
                               celldm_2=1.0,
                               celldm_3=4 / ((4 - x**2) * np.sqrt(f)),
                               celldm_4=x / f)
    electrons = ElectronsNamelist(conv_thr=1.e-9)
    ions      = IonsNamelist()
    species   = AtomicSpeciesCard.from_atoms(atoms, pseudos)
    positions = AtomicPositionsCard.from_atoms(atoms, units='crystal')
    kpoints   = KPointsAutoCard(12, nk1=nk, nk2=nk, nk3=nk)
    return PWInput(control, system, electrons, species, positions, kpoints, ions=ions)


# ---------------------------------------------------------------------------
# Fitting and analysis
# ---------------------------------------------------------------------------

def fit_elastic_quad(x_values, dE_values):
    """Fit ΔE(x) = coeff · x² (forced zero intercept, no linear term).

    Uses least-squares projection: coeff = (x²·ΔE) / (x²·x²).

    Parameters
    ----------
    x_values  : array of strain amplitudes
    dE_values : array of ΔE = E(x) − E(0)  [eV/cell]

    Returns
    -------
    coeff : float  [eV]
    """
    x2 = np.asarray(x_values, dtype=float) ** 2
    dE = np.asarray(dE_values, dtype=float)
    return float(np.dot(x2, dE) / np.dot(x2, x2))


def extract_c11_c12(coeff_ev, V0_conv_ang3, B_gpa):
    """Extract C₁₁ and C₁₂ (GPa) from the orthorhombic fit.

    Relations used:
        ΔE = V₀ (C₁₁ − C₁₂) x²   →   C₁₁ − C₁₂ = coeff / V₀
        B  = (C₁₁ + 2 C₁₂) / 3   →   C₁₁ + 2 C₁₂ = 3 B

    Parameters
    ----------
    coeff_ev      : float — fit coefficient from fit_elastic_quad [eV]
    V0_conv_ang3  : float — equilibrium volume of the 8-atom conventional cell [Å³]
                    (= 4 × primitive-cell volume from BM fit)
    B_gpa         : float — bulk modulus from BM EOS [GPa]
    """
    C11_minus_C12 = coeff_ev / V0_conv_ang3 * EV_ANG3_TO_GPA
    C11 = (C11_minus_C12 + 3 * B_gpa) / 3
    C12 = (3 * B_gpa - C11_minus_C12) / 3
    return C11, C12


def extract_c44(coeff_ev, V0_conv_ang3):
    """Extract C₄₄ (GPa) from the monoclinic shear fit.

    Relation used:
        ΔE = ½ V₀ C₄₄ x²   →   C₄₄ = 2 coeff / V₀

    Parameters
    ----------
    coeff_ev     : float — fit coefficient [eV]
    V0_conv_ang3 : float — equilibrium volume of the 8-atom conventional cell [Å³]
    """
    return 2 * coeff_ev / V0_conv_ang3 * EV_ANG3_TO_GPA
