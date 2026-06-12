"""
Quantum ESPRESSO phonon workflow — namelists reference dictionaries.

Codes covered:
    ph.x      — phonon calculations via density-functional perturbation theory (DFPT)
    q2r.x     — Fourier transform of dynamical matrices to real-space force constants
    matdyn.x  — phonon dispersion and DOS by Fourier interpolation of force constants

These three codes form the standard QE phonon workflow:
    pw.x (scf) → ph.x → q2r.x → matdyn.x

Same schema as pw_namelists.py:
    default     : default value (None if required / no default)
    type        : Python type hint string
    unit        : physical unit or \'\' if dimensionless/not applicable
    description : concise human-readable explanation
    valid       : list of valid choices (empty list = free-form)
"""

PH_INPUTPH = {

    # ── I/O and job control ──────────────────────────────────────────────────
    "prefix": {
        "default": "pwscf",
        "type": "str",
        "unit": "",
        "description": "Must match the prefix used in the pw.x ground-state calculation.",
        "valid": [],
    },
    "outdir": {
        "default": "./",
        "type": "str",
        "unit": "",
        "description": "Directory containing pw.x save files (tmp_dir in older versions).",
        "valid": [],
    },
    "fildyn": {
        "default": "matdyn",
        "type": "str",
        "unit": "",
        "description": (
            "Root name for dynamical matrix output files. "
            "For ldisp=.true., files are named <fildyn>1, <fildyn>2, … for each q-point."
        ),
        "valid": [],
    },
    "fildrho": {
        "default": "drho",
        "type": "str",
        "unit": "",
        "description": "Output file name for the first-order change in charge density δρ.",
        "valid": [],
    },
    "fildvscf": {
        "default": "dvscf",
        "type": "str",
        "unit": "",
        "description": (
            "Root name for the first-order change in the SCF potential δV_scf. "
            "Required by electron-phonon and dvscf_star workflows."
        ),
        "valid": [],
    },
    "verbosity": {
        "default": "default",
        "type": "str",
        "unit": "",
        "description": "Output verbosity: 'default' prints standard info; 'high' adds debug output.",
        "valid": ["default", "high"],
    },
    "max_seconds": {
        "default": 1.0e7,
        "type": "float",
        "unit": "s",
        "description": "Wall-clock time limit. ph.x saves and exits cleanly before this threshold.",
        "valid": [],
    },
    "recover": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Restart from a previous incomplete run. "
            "ph.x reads existing partial results from the _ph0 directory."
        ),
        "valid": [],
    },
    "reduce_io": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Keep wavefunctions in memory to reduce disk I/O (costs more RAM).",
        "valid": [],
    },
    "low_directory_check": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Use a faster but less rigorous check for existing output directories.",
        "valid": [],
    },
    "xml_file": {
        "default": "",
        "type": "str",
        "unit": "",
        "description": "Override the XML data file read from the pw.x save directory.",
        "valid": [],
    },
    "lqdir": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Store each q-point's partial results in a separate sub-directory "
            "<outdir>/_ph0/<prefix>.q_N/. Useful for parallel workflows."
        ),
        "valid": [],
    },

    # ── q-point specification ────────────────────────────────────────────────
    "ldisp": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Compute phonons on a uniform Monkhorst-Pack q-grid defined by nq1×nq2×nq3. "
            "Generates one dynamical matrix file per irreducible q-point."
        ),
        "valid": [],
    },
    "nq1": {
        "default": 1,
        "type": "int",
        "unit": "",
        "description": "q-mesh dimension along the first reciprocal lattice vector (ldisp=.true.).",
        "valid": [],
    },
    "nq2": {
        "default": 1,
        "type": "int",
        "unit": "",
        "description": "q-mesh dimension along the second reciprocal lattice vector.",
        "valid": [],
    },
    "nq3": {
        "default": 1,
        "type": "int",
        "unit": "",
        "description": "q-mesh dimension along the third reciprocal lattice vector.",
        "valid": [],
    },
    "qplot": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Compute phonons at an explicit list of q-points supplied after the namelist. "
            "Format: nqs \\n q1x q1y q1z nq1_points \\n …"
        ),
        "valid": [],
    },
    "q_in_band_form": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "When qplot=.true., interpret q-point list as a band path: "
            "each line gives a high-symmetry endpoint and the number of points to the next."
        ),
        "valid": [],
    },
    "q_in_cryst_coord": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Express q-points in crystallographic (reduced) reciprocal coordinates.",
        "valid": [],
    },
    "start_q": {
        "default": 1,
        "type": "int",
        "unit": "",
        "description": (
            "Index of the first q-point to compute in the ldisp grid. "
            "Used together with last_q to split a large run across jobs."
        ),
        "valid": [],
    },
    "last_q": {
        "default": -1000,
        "type": "int",
        "unit": "",
        "description": "Index of the last q-point to compute (-1000 = all).",
        "valid": [],
    },
    "lshift_q": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Shift the q-grid so that Γ is not included. "
            "Useful when a separate Γ calculation is done with a coarser grid."
        ),
        "valid": [],
    },

    # ── SCF convergence ──────────────────────────────────────────────────────
    "tr2_ph": {
        "default": 1.0e-12,
        "type": "float",
        "unit": "",
        "description": (
            "Convergence threshold for the self-consistent phonon response. "
            "The SCF loop stops when the RMS change in δV_scf falls below tr2_ph."
        ),
        "valid": [],
    },
    "niter_ph": {
        "default": 100,
        "type": "int",
        "unit": "",
        "description": "Maximum number of phonon SCF iterations per irreducible representation.",
        "valid": [],
    },
    "nmix_ph": {
        "default": 4,
        "type": "int",
        "unit": "",
        "description": "Number of previous iterations kept in the Broyden/Pulay mixing history.",
        "valid": [],
    },
    "alpha_mix": {
        "default": 0.7,
        "type": "float",
        "unit": "",
        "description": "Linear mixing parameter for the phonon SCF potential update (0 < α ≤ 1).",
        "valid": [],
    },
    "wpot_order": {
        "default": 2,
        "type": "int",
        "unit": "",
        "description": (
            "Order of the polynomial used to extrapolate the SCF potential between "
            "phonon perturbations (improves convergence for systems with many modes)."
        ),
        "valid": [0, 1, 2],
    },

    # ── What to compute ──────────────────────────────────────────────────────
    "trans": {
        "default": True,
        "type": "bool",
        "unit": "",
        "description": "Compute the dynamical matrix and phonon frequencies/eigenvectors.",
        "valid": [],
    },
    "epsil": {
        "default": None,
        "type": "bool",
        "unit": "",
        "description": (
            "Compute the macroscopic dielectric tensor ε∞ (electronic contribution). "
            "Default: same as trans. Only at q=Γ."
        ),
        "valid": [],
    },
    "zeu": {
        "default": None,
        "type": "bool",
        "unit": "",
        "description": (
            "Compute Born effective charges Z* via the electric-field perturbation. "
            "Default: same as trans. Only at q=Γ."
        ),
        "valid": [],
    },
    "zue": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Compute Born effective charges from the unsymmetrised mixed second derivative "
            "(uses atomic displacements at q=Γ instead of the electric-field perturbation). "
            "Results should agree with zeu=.true. as a consistency check."
        ),
        "valid": [],
    },
    "lraman": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Compute non-resonant Raman tensors χ^{(2)}_{αβγ} via third-order DFPT. "
            "Requires epsil=.true. and a Γ-point calculation."
        ),
        "valid": [],
    },
    "elop": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Compute the electro-optic tensor (requires lraman=.true.).",
        "valid": [],
    },
    "fpol": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Compute the frequency-dependent dielectric function ε(ω) "
            "at a set of frequencies specified in the input after the namelist."
        ),
        "valid": [],
    },
    "lnscf": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Perform an internal nscf step before the phonon calculation "
            "(useful when the pw.x save has only SCF k-points and more are needed)."
        ),
        "valid": [],
    },
    "ldiag": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Diagonalise the dynamical matrix and print frequencies at the end of the run.",
        "valid": [],
    },
    "nogg": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Skip the non-analytic (Γ-point) correction to the dynamical matrix. "
            "Use when epsil/zeu are not computed or for non-polar materials."
        ),
        "valid": [],
    },
    "asr": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Enforce the acoustic sum rule on the dynamical matrix before printing frequencies. "
            "Equivalent to zasr='simple' in matdyn.x."
        ),
        "valid": [],
    },
    "lrpa": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Use the random-phase approximation (no exchange-correlation kernel in the response).",
        "valid": [],
    },
    "lnoloc": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Neglect local-field effects in the dielectric response.",
        "valid": [],
    },
    "only_init": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Initialise all internal data structures and write them to disk, then stop. "
            "Used to prepare a parallel run where individual irreps are computed separately."
        ),
        "valid": [],
    },
    "search_sym": {
        "default": True,
        "type": "bool",
        "unit": "",
        "description": (
            "Search for the small-group symmetry at each q-point to reduce the "
            "number of irreducible representations to compute."
        ),
        "valid": [],
    },

    # ── Selecting representations / atoms ────────────────────────────────────
    "start_irr": {
        "default": 1,
        "type": "int",
        "unit": "",
        "description": (
            "Index of the first irreducible representation to compute at each q-point. "
            "Combine with last_irr to parallelise over representations across jobs."
        ),
        "valid": [],
    },
    "last_irr": {
        "default": -1000,
        "type": "int",
        "unit": "",
        "description": "Index of the last irreducible representation to compute (-1000 = all).",
        "valid": [],
    },
    "nat_todo": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": (
            "Number of atoms for which perturbations are computed (0 = all). "
            "The indices of the selected atoms follow the namelist."
        ),
        "valid": [],
    },
    "modenum": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "Compute only a single phonon mode (0 = all modes).",
        "valid": [],
    },

    # ── Electron-phonon coupling ─────────────────────────────────────────────
    "electron_phonon": {
        "default": "",
        "type": "str",
        "unit": "",
        "description": (
            "Activate electron-phonon coupling and select the method:\n"
            "  ''             = disabled\n"
            "  'simple'       = e-ph matrix elements on the coarse k-mesh\n"
            "  'interpolated' = Wannier-based interpolation (needs EPW)\n"
            "  'lambda_tetra' = λ and ω_log via tetrahedron integration\n"
            "  'gamma_tetra'  = phonon linewidths via tetrahedra\n"
            "  'yambo'        = write e-ph matrix elements for YAMBO\n"
            "  'ahc'          = Allen-Heine-Cardona self-energy (needs ahc_* params)"
        ),
        "valid": ["", "simple", "interpolated", "lambda_tetra", "gamma_tetra", "yambo", "ahc"],
    },
    "el_ph_sigma": {
        "default": 0.02,
        "type": "float",
        "unit": "Ry",
        "description": (
            "Gaussian smearing width for the Fermi-surface delta function "
            "in electron-phonon calculations."
        ),
        "valid": [],
    },
    "el_ph_nsigma": {
        "default": 10,
        "type": "int",
        "unit": "",
        "description": "Number of smearing values sampled between el_ph_sigma and 10×el_ph_sigma.",
        "valid": [],
    },
    "dek": {
        "default": 1.0e-3,
        "type": "float",
        "unit": "2π/a",
        "description": (
            "Finite-difference step in k-space for computing the e-ph coupling "
            "with electron_phonon='simple' or 'interpolated'."
        ),
        "valid": [],
    },
    "nk1": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": (
            "Dense k-mesh for electron-phonon (lambda_tetra / gamma_tetra). "
            "0 = use the k-mesh from the pw.x calculation."
        ),
        "valid": [],
    },
    "nk2": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "Dense k-mesh along b2 for electron-phonon (0 = from pw.x).",
        "valid": [],
    },
    "nk3": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "Dense k-mesh along b3 for electron-phonon (0 = from pw.x).",
        "valid": [],
    },
    "nk1_epa": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": (
            "k-mesh along b1 for the electron-phonon averaged (EPA) approximation. "
            "0 disables the EPA output."
        ),
        "valid": [],
    },
    "nk2_epa": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "k-mesh along b2 for EPA.",
        "valid": [],
    },
    "nk3_epa": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "k-mesh along b3 for EPA.",
        "valid": [],
    },

    # ── AHC (Allen-Heine-Cardona) self-energy ────────────────────────────────
    "ahc_dir": {
        "default": "./ahc_dir",
        "type": "str",
        "unit": "",
        "description": "Output directory for AHC electron-phonon self-energy files.",
        "valid": [],
    },
    "ahc_nbnd": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": (
            "Number of bands included in the AHC self-energy summation. "
            "0 = all bands up to nbnd from pw.x."
        ),
        "valid": [],
    },
    "ahc_nbndskip": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "Number of lowest bands skipped (excluded from the AHC sum).",
        "valid": [],
    },
    "skip_upperfan": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Skip the upper Fan self-energy contribution in the AHC calculation "
            "(includes only the Debye-Waller term and lower Fan)."
        ),
        "valid": [],
    },
    "ldvscf_interpolate": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Interpolate δV_scf from a coarse q-grid to the target q-points "
            "using a Fourier method (requires fildvscf and the dvscf_star data)."
        ),
        "valid": [],
    },
    "read_dns_bare": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Read the bare (non-self-consistent) density-matrix response from disk "
            "instead of recomputing it (speeds up restarts of AHC calculations)."
        ),
        "valid": [],
    },
    "ldoubledelta": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Use a double-delta approximation for the imaginary part of the "
            "AHC self-energy (equivalent to the Fan-Migdal approximation on the Fermi surface)."
        ),
        "valid": [],
    },
    "dw_aw": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": (
            "Include the Debye-Waller term in the spectral function A(ω) "
            "for the AHC self-energy output."
        ),
        "valid": [],
    },

    # ── dvscf_star / drho_star interpolation ─────────────────────────────────
    "dvscf_star": {
        "default": {
            "open": False,
            "dir": "./dvscf_star",
            "ext": "dvscf",
            "basis": "cartesian",
            "pat": False,
        },
        "type": "dict",
        "unit": "",
        "description": (
            "Controls saving/reading of the symmetry-unfolded δV_scf in the star of q. "
            "Sub-keys:\n"
            "  open    (bool)  — activate dvscf_star I/O\n"
            "  dir     (str)   — directory for dvscf_star files\n"
            "  ext     (str)   — filename extension\n"
            "  basis   (str)   — 'cartesian' or 'modes'\n"
            "  pat     (bool)  — write pattern files"
        ),
        "valid": [],
    },
    "drho_star": {
        "default": {
            "open": False,
            "dir": "./drho_star",
            "ext": "drho",
            "basis": "cartesian",
            "pat": False,
        },
        "type": "dict",
        "unit": "",
        "description": (
            "Controls saving/reading of the symmetry-unfolded δρ in the star of q. "
            "Same sub-keys as dvscf_star."
        ),
        "valid": [],
    },
}

PH_NAMELISTS = {"INPUTPH": PH_INPUTPH}


# ============================================================================
# q2r.x
# ============================================================================

Q2R_INPUT = {
    "fildyn": {
        "default": None,
        "type": "str",
        "unit": "",
        "description": "Root name of the dynamical matrix files produced by ph.x (required).",
        "valid": [],
    },
    "flfrc": {
        "default": None,
        "type": "str",
        "unit": "",
        "description": "Output file name for the interatomic force constants (required).",
        "valid": [],
    },
    "zasr": {
        "default": "no",
        "type": "str",
        "unit": "",
        "description": (
            "Type of acoustic sum rule correction applied to the force constants:\n"
            "  'no'        = none\n"
            "  'simple'    = diagonal correction\n"
            "  'crystal'   = impose ASR in crystal coordinates\n"
            "  'one-dim'   = for 1D systems\n"
            "  'zero-dim'  = for 0D (isolated) systems"
        ),
        "valid": ["no", "simple", "crystal", "one-dim", "zero-dim"],
    },
    "loto_2d": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Use 2D LO-TO splitting treatment for polar 2D materials.",
        "valid": [],
    },
}

Q2R_NAMELISTS = {"INPUT": Q2R_INPUT}


# ============================================================================
# matdyn.x
# ============================================================================

MATDYN_INPUT = {
    "flfrc": {
        "default": None,
        "type": "str",
        "unit": "",
        "description": "Input file of interatomic force constants produced by q2r.x (required).",
        "valid": [],
    },
    "asr": {
        "default": "no",
        "type": "str",
        "unit": "",
        "description": "Acoustic sum rule correction applied during interpolation.",
        "valid": ["no", "simple", "crystal", "one-dim", "zero-dim"],
    },
    "dos": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Compute the phonon density of states (requires q-mesh).",
        "valid": [],
    },
    "nk1": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "q-mesh dimension along b1 for DOS calculation.",
        "valid": [],
    },
    "nk2": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "q-mesh dimension along b2.",
        "valid": [],
    },
    "nk3": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "q-mesh dimension along b3.",
        "valid": [],
    },
    "nq_plot": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "Total number of q-points for dispersion plot (0 = read from file).",
        "valid": [],
    },
    "q_in_band_form": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Specify q-path in band form (high-symmetry points + nq).",
        "valid": [],
    },
    "q_in_cryst_coord": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Specify q-points in crystallographic coordinates.",
        "valid": [],
    },
    "flfrq": {
        "default": "matdyn.freq",
        "type": "str",
        "unit": "",
        "description": "Output file for phonon frequencies along the q-path.",
        "valid": [],
    },
    "flvec": {
        "default": "matdyn.modes",
        "type": "str",
        "unit": "",
        "description": "Output file for phonon eigenvectors (polarisation vectors).",
        "valid": [],
    },
    "fldos": {
        "default": "matdyn.dos",
        "type": "str",
        "unit": "",
        "description": "Output file for the phonon density of states.",
        "valid": [],
    },
    "fldyn": {
        "default": "",
        "type": "str",
        "unit": "",
        "description": "Output file for dynamical matrices at interpolated q-points.",
        "valid": [],
    },
    "fleig": {
        "default": "",
        "type": "str",
        "unit": "",
        "description": "Output file for phonon group velocities.",
        "valid": [],
    },
    "fltau": {
        "default": "",
        "type": "str",
        "unit": "",
        "description": "Output file for force constants in real space.",
        "valid": [],
    },
    "ndos": {
        "default": 50,
        "type": "int",
        "unit": "",
        "description": "Number of frequency bins in the phonon DOS.",
        "valid": [],
    },
    "deltaE": {
        "default": 1.0,
        "type": "float",
        "unit": "cm^{-1}",
        "description": "Energy resolution (bin width) of the DOS.",
        "valid": [],
    },
    "degauss": {
        "default": 0.0,
        "type": "float",
        "unit": "cm^{-1}",
        "description": "Gaussian smearing width for the DOS.",
        "valid": [],
    },
    "la2F": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Compute the Eliashberg spectral function a²F(ω).",
        "valid": [],
    },
    "eigen_similarity": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Sort phonon branches using eigenvector similarity.",
        "valid": [],
    },
    "fd": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Use finite-difference force constants (for 2D LO-TO).",
        "valid": [],
    },
    "na_ifc": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Include non-analytic correction to interatomic force constants.",
        "valid": [],
    },
    "nosym": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Disable symmetry reduction of the q-mesh.",
        "valid": [],
    },
    "loto_2d": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Use 2D LO-TO splitting for polar 2D materials.",
        "valid": [],
    },
    "loto_disable": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Disable LO-TO splitting correction entirely.",
        "valid": [],
    },
    "at": {
        "default": None,
        "type": "list[float]",
        "unit": "bohr",
        "description": "Lattice vectors (alternative to reading from flfrc).",
        "valid": [],
    },
    "readtau": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Read atomic positions from input (overrides those in flfrc).",
        "valid": [],
    },
    "amass": {
        "default": None,
        "type": "list[float]",
        "unit": "amu",
        "description": "Atomic masses (override values stored in flfrc).",
        "valid": [],
    },
}

MATDYN_NAMELISTS = {"INPUT": MATDYN_INPUT}


# ============================================================================
# Convenience bundle
# ============================================================================
PH_WORKFLOW_NAMELISTS = {
    "ph.x":     {"namelists": PH_NAMELISTS},
    "q2r.x":    {"namelists": Q2R_NAMELISTS},
    "matdyn.x": {"namelists": MATDYN_NAMELISTS},
}


def describe(namelist: dict, key: str) -> None:
    """Print a formatted description of a single parameter."""
    p = namelist[key]
    print(f"  {key}")
    print(f"    type    : {p.get('type', '')}")
    if p.get("unit"):
        print(f"    unit    : {p['unit']}")
    print(f"    default : {p.get('default', '—')}")
    if p.get("valid"):
        print(f"    valid   : {p['valid']}")
    print(f"    info    : {p['description']}")


def defaults_from(namelist: dict) -> dict:
    """Return {param: default} for every param that has a non-None default."""
    return {k: v["default"] for k, v in namelist.items() if v["default"] is not None}
