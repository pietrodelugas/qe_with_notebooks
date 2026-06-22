"""
Quantum ESPRESSO pw.x — namelists and cards reference dictionaries.

Each entry is a dict with:
    default   : default value (None if required / no default)
    type      : Python type hint string
    unit      : physical unit or '' if dimensionless
    description: concise human-readable explanation
    valid     : list of valid choices (empty list = free-form value)

Usage example
-------------
>>> from pw_namelists import CONTROL, SYSTEM, ELECTRONS, IONS, CELL
>>> print(CONTROL['calculation']['description'])
>>> # Build a minimal pw.x input dict
>>> inp = {k: v['default'] for k, v in CONTROL.items() if v['default'] is not None}
"""

# ---------------------------------------------------------------------------
# &CONTROL
# ---------------------------------------------------------------------------
CONTROL = {
    "calculation": {
        "default": "scf",
        "type": "str",
        "unit": "",
        "description": "Type of calculation to perform.",
        "valid": ["scf", "nscf", "bands", "relax", "md", "vc-relax", "vc-md"],
    },
    "restart_mode": {
        "default": "from_scratch",
        "type": "str",
        "unit": "",
        "description": "Start a fresh calculation or restart from a checkpoint.",
        "valid": ["from_scratch", "restart"],
    },
    "prefix": {
        "default": "pwscf",
        "type": "str",
        "unit": "",
        "description": "Prefix for all output file names.",
        "valid": [],
    },
    "outdir": {
        "default": "./",
        "type": "str",
        "unit": "",
        "description": "Directory for temporary and output files.",
        "valid": [],
    },
    "pseudo_dir": {
        "default": "./",
        "type": "str",
        "unit": "",
        "description": "Directory where pseudopotential files are stored.",
        "valid": [],
    },
    "verbosity": {
        "default": "low",
        "type": "str",
        "unit": "",
        "description": "Level of output verbosity.",
        "valid": ["low", "medium", "high", "debug"],
    },
    "tprnfor": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Print forces on atoms at the end of each SCF step.",
        "valid": [],
    },
    "tstress": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Print stress tensor at the end of each SCF step.",
        "valid": [],
    },
    "nstep": {
        "default": 50,
        "type": "int",
        "unit": "",
        "description": "Number of ionic/molecular dynamics steps.",
        "valid": [],
    },
    "iprint": {
        "default": None,
        "type": "int",
        "unit": "",
        "description": "Band energies are written every iprint steps.",
        "valid": [],
    },
    "forc_conv_thr": {
        "default": 1e-3,
        "type": "float",
        "unit": "Ry/bohr",
        "description": "Convergence threshold for forces in relax/vc-relax.",
        "valid": [],
    },
    "etot_conv_thr": {
        "default": 1e-4,
        "type": "float",
        "unit": "Ry",
        "description": "Convergence threshold for total energy in relax.",
        "valid": [],
    },
    "disk_io": {
        "default": "low",
        "type": "str",
        "unit": "",
        "description": "Amount of data written to disk during the run.",
        "valid": ["none", "nowf", "low", "medium", "high"],
    },
    "wf_collect": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Collect wavefunctions into a single file at the end.",
        "valid": [],
    },
    "max_seconds": {
        "default": 1e7,
        "type": "float",
        "unit": "s",
        "description": "Maximum allowed wall-clock time before a clean stop.",
        "valid": [],
    },
    "dt": {
        "default": 20.0,
        "type": "float",
        "unit": "Ry·a.u.^{-1}",
        "description": "MD time step (in Ry atomic units).",
        "valid": [],
    },
    "lkpoint_dir": {
        "default": True,
        "type": "bool",
        "unit": "",
        "description": "Store k-point wavefunctions in separate sub-directories.",
        "valid": [],
    },
    "lforcet": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "When restarting, read and use old atomic velocities (MD).",
        "valid": [],
    },
}

# ---------------------------------------------------------------------------
# &SYSTEM
# ---------------------------------------------------------------------------
SYSTEM = {
    "ibrav": {
        "default": None,
        "type": "int",
        "unit": "",
        "description": (
            "Bravais-lattice index. 0 = free cell given in CELL_PARAMETERS. "
            "1–14 = standard lattice types."
        ),
        "valid": list(range(0, 15)),
    },
    "celldm": {
        "default": None,
        "type": "list[float]",
        "unit": "bohr",
        "description": (
            "Crystallographic constants: celldm(1) = a, celldm(2) = b/a, etc. "
            "Required unless A, B, C, cosBC … are specified."
        ),
        "valid": [],
    },
    "A": {
        "default": None,
        "type": "float",
        "unit": "Angstrom",
        "description": "Lattice parameter a (alternative to celldm).",
        "valid": [],
    },
    "B": {
        "default": None,
        "type": "float",
        "unit": "Angstrom",
        "description": "Lattice parameter b.",
        "valid": [],
    },
    "C": {
        "default": None,
        "type": "float",
        "unit": "Angstrom",
        "description": "Lattice parameter c.",
        "valid": [],
    },
    "cosBC": {
        "default": None,
        "type": "float",
        "unit": "",
        "description": "cos(alpha): cosine of the angle between b and c.",
        "valid": [],
    },
    "cosAC": {
        "default": None,
        "type": "float",
        "unit": "",
        "description": "cos(beta): cosine of the angle between a and c.",
        "valid": [],
    },
    "cosAB": {
        "default": None,
        "type": "float",
        "unit": "",
        "description": "cos(gamma): cosine of the angle between a and b.",
        "valid": [],
    },
    "nat": {
        "default": None,
        "type": "int",
        "unit": "",
        "description": "Number of atoms in the unit cell (required).",
        "valid": [],
    },
    "ntyp": {
        "default": None,
        "type": "int",
        "unit": "",
        "description": "Number of distinct atomic species (required).",
        "valid": [],
    },
    "nbnd": {
        "default": None,
        "type": "int",
        "unit": "",
        "description": "Number of electronic bands. Default: max(nelec/2+4, 1.2*nelec/2).",
        "valid": [],
    },
    "ecutwfc": {
        "default": None,
        "type": "float",
        "unit": "Ry",
        "description": "Kinetic-energy cutoff for wavefunctions (required).",
        "valid": [],
    },
    "ecutrho": {
        "default": None,
        "type": "float",
        "unit": "Ry",
        "description": "Kinetic-energy cutoff for charge density. Default: 4×ecutwfc (NC) or 12×ecutwfc (US/PAW).",
        "valid": [],
    },
    "occupations": {
        "default": "fixed",
        "type": "str",
        "unit": "",
        "description": "Occupation function for electronic states.",
        "valid": ["smearing", "tetrahedra", "tetrahedra_opt", "tetrahedra_lin", "fixed", "from_input"],
    },
    "smearing": {
        "default": "gaussian",
        "type": "str",
        "unit": "",
        "description": "Smearing type used with occupations='smearing'.",
        "valid": ["gaussian", "methfessel-paxton", "marzari-vanderbilt", "fermi-dirac"],
    },
    "degauss": {
        "default": 0.0,
        "type": "float",
        "unit": "Ry",
        "description": "Smearing width (Gaussian broadening parameter).",
        "valid": [],
    },
    "nspin": {
        "default": 1,
        "type": "int",
        "unit": "",
        "description": "Spin polarisation: 1=unpolarised, 2=collinear spin, 4=non-collinear.",
        "valid": [1, 2, 4],
    },
    "noncolin": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Enable non-collinear magnetism calculation.",
        "valid": [],
    },
    "lspinorb": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Include spin-orbit coupling (requires noncolin=.true.).",
        "valid": [],
    },
    "starting_magnetization": {
        "default": None,
        "type": "list[float]",
        "unit": "",
        "description": "Starting magnetisation for each species, between -1 and 1.",
        "valid": [],
    },
    "tot_charge": {
        "default": 0.0,
        "type": "float",
        "unit": "e",
        "description": "Total charge of the system (positive = electrons removed).",
        "valid": [],
    },
    "tot_magnetization": {
        "default": -1.0,
        "type": "float",
        "unit": "μ_B",
        "description": "Constrained total magnetisation (-1 = unconstrained).",
        "valid": [],
    },
    "input_dft": {
        "default": "",
        "type": "str",
        "unit": "",
        "description": "Override the DFT functional read from pseudopotentials.",
        "valid": [],
    },
    "exx_fraction": {
        "default": None,
        "type": "float",
        "unit": "",
        "description": "Fraction of exact exchange for hybrid functionals.",
        "valid": [],
    },
    "screening_parameter": {
        "default": 0.106,
        "type": "float",
        "unit": "bohr^{-1}",
        "description": "Screening parameter for range-separated hybrids (e.g. HSE06).",
        "valid": [],
    },
    "vdw_corr": {
        "default": "none",
        "type": "str",
        "unit": "",
        "description": "Type of van der Waals correction to apply.",
        "valid": ["none", "grimme-d2", "grimme-d3", "grimme-d3bj", "ts-vdw", "xdm", "dft-d", "dft-d3"],
    },
    "assume_isolated": {
        "default": "none",
        "type": "str",
        "unit": "",
        "description": "Method for treating isolated systems (removes PBC corrections).",
        "valid": ["none", "makov-payne", "martyna-tuckerman", "esm", "2D"],
    },
    "esm_bc": {
        "default": "pbc",
        "type": "str",
        "unit": "",
        "description": "Boundary condition for effective screening medium (ESM).",
        "valid": ["pbc", "bc1", "bc2", "bc3"],
    },
    "la2F": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Compute the electron-phonon coupling function a^2F for superconductivity.",
        "valid": [],
    },
    "lda_plus_u": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Enable DFT+U (Hubbard U) correction.",
        "valid": [],
    },
    "Hubbard_U": {
        "default": None,
        "type": "list[float]",
        "unit": "eV",
        "description": "Hubbard U parameter for each species with lda_plus_u=.true.",
        "valid": [],
    },
    "Hubbard_J0": {
        "default": None,
        "type": "list[float]",
        "unit": "eV",
        "description": "Hubbard J0 parameter for each species.",
        "valid": [],
    },
    "nosym": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Disable all symmetry operations.",
        "valid": [],
    },
    "noinv": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Disable inversion symmetry.",
        "valid": [],
    },
    "no_t_rev": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Disable time-reversal symmetry.",
        "valid": [],
    },
    "force_symmorphic": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Force use of symmorphic operations only.",
        "valid": [],
    },
    "nr1": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "FFT mesh size along x (0 = automatic).",
        "valid": [],
    },
    "nr2": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "FFT mesh size along y (0 = automatic).",
        "valid": [],
    },
    "nr3": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "FFT mesh size along z (0 = automatic).",
        "valid": [],
    },
    "nr1s": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "Smooth FFT mesh size along x (0 = automatic).",
        "valid": [],
    },
    "nr2s": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "Smooth FFT mesh size along y (0 = automatic).",
        "valid": [],
    },
    "nr3s": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "Smooth FFT mesh size along z (0 = automatic).",
        "valid": [],
    },
}

# ---------------------------------------------------------------------------
# &ELECTRONS
# ---------------------------------------------------------------------------
ELECTRONS = {
    "electron_maxstep": {
        "default": 100,
        "type": "int",
        "unit": "",
        "description": "Maximum number of SCF iterations.",
        "valid": [],
    },
    "conv_thr": {
        "default": 1e-6,
        "type": "float",
        "unit": "Ry",
        "description": "Convergence threshold for the SCF total energy.",
        "valid": [],
    },
    "mixing_mode": {
        "default": "plain",
        "type": "str",
        "unit": "",
        "description": "Charge density mixing scheme for SCF convergence.",
        "valid": ["plain", "TF", "local-TF"],
    },
    "mixing_beta": {
        "default": 0.7,
        "type": "float",
        "unit": "",
        "description": "Mixing parameter (0 < β ≤ 1). Smaller values → more stable but slower.",
        "valid": [],
    },
    "mixing_ndim": {
        "default": 8,
        "type": "int",
        "unit": "",
        "description": "Number of iterations used in Broyden/Pulay mixing.",
        "valid": [],
    },
    "diagonalization": {
        "default": "david",
        "type": "str",
        "unit": "",
        "description": "Diagonalisation algorithm for the Kohn-Sham Hamiltonian.",
        "valid": ["david", "cg", "ppcg", "paro", "rmm-davidson", "rmm-paro"],
    },
    "diago_thr_init": {
        "default": 0.0,
        "type": "float",
        "unit": "Ry",
        "description": "Threshold for diagonalisation at the first SCF step.",
        "valid": [],
    },
    "diago_cg_maxiter": {
        "default": 20,
        "type": "int",
        "unit": "",
        "description": "Maximum CG iterations in diagonalisation.",
        "valid": [],
    },
    "diago_david_ndim": {
        "default": 4,
        "type": "int",
        "unit": "",
        "description": "Subspace dimension per band in David diagonalisation.",
        "valid": [],
    },
    "diago_full_acc": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Force all empty states to be converged to full accuracy.",
        "valid": [],
    },
    "efield": {
        "default": 0.0,
        "type": "float",
        "unit": "Ry/bohr/e",
        "description": "Homogeneous electric field along efield_cart (a.u.).",
        "valid": [],
    },
    "efield_cart": {
        "default": [0.0, 0.0, 0.0],
        "type": "list[float]",
        "unit": "",
        "description": "Direction of the electric field in Cartesian coordinates.",
        "valid": [],
    },
    "startingpot": {
        "default": "atomic",
        "type": "str",
        "unit": "",
        "description": "Starting potential: superposition of atomic densities or read from file.",
        "valid": ["atomic", "file"],
    },
    "startingwfc": {
        "default": "atomic+random",
        "type": "str",
        "unit": "",
        "description": "Starting wavefunctions.",
        "valid": ["atomic", "atomic+random", "random", "file"],
    },
    "tq_smear": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Enable occupation smearing for exact-exchange calculations.",
        "valid": [],
    },
    "scf_must_converge": {
        "default": True,
        "type": "bool",
        "unit": "",
        "description": "If .false., do not stop on non-convergence of SCF (continue relax).",
        "valid": [],
    },
}

# ---------------------------------------------------------------------------
# &IONS
# ---------------------------------------------------------------------------
IONS = {
    "ion_dynamics": {
        "default": "bfgs",
        "type": "str",
        "unit": "",
        "description": "Algorithm used for ionic minimisation or molecular dynamics.",
        "valid": [
            "bfgs", "damp",
            "verlet", "langevin", "langevin-smc",
            "beeman", "isokinetic",
        ],
    },
    "ion_positions": {
        "default": "default",
        "type": "str",
        "unit": "",
        "description": "How atomic positions are updated in MD.",
        "valid": ["default", "from_input"],
    },
    "pot_extrapolation": {
        "default": "atomic",
        "type": "str",
        "unit": "",
        "description": "Method for extrapolating the potential between MD steps.",
        "valid": ["none", "atomic", "first_order", "second_order"],
    },
    "wfc_extrapolation": {
        "default": "none",
        "type": "str",
        "unit": "",
        "description": "Method for extrapolating wavefunctions between MD steps.",
        "valid": ["none", "first_order", "second_order"],
    },
    "remove_rigid_rot": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Remove rigid body rotation contribution in MD.",
        "valid": [],
    },
    "ion_temperature": {
        "default": "not_controlled",
        "type": "str",
        "unit": "",
        "description": "Thermostat type for MD simulations.",
        "valid": [
            "rescaling", "rescale-v", "rescale-T", "reduce-T",
            "berendsen", "andersen", "svr", "initial", "not_controlled",
        ],
    },
    "tempw": {
        "default": 300.0,
        "type": "float",
        "unit": "K",
        "description": "Target temperature for MD thermostat.",
        "valid": [],
    },
    "tolp": {
        "default": 100.0,
        "type": "float",
        "unit": "K",
        "description": "Tolerance for temperature control (rescaling thermostat).",
        "valid": [],
    },
    "delta_t": {
        "default": 1.0,
        "type": "float",
        "unit": "",
        "description": "Factor for temperature rescaling each nraise steps.",
        "valid": [],
    },
    "nraise": {
        "default": 1,
        "type": "int",
        "unit": "",
        "description": "Rescale temperature every nraise steps.",
        "valid": [],
    },
    "refold_pos": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Refold atomic positions into the first unit cell after each step.",
        "valid": [],
    },
    "upscale": {
        "default": 100.0,
        "type": "float",
        "unit": "",
        "description": "Max factor by which BFGS step can be scaled up.",
        "valid": [],
    },
    "bfgs_ndim": {
        "default": 1,
        "type": "int",
        "unit": "",
        "description": "Number of vectors in the BFGS inverse Hessian update.",
        "valid": [],
    },
    "trust_radius_max": {
        "default": 0.8,
        "type": "float",
        "unit": "bohr",
        "description": "Maximum step allowed in BFGS line search.",
        "valid": [],
    },
    "trust_radius_min": {
        "default": 1e-3,
        "type": "float",
        "unit": "bohr",
        "description": "Minimum step below which BFGS is restarted.",
        "valid": [],
    },
    "trust_radius_ini": {
        "default": 0.5,
        "type": "float",
        "unit": "bohr",
        "description": "Initial BFGS trust radius.",
        "valid": [],
    },
    "w_1": {
        "default": 0.01,
        "type": "float",
        "unit": "",
        "description": "Parameters for the BFGS line search (Wolfe condition 1).",
        "valid": [],
    },
    "w_2": {
        "default": 0.5,
        "type": "float",
        "unit": "",
        "description": "Parameters for the BFGS line search (Wolfe condition 2).",
        "valid": [],
    },
}

# ---------------------------------------------------------------------------
# &CELL
# ---------------------------------------------------------------------------
CELL = {
    "cell_dynamics": {
        "default": "none",
        "type": "str",
        "unit": "",
        "description": "Algorithm for cell relaxation (requires calculation='vc-relax' or 'vc-md').",
        "valid": ["none", "sd", "damp-pr", "damp-w", "bfgs", "pr", "w"],
    },
    "press": {
        "default": 0.0,
        "type": "float",
        "unit": "kbar",
        "description": "Target pressure for variable-cell calculations.",
        "valid": [],
    },
    "wmass": {
        "default": None,
        "type": "float",
        "unit": "a.u.",
        "description": "Fictitious cell mass for Parrinello-Rahman MD.",
        "valid": [],
    },
    "cell_factor": {
        "default": 1.2,
        "type": "float",
        "unit": "",
        "description": "Factor multiplying the initial cell to set the real-space FFT grid.",
        "valid": [],
    },
    "press_conv_thr": {
        "default": 0.5,
        "type": "float",
        "unit": "kbar",
        "description": "Convergence threshold on pressure for vc-relax.",
        "valid": [],
    },
    "cell_dofree": {
        "default": "all",
        "type": "str",
        "unit": "",
        "description": "Which degrees of freedom are free in cell relaxation.",
        "valid": [
            "all", "ibrav", "x", "y", "z", "xy", "xz", "yz",
            "xyz", "shape", "volume", "2Dxy", "2Dshape", "epitaxial_ab",
            "epitaxial_ac", "epitaxial_bc",
        ],
    },
    "cell_velocities": {
        "default": "default",
        "type": "str",
        "unit": "",
        "description": "Whether to reinitialise cell velocities when restarting.",
        "valid": ["default", "zero"],
    },
}

# ---------------------------------------------------------------------------
# CARDS  (each entry is a card template dict, not a flat-parameter dict)
# ---------------------------------------------------------------------------

ATOMIC_SPECIES_CARD = {
    "__description__": (
        "ATOMIC_SPECIES\n"
        "  label  mass  pseudo_file\n"
        "One line per species: symbol, atomic mass (amu), pseudopotential filename."
    ),
    "label": {
        "type": "str",
        "description": "Element symbol or custom label (up to 3 characters).",
    },
    "mass": {
        "type": "float",
        "unit": "amu",
        "description": "Atomic mass in atomic mass units.",
    },
    "pseudo_file": {
        "type": "str",
        "description": "Name of the pseudopotential file (must exist in pseudo_dir).",
    },
}

ATOMIC_POSITIONS_CARD = {
    "__description__": (
        "ATOMIC_POSITIONS { alat | bohr | angstrom | crystal | crystal_sg }\n"
        "  label  x  y  z  [if1  if2  if3]\n"
        "if1/if2/if3 = 1 (free) or 0 (fixed) for constrained relaxation."
    ),
    "units": {
        "type": "str",
        "description": "Units for atomic positions.",
        "valid": ["alat", "bohr", "angstrom", "crystal", "crystal_sg"],
    },
    "label": {"type": "str", "description": "Atomic species label (matches ATOMIC_SPECIES)."},
    "x": {"type": "float", "description": "x coordinate."},
    "y": {"type": "float", "description": "y coordinate."},
    "z": {"type": "float", "description": "z coordinate."},
    "if1": {"type": "int", "description": "1 = free, 0 = fixed along x.", "valid": [0, 1]},
    "if2": {"type": "int", "description": "1 = free, 0 = fixed along y.", "valid": [0, 1]},
    "if3": {"type": "int", "description": "1 = free, 0 = fixed along z.", "valid": [0, 1]},
}

K_POINTS_CARD = {
    "__description__": (
        "K_POINTS { tpiba | automatic | crystal | gamma | tpiba_b | crystal_b | tpiba_c | crystal_c }\n"
        "Defines the k-point sampling mesh."
    ),
    "units": {
        "type": "str",
        "description": "K-point specification mode.",
        "valid": ["tpiba", "automatic", "crystal", "gamma",
                  "tpiba_b", "crystal_b", "tpiba_c", "crystal_c"],
    },
    "nk1": {"type": "int", "description": "# k-points along b1 (automatic mode)."},
    "nk2": {"type": "int", "description": "# k-points along b2 (automatic mode)."},
    "nk3": {"type": "int", "description": "# k-points along b3 (automatic mode)."},
    "sk1": {"type": "int", "description": "Shift along b1: 0 = Γ-centred, 1 = shifted."},
    "sk2": {"type": "int", "description": "Shift along b2."},
    "sk3": {"type": "int", "description": "Shift along b3."},
}

CELL_PARAMETERS_CARD = {
    "__description__": (
        "CELL_PARAMETERS { bohr | angstrom | alat }\n"
        "  v1x  v1y  v1z\n"
        "  v2x  v2y  v2z\n"
        "  v3x  v3y  v3z\n"
        "Required when ibrav=0. Rows are the three lattice vectors."
    ),
    "units": {
        "type": "str",
        "description": "Units for the lattice vectors.",
        "valid": ["bohr", "angstrom", "alat"],
    },
}

CONSTRAINTS_CARD = {
    "__description__": (
        "CONSTRAINTS\n"
        "  nconstr  [ constr_tol ]\n"
        "  type  atom1  atom2  ...  target\n"
        "Defines geometrical constraints during structural relaxation or MD."
    ),
    "nconstr": {"type": "int", "description": "Number of constraints."},
    "constr_tol": {
        "type": "float",
        "unit": "bohr",
        "description": "Tolerance for constraint satisfaction (default 1e-6).",
    },
    "type": {
        "type": "str",
        "description": "Constraint type.",
        "valid": ["'distance'", "'planar_angle'", "'torsional_angle'",
                  "'bennett_constraint'", "'potential_wall'"],
    },
}

OCCUPATIONS_CARD = {
    "__description__": (
        "OCCUPATIONS\n"
        "  f_1  f_2  ...  f_nbnd\n"
        "Explicit occupation numbers when SYSTEM::occupations='from_input'."
    ),
    "f_n": {
        "type": "float",
        "description": "Occupation of band n (between 0 and 1 for collinear spin).",
    },
}

ATOMIC_FORCES_CARD = {
    "__description__": (
        "ATOMIC_FORCES\n"
        "  label  Fx  Fy  Fz\n"
        "External forces applied to atoms (one line per atom, in Ry/bohr)."
    ),
    "label": {"type": "str", "description": "Atomic species label."},
    "Fx": {"type": "float", "unit": "Ry/bohr", "description": "External force component x."},
    "Fy": {"type": "float", "unit": "Ry/bohr", "description": "External force component y."},
    "Fz": {"type": "float", "unit": "Ry/bohr", "description": "External force component z."},
}

# ---------------------------------------------------------------------------
# Convenience: all namelists and cards in one place
# ---------------------------------------------------------------------------
PW_NAMELISTS = {
    "CONTROL": CONTROL,
    "SYSTEM": SYSTEM,
    "ELECTRONS": ELECTRONS,
    "IONS": IONS,
    "CELL": CELL,
}

PW_CARDS = {
    "ATOMIC_SPECIES": ATOMIC_SPECIES_CARD,
    "ATOMIC_POSITIONS": ATOMIC_POSITIONS_CARD,
    "K_POINTS": K_POINTS_CARD,
    "CELL_PARAMETERS": CELL_PARAMETERS_CARD,
    "CONSTRAINTS": CONSTRAINTS_CARD,
    "OCCUPATIONS": OCCUPATIONS_CARD,
    "ATOMIC_FORCES": ATOMIC_FORCES_CARD,
}


PROJWFC = {
    "prefix": {
        "default": "pwscf", "type": "str", "unit": "",
        "description": "Prefix matching the pw.x run whose wavefunctions are projected.",
    },
    "outdir": {
        "default": "./", "type": "str", "unit": "",
        "description": "Directory containing the pw.x save files.",
    },
    "filpdos": {
        "default": "", "type": "str", "unit": "",
        "description": "Prefix for PDOS output files (projwfc.x writes one file per orbital).",
    },
    "filproj": {
        "default": "", "type": "str", "unit": "",
        "description": "Output file prefix for raw |<phi|psi>|^2 projections (fat bands).",
    },
    "Emin": {
        "default": -1000.0, "type": "float", "unit": "eV",
        "description": "Lower bound of the energy window for the PDOS output grid.",
    },
    "Emax": {
        "default": 1000.0, "type": "float", "unit": "eV",
        "description": "Upper bound of the energy window for the PDOS output grid.",
    },
    "DeltaE": {
        "default": 0.01, "type": "float", "unit": "eV",
        "description": "Energy resolution of the PDOS output grid.",
    },
    "kresolveddos": {
        "default": False, "type": "bool", "unit": "",
        "description": "Write k-resolved PDOS (different from fat-band projections).",
    },
}


def describe(namelist: dict, key: str) -> None:
    """Print a formatted description of a single parameter."""
    p = namelist[key]
    print(f"  {key}")
    print(f"    type    : {p.get('type','')}")
    if p.get("unit"):
        print(f"    unit    : {p['unit']}")
    print(f"    default : {p.get('default','—')}")
    if p.get("valid"):
        print(f"    valid   : {p['valid']}")
    print(f"    info    : {p['description']}")
