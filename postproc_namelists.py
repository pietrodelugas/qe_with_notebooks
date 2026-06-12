"""
Quantum ESPRESSO post-processing codes — namelists and cards reference dictionaries.

Codes covered:
    pp.x      — charge density, potential, and wavefunction post-processing
    bands.x   — band structure post-processing and symmetry analysis

For the phonon workflow codes (ph.x, q2r.x, matdyn.x) see ph_namelists.py.
For PDOS and total DOS codes (projwfc.x, dos.x) see pdos_dos_namelists.py.

Same schema as pw_namelists.py:
    default     : default value (None if required / no default)
    type        : Python type hint string
    unit        : physical unit or \'\' if dimensionless/not applicable
    description : concise human-readable explanation
    valid       : list of valid choices (empty list = free-form)
"""

# ============================================================================
# pp.x
# ============================================================================

PP_INPUTPP = {
    "prefix": {
        "default": "pwscf",
        "type": "str",
        "unit": "",
        "description": "Must match prefix used in the pw.x calculation.",
        "valid": [],
    },
    "outdir": {
        "default": "./",
        "type": "str",
        "unit": "",
        "description": "Directory containing pw.x save files.",
        "valid": [],
    },
    "filplot": {
        "default": "tmp_pp",
        "type": "str",
        "unit": "",
        "description": "Output file name for the intermediate plot data.",
        "valid": [],
    },
    "plot_num": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": (
            "Selects the quantity to plot:\n"
            "  0  = charge density\n"
            "  1  = total potential (V_bare + V_H + V_xc)\n"
            "  2  = local ionic potential\n"
            "  3  = local density of states at E_F\n"
            "  4  = local density of electronic entropy\n"
            "  5  = STM image (Tersoff-Hamann)\n"
            "  6  = spin polarisation (rho_up - rho_dn)\n"
            "  7  = |psi|^2 for a given band/k-point\n"
            "  8  = electron localisation function (ELF)\n"
            "  9  = charge density minus superposition of atomic densities\n"
            " 10  = integrated LDOS (Tersoff-Hamann, up to E_F)\n"
            " 11  = V_bare (local pseudopotential)\n"
            " 12  = SAHAR electrostatic potential (ESM)\n"
            " 13  = noncollinear magnetisation (vector)\n"
            " 17  = all-electron charge density (PAW)\n"
            " 18  = exchange-correlation potential\n"
            " 19  = reduced gradient |∇ρ|/(2*(3π²ρ)^{1/3}ρ) for NCI\n"
            " 20  = product s²·sign(λ2)·ρ for NCI isosurface\n"
            " 21  = all-electron + valence charge (PAW)\n"
        ),
        "valid": list(range(0, 22)),
    },
    "spin_component": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "Spin component: 0=total, 1=up, 2=down (nspin=2); or 1/2/3 for x/y/z (nspin=4).",
        "valid": [0, 1, 2, 3],
    },
    "kband": {
        "default": 1,
        "type": "int",
        "unit": "",
        "description": "Band index for plot_num=7 (|ψ|²).",
        "valid": [],
    },
    "kpoint": {
        "default": 1,
        "type": "int",
        "unit": "",
        "description": "K-point index for plot_num=7 (|ψ|²).",
        "valid": [],
    },
    "lsign": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "If .true., save the sign of the wavefunction (plot_num=7).",
        "valid": [],
    },
    "emin": {
        "default": None,
        "type": "float",
        "unit": "eV",
        "description": "Lower energy bound for LDOS integration (plot_num=3,10).",
        "valid": [],
    },
    "emax": {
        "default": None,
        "type": "float",
        "unit": "eV",
        "description": "Upper energy bound for LDOS integration (plot_num=3,10).",
        "valid": [],
    },
    "sample_bias": {
        "default": -0.01,
        "type": "float",
        "unit": "Ry",
        "description": "STM sample bias voltage (plot_num=5,10).",
        "valid": [],
    },
    "stm_wfc_matching": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Perform STM wavefunction matching to the vacuum.",
        "valid": [],
    },
    "z": {
        "default": 1.0,
        "type": "float",
        "unit": "alat",
        "description": "Z-plane height for STM matching (plot_num=5).",
        "valid": [],
    },
    "dz": {
        "default": 0.05,
        "type": "float",
        "unit": "alat",
        "description": "Step size for vacuum extrapolation in STM (plot_num=5).",
        "valid": [],
    },
    "weight_factor": {
        "default": 1.0,
        "type": "float",
        "unit": "",
        "description": "Weighting factor applied to the plotted quantity.",
        "valid": [],
    },
    "iflag": {
        "default": None,
        "type": "int",
        "unit": "",
        "description": "Internal flag for spin/orbital decomposition.",
        "valid": [],
    },
}

PP_PLOT = {
    "nfile": {
        "default": 1,
        "type": "int",
        "unit": "",
        "description": "Number of intermediate plot files to combine.",
        "valid": [],
    },
    "filepp": {
        "default": None,
        "type": "list[str]",
        "unit": "",
        "description": "List of filplot names to combine (length nfile).",
        "valid": [],
    },
    "weight": {
        "default": None,
        "type": "list[float]",
        "unit": "",
        "description": "Weights for each filplot file when combining (length nfile).",
        "valid": [],
    },
    "iflag": {
        "default": 3,
        "type": "int",
        "unit": "",
        "description": (
            "Dimensionality of the plot:\n"
            "  0 = 1D along a line\n"
            "  1 = 1D spherical average\n"
            "  2 = 2D contour plot\n"
            "  3 = 3D volumetric (cube/xsf/xcrysden)\n"
            "  4 = 2D polar plot on a sphere"
        ),
        "valid": [0, 1, 2, 3, 4],
    },
    "output_format": {
        "default": 5,
        "type": "int",
        "unit": "",
        "description": (
            "Output file format:\n"
            "  0 = gnuplot (1D/2D)\n"
            "  1 = contour.x (2D)\n"
            "  2 = plotrho.x (2D)\n"
            "  3 = XCrysDen (.xsf)\n"
            "  4 = gOpenMol (.plt)\n"
            "  5 = Gaussian cube (.cube)\n"
            "  6 = gnuplot (2D polar)\n"
            "  7 = gnuplot (2D, arbitrary plane)\n"
        ),
        "valid": [0, 1, 2, 3, 4, 5, 6, 7],
    },
    "fileout": {
        "default": "tmp_plot",
        "type": "str",
        "unit": "",
        "description": "Output file name for the plot.",
        "valid": [],
    },
    "interpolation": {
        "default": "fourier",
        "type": "str",
        "unit": "",
        "description": "Interpolation method for the plot.",
        "valid": ["fourier", "bspline"],
    },
    "x0": {
        "default": [0.0, 0.0, 0.0],
        "type": "list[float]",
        "unit": "alat",
        "description": "Origin of the plot (3-vector in alat units, 2D/3D).",
        "valid": [],
    },
    "e1": {
        "default": None,
        "type": "list[float]",
        "unit": "alat",
        "description": "First plot direction vector.",
        "valid": [],
    },
    "e2": {
        "default": None,
        "type": "list[float]",
        "unit": "alat",
        "description": "Second plot direction vector (2D/3D).",
        "valid": [],
    },
    "e3": {
        "default": None,
        "type": "list[float]",
        "unit": "alat",
        "description": "Third plot direction vector (3D).",
        "valid": [],
    },
    "nx": {
        "default": None,
        "type": "int",
        "unit": "",
        "description": "Number of grid points along e1.",
        "valid": [],
    },
    "ny": {
        "default": None,
        "type": "int",
        "unit": "",
        "description": "Number of grid points along e2.",
        "valid": [],
    },
    "nz": {
        "default": None,
        "type": "int",
        "unit": "",
        "description": "Number of grid points along e3.",
        "valid": [],
    },
    "radius": {
        "default": 1.0,
        "type": "float",
        "unit": "alat",
        "description": "Radius for spherical/polar plot (iflag=1 or 4).",
        "valid": [],
    },
}

PP_NAMELISTS = {"INPUTPP": PP_INPUTPP, "PLOT": PP_PLOT}


# ============================================================================
# bands.x
# ============================================================================

BANDS_INPUTBANDS = {
    "prefix": {
        "default": "pwscf",
        "type": "str",
        "unit": "",
        "description": "Must match prefix of the pw.x bands calculation.",
        "valid": [],
    },
    "outdir": {
        "default": "./",
        "type": "str",
        "unit": "",
        "description": "Directory containing pw.x save files.",
        "valid": [],
    },
    "filband": {
        "default": "bands.dat",
        "type": "str",
        "unit": "",
        "description": "Output file for band energies (processed by plotband.x).",
        "valid": [],
    },
    "spin_component": {
        "default": 1,
        "type": "int",
        "unit": "",
        "description": "Spin component to plot: 1=up (or unpolarised), 2=down.",
        "valid": [1, 2],
    },
    "firstk": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "First k-point to process (0 = start from the beginning).",
        "valid": [],
    },
    "lastk": {
        "default": 0,
        "type": "int",
        "unit": "",
        "description": "Last k-point to process (0 = all k-points).",
        "valid": [],
    },
    "no_overlap": {
        "default": True,
        "type": "bool",
        "unit": "",
        "description": "Do not compute wavefunction overlaps for band sorting.",
        "valid": [],
    },
    "plot_2d": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Output for 2D colour-coded band plot.",
        "valid": [],
    },
    "lsym": {
        "default": True,
        "type": "bool",
        "unit": "",
        "description": "Compute band symmetry character (requires symmetry ops).",
        "valid": [],
    },
    "filp": {
        "default": "p_avg.dat",
        "type": "str",
        "unit": "",
        "description": "File for the projected (fat) bands output.",
        "valid": [],
    },
    "lp": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Compute the average over k-point pairs (for group velocity).",
        "valid": [],
    },
    "spin_mask": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Apply a spin mask when computing expectation values.",
        "valid": [],
    },
    "rhombohedral_to_hexagonal": {
        "default": False,
        "type": "bool",
        "unit": "",
        "description": "Convert band path from rhombohedral to hexagonal axes.",
        "valid": [],
    },
}

BANDS_NAMELISTS = {"BANDS": BANDS_INPUTBANDS}



# ============================================================================
# Convenience bundle
# ============================================================================
POSTPROC_NAMELISTS = {
    "pp.x":    {"namelists": PP_NAMELISTS},
    "bands.x": {"namelists": BANDS_NAMELISTS},
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
