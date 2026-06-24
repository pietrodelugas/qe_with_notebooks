#!/usr/bin/env python3
"""
create_variants.py — generate local/, jupyterlab/, and colab/ notebook variants.

Run from the repo root:
    python create_variants.py

Effect:
  local/         — verbatim copies of the originals (run JupyterLab from repo root)
  jupyterlab/    — conda-env variants: auto-detect modules via sys.path, shutil.which for QE
  colab/         — Colab variants: Drive mount + QE env restore, all paths from DRIVE_DIR
"""

import json
import copy
from pathlib import Path

ROOT = Path(__file__).parent.resolve()

NOTEBOOKS = [
    'qe_pw_input_intro.ipynb',
    'qe_convergence_tests.ipynb',
    'qe_eos_bulkmodulus.ipynb',
    'qe_elastic_constants.ipynb',
    'qe_bandstructure.ipynb',
]

# First markdown cell id in each notebook (bootstrap cell inserted after this)
INTRO_CELL_IDS = {
    'qe_pw_input_intro.ipynb':    'a1b2c3d4',
    'qe_convergence_tests.ipynb': '8266a7d3',
    'qe_eos_bulkmodulus.ipynb':   'intro-cell',
    'qe_elastic_constants.ipynb': 'intro-cell',
    'qe_bandstructure.ipynb':     '06b2da29',
}

# Notebooks that run QE (need PW_CMD, PSEUDO_DIR modifications)
HAS_QE = {nb for nb in NOTEBOOKS if nb != 'qe_pw_input_intro.ipynb'}

# ── helpers ──────────────────────────────────────────────────────────────────

def load_nb(name):
    return json.loads((ROOT / name).read_text())

def save_nb(nb, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(nb, indent=1, ensure_ascii=False) + '\n')

def src_to_list(src):
    """Convert a multiline string to the list-of-lines notebook source format."""
    lines = src.split('\n')
    return [l + '\n' for l in lines[:-1]] + [lines[-1]]

def code_cell(src, cell_id):
    return {
        'cell_type': 'code',
        'execution_count': None,
        'id': cell_id,
        'metadata': {},
        'outputs': [],
        'source': src_to_list(src),
    }

def find_idx(nb, cell_id):
    for i, c in enumerate(nb['cells']):
        if c.get('id') == cell_id:
            return i
    raise KeyError(f'cell id {cell_id!r} not found')

def insert_after(nb, cell_id, new_cell):
    nb['cells'].insert(find_idx(nb, cell_id) + 1, new_cell)

def replace_source(nb, cell_id, src):
    i = find_idx(nb, cell_id)
    nb['cells'][i]['source'] = src_to_list(src)
    nb['cells'][i]['outputs'] = []
    nb['cells'][i]['execution_count'] = None

# ── shared bootstrap cell content ────────────────────────────────────────────

JL_SYSPATH = """\
import sys
from pathlib import Path

# Modules live at the repo root.
# When running from jupyterlab/ in the repo: parent dir contains pw_input.py.
# When running from the self-contained drive package: current dir has them.
_parent = Path('..').resolve()
_here   = Path('.').resolve()
_module_dir = _parent if (_parent / 'pw_input.py').exists() else _here
if str(_module_dir) not in sys.path:
    sys.path.insert(0, str(_module_dir))"""

COLAB_BOOTSTRAP = """\
from google.colab import drive
import sys, subprocess
from pathlib import Path

drive.mount('/content/drive')

DRIVE_DIR  = Path('/content/drive/MyDrive/qe_tutorial_colab')
QE_ENV_DIR = Path('/content/qe_env')

if str(DRIVE_DIR) not in sys.path:
    sys.path.insert(0, str(DRIVE_DIR))

# Restore QE environment (built once by qe_environment_setup.ipynb, saved to Drive)
if not QE_ENV_DIR.exists():
    print('Restoring QE environment (~1 min)...')
    subprocess.run(
        ['tar', '-xzf', str(DRIVE_DIR / 'qe_env.tar.gz'), '-C', '/content'],
        check=True,
    )
    print('Done.')
else:
    print('QE env already extracted.')

QE_BIN = QE_ENV_DIR / 'bin'"""

COLAB_BOOTSTRAP_NO_QE = """\
from google.colab import drive
import sys
from pathlib import Path

drive.mount('/content/drive')

DRIVE_DIR = Path('/content/drive/MyDrive/qe_tutorial_colab')

if str(DRIVE_DIR) not in sys.path:
    sys.path.insert(0, str(DRIVE_DIR))

print('Setup complete.')"""

# ── per-notebook setup cell replacements ────────────────────────────────────

SETUP = {
    'qe_convergence_tests.ipynb': {
        'cell_id': 'e0b9b4c1',

        'jupyterlab': """\
from pathlib import Path
import shutil, os

os.environ['OMP_NUM_THREADS'] = '1'

RUN_ROOT = Path('.').resolve()

_pw = shutil.which('pw.x')
if _pw is None:
    raise RuntimeError('pw.x not found — activate the qe_env conda environment.')
PW_CMD = [_pw]

PSEUDO_DIR = _module_dir / 'pseudo'
OUT_DIR    = RUN_ROOT / 'out'
ECUT_DIR   = RUN_ROOT / 'convergence' / 'ecut'
KPTS_DIR   = RUN_ROOT / 'convergence' / 'kpoints'

for d in [OUT_DIR, ECUT_DIR, KPTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

print('pw.x:', PW_CMD[0])""",

        'colab': """\
import os
from pathlib import Path

os.environ['OMP_NUM_THREADS'] = '1'

RUN_ROOT = Path('/content')
PW_CMD   = [str(QE_BIN / 'pw.x')]

PSEUDO_DIR = DRIVE_DIR / 'pseudo'
OUT_DIR    = RUN_ROOT / 'out'
ECUT_DIR   = RUN_ROOT / 'convergence' / 'ecut'
KPTS_DIR   = RUN_ROOT / 'convergence' / 'kpoints'

for d in [OUT_DIR, ECUT_DIR, KPTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

print('pw.x:', PW_CMD[0])""",
    },

    'qe_eos_bulkmodulus.ipynb': {
        'cell_id': '645e57e1',

        'jupyterlab': """\
from pathlib import Path
import shutil, os
import numpy as np
from ase.build import bulk

from pw_input import (
    ControlNamelist, SystemNamelist, ElectronsNamelist,
    AtomicSpeciesCard, AtomicPositionsCard, KPointsAutoCard, PWInput,
)
from convergence_runner import QERunner, RY_TO_EV, BOHR_TO_ANG

os.environ['OMP_NUM_THREADS'] = '1'

RUN_ROOT   = Path('.').resolve()
PSEUDO_DIR = _module_dir / 'pseudo'
EOS_DIR    = RUN_ROOT / 'out' / 'eos'
EOS_DIR.mkdir(parents=True, exist_ok=True)

_pw = shutil.which('pw.x')
if _pw is None:
    raise RuntimeError('pw.x not found — activate the qe_env conda environment.')
PW_CMD = [_pw]

PSEUDOS = {
    'Mg': 'Mg.upf',
    'O':  'O.upf',
}""",

        'colab': """\
from pathlib import Path
import os
import numpy as np
from ase.build import bulk

from pw_input import (
    ControlNamelist, SystemNamelist, ElectronsNamelist,
    AtomicSpeciesCard, AtomicPositionsCard, KPointsAutoCard, PWInput,
)
from convergence_runner import QERunner, RY_TO_EV, BOHR_TO_ANG

os.environ['OMP_NUM_THREADS'] = '1'

RUN_ROOT   = Path('/content')
PSEUDO_DIR = DRIVE_DIR / 'pseudo'
EOS_DIR    = RUN_ROOT / 'out' / 'eos'
EOS_DIR.mkdir(parents=True, exist_ok=True)

PW_CMD = [str(QE_BIN / 'pw.x')]

PSEUDOS = {
    'Mg': 'Mg.upf',
    'O':  'O.upf',
}""",
    },

    'qe_elastic_constants.ipynb': {
        'cell_id': 'setup-cell',

        'jupyterlab': """\
from pathlib import Path
import shutil
import numpy as np
from ase.build import bulk

from convergence_runner import QERunner, RY_TO_EV
from elastic_tools import (
    build_ortho_input, build_mono_input,
    fit_elastic_quad, extract_c11_c12, extract_c44,
)
import os
os.environ['OMP_NUM_THREADS'] = '1'
RUN_ROOT   = Path('.').resolve()
PSEUDO_DIR = _module_dir / 'pseudo'
ELAST_DIR  = RUN_ROOT / 'out' / 'elastic'
ELAST_DIR.mkdir(parents=True, exist_ok=True)

_pw = shutil.which('pw.x')
if _pw is None:
    raise RuntimeError('pw.x not found — activate the qe_env conda environment.')
PW_CMD = [_pw]

PSEUDOS = {
    'Mg': 'Mg.upf',
    'O':  'O.upf',
}

# ---- Results from qe_eos_bulkmodulus.ipynb (Birch-Murnaghan fit) ----------
# Fill in your values here before running the elastic constant sweeps.
A0_ANG       = 4.21    # Å  — equilibrium lattice parameter
B0_GPA       = 160.0   # GPa — bulk modulus
V0_PRIM_ANG3 = (A0_ANG**3) / 4   # Å³ — primitive cell volume (2 atoms)
V0_CONV_ANG3 = 4 * V0_PRIM_ANG3  # Å³ — conventional cell volume (8 atoms)

ECUTWFC_CONV = 60   # Ry  — from convergence notebook
NK_ELASTIC   = 4    # k-grid along each direction (equivalent density to nk=8 for primitive cell)
FORCE_RERUN  = False

runner = QERunner(PW_CMD)

print(f'a₀ = {A0_ANG:.4f} Å,  B₀ = {B0_GPA:.1f} GPa')
print(f'V₀ (primitive) = {V0_PRIM_ANG3:.4f} Å³,  V₀ (conventional, 8 atoms) = {V0_CONV_ANG3:.4f} Å³')""",

        'colab': """\
from pathlib import Path
import numpy as np
from ase.build import bulk

from convergence_runner import QERunner, RY_TO_EV
from elastic_tools import (
    build_ortho_input, build_mono_input,
    fit_elastic_quad, extract_c11_c12, extract_c44,
)
import os
os.environ['OMP_NUM_THREADS'] = '1'
RUN_ROOT   = Path('/content')
PSEUDO_DIR = DRIVE_DIR / 'pseudo'
ELAST_DIR  = RUN_ROOT / 'out' / 'elastic'
ELAST_DIR.mkdir(parents=True, exist_ok=True)

PW_CMD = [str(QE_BIN / 'pw.x')]

PSEUDOS = {
    'Mg': 'Mg.upf',
    'O':  'O.upf',
}

# ---- Results from qe_eos_bulkmodulus.ipynb (Birch-Murnaghan fit) ----------
# Fill in your values here before running the elastic constant sweeps.
A0_ANG       = 4.21    # Å  — equilibrium lattice parameter
B0_GPA       = 160.0   # GPa — bulk modulus
V0_PRIM_ANG3 = (A0_ANG**3) / 4   # Å³ — primitive cell volume (2 atoms)
V0_CONV_ANG3 = 4 * V0_PRIM_ANG3  # Å³ — conventional cell volume (8 atoms)

ECUTWFC_CONV = 60   # Ry  — from convergence notebook
NK_ELASTIC   = 4    # k-grid along each direction (equivalent density to nk=8 for primitive cell)
FORCE_RERUN  = False

runner = QERunner(PW_CMD)

print(f'a₀ = {A0_ANG:.4f} Å,  B₀ = {B0_GPA:.1f} GPa')
print(f'V₀ (primitive) = {V0_PRIM_ANG3:.4f} Å³,  V₀ (conventional, 8 atoms) = {V0_CONV_ANG3:.4f} Å³')""",
    },

    'qe_bandstructure.ipynb': {
        'cell_id': '1aa83aa7',

        'jupyterlab': """\
from pathlib import Path
import shutil, subprocess, os
import numpy as np
import matplotlib.pyplot as plt
from ase.build import bulk

from pw_input import (
    ControlNamelist, SystemNamelist, ElectronsNamelist,
    AtomicSpeciesCard, AtomicPositionsCard, KPointsAutoCard, PWInput,
    ProjwfcNamelist,
)
from convergence_runner import QERunner, RY_TO_EV
from bandstructure_tools import (
    parse_fermi_energy, parse_hs_positions, parse_gamma_symmetries,
    read_bands_gnu, build_kpath_str, print_bands_summary,
    plot_band_structure, read_projwfc_weights,
)

os.environ['OMP_NUM_THREADS'] = '1'

RUN_ROOT   = Path('.').resolve()
PSEUDO_DIR = _module_dir / 'pseudo'
BS_DIR     = RUN_ROOT / 'out' / 'bandstructure'
BS_DIR.mkdir(parents=True, exist_ok=True)

def _find_cmd(name):
    p = shutil.which(name)
    if p is None:
        raise RuntimeError(f'{name} not found — activate the qe_env conda environment.')
    return [p]

PW_CMD      = _find_cmd('pw.x')
BANDS_CMD   = _find_cmd('bands.x')
DOS_CMD     = _find_cmd('dos.x')
PROJWFC_CMD = _find_cmd('projwfc.x')

PSEUDOS = {
    'Mg': 'Mg.upf',
    'O':  'O.upf',
}

# Converged parameters from qe_convergence_tests.ipynb and qe_eos_bulkmodulus.ipynb
ECUTWFC = 60      # Ry
NK_SCF  = 8       # nk×nk×nk for SCF (FCC primitive cell)
A0      = 4.212   # Å — equilibrium lattice parameter from BM EOS fit; update with your value
NBND    = 16      # bands to compute (≥2× the 4 occupied bands)""",

        'colab': """\
from pathlib import Path
import subprocess, os
import numpy as np
import matplotlib.pyplot as plt
from ase.build import bulk

from pw_input import (
    ControlNamelist, SystemNamelist, ElectronsNamelist,
    AtomicSpeciesCard, AtomicPositionsCard, KPointsAutoCard, PWInput,
    ProjwfcNamelist,
)
from convergence_runner import QERunner, RY_TO_EV
from bandstructure_tools import (
    parse_fermi_energy, parse_hs_positions, parse_gamma_symmetries,
    read_bands_gnu, build_kpath_str, print_bands_summary,
    plot_band_structure, read_projwfc_weights,
)

os.environ['OMP_NUM_THREADS'] = '1'

RUN_ROOT   = Path('/content')
PSEUDO_DIR = DRIVE_DIR / 'pseudo'
BS_DIR     = RUN_ROOT / 'out' / 'bandstructure'
BS_DIR.mkdir(parents=True, exist_ok=True)

PW_CMD      = [str(QE_BIN / 'pw.x')]
BANDS_CMD   = [str(QE_BIN / 'bands.x')]
DOS_CMD     = [str(QE_BIN / 'dos.x')]
PROJWFC_CMD = [str(QE_BIN / 'projwfc.x')]

PSEUDOS = {
    'Mg': 'Mg.upf',
    'O':  'O.upf',
}

# Converged parameters from qe_convergence_tests.ipynb and qe_eos_bulkmodulus.ipynb
ECUTWFC = 60      # Ry
NK_SCF  = 8       # nk×nk×nk for SCF (FCC primitive cell)
A0      = 4.212   # Å — equilibrium lattice parameter from BM EOS fit; update with your value
NBND    = 16      # bands to compute (≥2× the 4 occupied bands)""",
    },
}

# ── main ─────────────────────────────────────────────────────────────────────

for d in ['local', 'jupyterlab', 'colab']:
    (ROOT / d).mkdir(exist_ok=True)

for nb_name in NOTEBOOKS:
    print(f'{nb_name}')
    nb = load_nb(nb_name)
    intro_id = INTRO_CELL_IDS[nb_name]
    has_qe   = nb_name in HAS_QE

    # local/ — verbatim copy
    save_nb(copy.deepcopy(nb), ROOT / 'local' / nb_name)

    # jupyterlab/ — insert sys.path cell; update setup cell
    nb_jl = copy.deepcopy(nb)
    insert_after(nb_jl, intro_id, code_cell(JL_SYSPATH, 'jl-syspath-cell'))
    if has_qe:
        sc = SETUP[nb_name]
        replace_source(nb_jl, sc['cell_id'], sc['jupyterlab'])
    save_nb(nb_jl, ROOT / 'jupyterlab' / nb_name)

    # colab/ — insert bootstrap cell; update setup cell
    nb_co = copy.deepcopy(nb)
    bootstrap = COLAB_BOOTSTRAP if has_qe else COLAB_BOOTSTRAP_NO_QE
    insert_after(nb_co, intro_id, code_cell(bootstrap, 'colab-bootstrap-cell'))
    if has_qe:
        sc = SETUP[nb_name]
        replace_source(nb_co, sc['cell_id'], sc['colab'])
    save_nb(nb_co, ROOT / 'colab' / nb_name)

    print(f'  local/ jupyterlab/ colab/  ✓')

print('\nDone.')
