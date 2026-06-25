# QE Tutorial Notebooks — Context for Claude

Read this before creating any Quantum ESPRESSO tutorial notebook for Google Colab.

## Environment

- QE 7.5 installed in a conda environment named `qe_env` via condacolab
- Colab runs **Python 3.12**; `qe_env` uses Python 3.11
- `qe_env` is saved as a tar archive on Google Drive at:
  `MyDrive/conda_envs/qe_env.tar.gz`
- `ovito` must be installed via **pip** (not conda) because the conda version
  is compiled for Python 3.11 and conflicts with Colab's 3.12 interpreter

## Every tutorial notebook must start with these three cells

### Cell 1 — Configuration (the ONLY cell the user needs to edit)

```python
# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit this cell only
# ══════════════════════════════════════════════════════════════════════════════

# Option 1 — your own Google Drive path (default)
ENV_ARCHIVE = '/content/drive/MyDrive/conda_envs/qe_env.tar.gz'

# Option 2 — shared Google Drive folder
# ENV_ARCHIVE = '/content/drive/Shareddrives/QE_Tutorials/qe_env.tar.gz'

# Option 3 — shared via a Google Drive sharing link
# To use this:
#   1. The owner shares the file with 'Anyone with the link'
#   2. Copy the file ID from the URL:
#      https://drive.google.com/file/d/FILE_ID_HERE/view
#   3. Paste the FILE_ID below and uncomment these two lines:
# GDRIVE_FILE_ID = 'FILE_ID_HERE'
# ENV_ARCHIVE    = f'/content/qe_env.tar.gz'

# ══════════════════════════════════════════════════════════════════════════════
print(f'ENV_ARCHIVE set to: {ENV_ARCHIVE}')
```

### Cell 2 — Bootstrap condacolab

```python
# ── Bootstrap condacolab (triggers kernel restart on first run) ───────────────
# After the restart, re-run this cell once — it will skip the install.
try:
    import condacolab
    condacolab.check()
    print('✅ condacolab active — continuing to restore …')
except Exception:
    import subprocess, sys
    print('Installing condacolab …')
    subprocess.check_call(
        [sys.executable, '-m', 'pip', 'install', '-q', 'condacolab'],
        stdout=subprocess.DEVNULL
    )
    import condacolab
    condacolab.install()    # ← kernel restarts here; re-run this cell after restart
```

### Cell 3 — Restore environment, expose packages, install ovito

```python
# ── Restore qe_env from Drive ─────────────────────────────────────────────────
import subprocess, os

ENV_PATH = '/usr/local/envs/qe_env'

# If using a shared Google Drive link (Option 3), download with gdown first
if 'GDRIVE_FILE_ID' in dir() and not os.path.isfile(ENV_ARCHIVE):
    print('Downloading from Google Drive sharing link …')
    subprocess.check_call(
        ['pip', 'install', '-q', 'gdown'], stdout=subprocess.DEVNULL
    )
    import gdown
    gdown.download(id=GDRIVE_FILE_ID, output=ENV_ARCHIVE, quiet=False)

# If using own/shared Drive path, mount Drive first
elif ENV_ARCHIVE.startswith('/content/drive'):
    from google.colab import drive
    drive.mount('/content/drive')

# Restore the environment
if not os.path.isdir(ENV_PATH):
    print(f'Restoring qe_env from {ENV_ARCHIVE} …')
    os.makedirs(ENV_PATH, exist_ok=True)
    subprocess.run(
        ['tar', '-xzf', ENV_ARCHIVE, '-C', ENV_PATH],
        check=True
    )
    print('✅ Environment restored.')
else:
    print('✅ qe_env already present on disk.')

# ── Expose Python packages to this interpreter ────────────────────────────────
import sys, glob
for sp in glob.glob('/usr/local/envs/qe_env/lib/python*/site-packages'):
    if sp not in sys.path:
        sys.path.insert(0, sp)
        print(f'✅ Added to sys.path: {sp}')

# ── Install ovito via pip (conda version incompatible with Colab Python 3.12) -
print('Installing ovito via pip …')
subprocess.check_call(
    [sys.executable, '-m', 'pip', 'install', '-q', 'ovito'],
    stdout=subprocess.DEVNULL
)

# ── Verify ────────────────────────────────────────────────────────────────────
print('\nPackage availability:')
for pkg in ['numpy', 'matplotlib', 'ase', 'ovito']:
    try:
        __import__(pkg)
        print(f'  ✅  {pkg}')
    except ImportError as e:
        print(f'  ❌  {pkg} — {e}')

print('\nQE executables:')
for exe in ['pw.x', 'ph.x', 'pp.x', 'bands.x', 'dos.x']:
    found = subprocess.run(
        ['conda', 'run', '-n', 'qe_env', 'which', exe],
        capture_output=True, text=True
    )
    status = '✅' if found.returncode == 0 else '❌'
    print(f'  {status}  {exe:12s}  {found.stdout.strip()}')

print('\n🎉 Ready — run QE with: !conda run -n qe_env pw.x < input.in')
```

## Running QE in tutorial cells

Always run QE executables via `conda run`:

```python
# From Python
result = subprocess.run(
    ['conda', 'run', '-n', 'qe_env', 'pw.x', '-input', 'scf.in'],
    capture_output=True, text=True
)

# From shell cell
!conda run -n qe_env pw.x < scf.in | tee scf.out
```

## Pseudopotentials

- Store in `/content/pseudo/`
- Download at the start of each tutorial notebook
- Si pseudopotential URL (confirmed working):
  `https://pseudopotentials.quantum-espresso.org/upf_files/Si.pbe-n-kjpaw_psl.1.0.0.UPF`
- Always set `pseudo_dir = '/content/pseudo'` in the `&CONTROL` namelist

## Important notes

- `pw.x --version` triggers an MPI abort — do not use it for version checks; use `which pw.x` instead
- `condacolab` has no `activate()` method — use `conda run -n qe_env` or `sys.path` manipulation
- The env is located at `/usr/local/envs/qe_env/` (not `/opt/conda/envs/`)
- After `condacolab.install()` the kernel restarts — Cell 2 must be re-run once after restart
