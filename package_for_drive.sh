#!/usr/bin/env bash
# package_for_drive.sh — assemble the Google Drive distribution packages.
#
# Produces two self-contained directories under drive_package/:
#   drive_package/jupyterlab/  → share as 'qe_tutorial_jupyterlab' on Google Drive
#   drive_package/colab/       → share as 'qe_tutorial_colab'       on Google Drive
#
# Students set DRIVE_DIR to the mounted path of whichever folder they download.
# Colab students also need qe_env.tar.gz in the same Drive folder (built by
# qe_environment_setup.ipynb).
#
# Run from the repo root after create_variants.py has been executed:
#   bash package_for_drive.sh

set -euo pipefail

REPO="$(cd "$(dirname "$0")" && pwd)"
PKG="$REPO/drive_package"

MODULES=(
    pw_input.py
    pw_namelists.py
    convergence_runner.py
    convergence_plotting.py
    convergence_analysis.py
    eos_tools.py
    eos_plotting.py
    elastic_tools.py
    bandstructure_tools.py
    postproc_namelists.py
)

TUTORIAL_NOTEBOOKS=(
    qe_pw_input_intro
    qe_convergence_tests
    qe_scf_conv_thr
    qe_eos_bulkmodulus
    qe_elastic_constants
    qe_bandstructure
)

for variant in jupyterlab colab; do
    dst="$PKG/$variant"
    rm -rf "$dst"
    mkdir -p "$dst"

    # Python modules
    for mod in "${MODULES[@]}"; do
        cp "$REPO/$mod" "$dst/"
    done

    # Pseudopotentials
    cp -r "$REPO/pseudo" "$dst/pseudo"

    # Notebooks from the appropriate variant directory
    for nb in "${TUTORIAL_NOTEBOOKS[@]}"; do
        cp "$REPO/$variant/${nb}.ipynb" "$dst/"
    done

    echo "Built $dst"
done

# environment.yml goes in the jupyterlab package only (conda-specific)
cp "$REPO/environment.yml" "$PKG/jupyterlab/"

echo ""
echo "Packages ready in drive_package/:"
echo "  jupyterlab/  →  share as 'qe_tutorial_jupyterlab' on Google Drive"
echo "  colab/       →  share as 'qe_tutorial_colab'       on Google Drive"
echo ""
echo "Colab students also need qe_env.tar.gz in 'qe_tutorial_colab/' on Drive."
echo "  (Built once by qe_environment_setup.ipynb — see that notebook for instructions.)"
echo ""
echo "Drive path used in colab notebooks:"
echo "  DRIVE_DIR = Path('/content/drive/MyDrive/qe_tutorial_colab')"
