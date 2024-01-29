#!/bin/bash
# __version__ = "2024.1.22"

CONDA_BASE=$(conda run -n base bash -c "echo \${CONDA_PREFIX}")

for FP in ${CONDA_BASE}/envs/VEBA ${CONDA_BASE}/envs/VEBA-*_env; do 
	ENV_NAME=$(basename $FP)
	mamba env remove -n ${ENV_NAME}
	done

echo -e " _    _ _______ ______  _______\n  \  /  |______ |_____] |_____|\n   \/   |______ |_____] |     |"
echo -e "..............................."
echo -e "     Uninstall Complete     "
echo -e "..............................."
echo -e "Don't forget to remove the VEBA database directory if you don't need it anymore.  If you're doing a reinstall, then think twice about this."
