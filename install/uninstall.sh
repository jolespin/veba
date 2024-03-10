#!/bin/bash
# __version__ = "2024.2.19"

#CONDA_BASE=$(conda run -n base bash -c "echo \${CONDA_PREFIX}")
CONDA_ENVS_PATH=${1:-"$(conda info --base)/envs/"}


for FP in ${CONDA_ENVS_PATH}/VEBA ${CONDA_ENVS_PATH}/VEBA-*_env; do 
	ENV_NAME=$(basename $FP)
	mamba env remove -n ${ENV_NAME}
	done

echo -e " _    _ _______ ______  _______\n  \  /  |______ |_____] |_____|\n   \/   |______ |_____] |     |"
echo -e "..............................."
echo -e "     Uninstall Complete     "
echo -e "..............................."
echo -e "Don't forget to remove the VEBA database directory if you don't need it anymore.  If you're doing a reinstall, then think twice about this."
