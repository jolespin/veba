#!/usr/bin/env bash
# __version__ = "2023.01.05"

# Usage: git clone https://github.com/jolespin/veba && update_environment_scripts.sh /path/to/veba_repository
echo "-----------------------------------------------------------------------------------------------------"
echo " Updating module and script files for VEBA"
echo "-----------------------------------------------------------------------------------------------------"

if [ $# -eq 0 ]; then
    echo " * No VEBA repository directory provided.  Downloading repository from: https://github.com/jolespin/veba"
    echo "-----------------------------------------------------------------------------------------------------"
    git clone https://github.com/jolespin/veba
    VEBA_REPOSITORY_DIRECTORY=${PWD}/veba

    # Update permissions
    echo " * Updating permissions: ${VEBA_REPOSITORY_DIRECTORY}"
    chmod 775 ${VEBA_REPOSITORY_DIRECTORY}/src/*
    chmod 775 ${VEBA_REPOSITORY_DIRECTORY}/src/scripts/*

                 else

    VEBA_REPOSITORY_DIRECTORY=$1

fi

CONDA_BASE=$(conda run -n base bash -c "echo \${CONDA_PREFIX}")

echo "-----------------------------------------------------------------------------------------------------"
echo " * Source VEBA: ${VEBA_REPOSITORY_DIRECTORY}"
echo " * Destination VEBA environments CONDA_BASE: ${CONDA_BASE}"
echo "-----------------------------------------------------------------------------------------------------"

# Environemnts
for ENV_PREFIX in ${CONDA_BASE}/envs/VEBA-*; do
    echo $ENV_PREFIX
    cp ${VEBA_REPOSITORY_DIRECTORY}/src/*.py ${ENV_PREFIX}/bin/
    cp -r ${VEBA_REPOSITORY_DIRECTORY}/src/scripts/ ${ENV_PREFIX}/bin/
    ln -sf ${ENV_PREFIX}/bin/scripts/* ${ENV_PREFIX}/bin/
    done
