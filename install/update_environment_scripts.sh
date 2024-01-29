#!/usr/bin/env bash
# __version__ = "2024.1.23"

# Usage: git clone https://github.com/jolespin/veba && update_environment_scripts.sh /path/to/veba_repository
echo "-----------------------------------------------------------------------------------------------------"
echo " Updating module and script files for VEBA"
echo "-----------------------------------------------------------------------------------------------------"

if [ $# -eq 0 ]; then
    SCRIPT_DIRECTORY=$(dirname $0)
    VEBA_REPOSITORY_DIRECTORY=$(realpath ${SCRIPT_DIRECTORY}/../)
    echo " * No VEBA repository directory provided.  Using ${VEBA_REPOSITORY_DIRECTORY}"

    # Update permissions
    echo " * Updating permissions: ${VEBA_REPOSITORY_DIRECTORY}"
    chmod 755 ${VEBA_REPOSITORY_DIRECTORY}/bin/*
    chmod 755 ${VEBA_REPOSITORY_DIRECTORY}/bin/scripts/*

    else

    VEBA_REPOSITORY_DIRECTORY=$1

fi

# CONDA_BASE=$(conda run -n base bash -c "echo \${CONDA_PREFIX}")
CONDA_BASE=$(conda info --base)

echo "-----------------------------------------------------------------------------------------------------"
echo " * Source VEBA: ${VEBA_REPOSITORY_DIRECTORY}"
echo " * Destination VEBA environments CONDA_BASE: ${CONDA_BASE}"
echo "-----------------------------------------------------------------------------------------------------"

# Environemnts
for ENV_PREFIX in ${CONDA_BASE}/envs/VEBA ${CONDA_BASE}/envs/VEBA-*; 
do
    echo $ENV_PREFIX
    cp ${VEBA_REPOSITORY_DIRECTORY}/bin/*.py ${ENV_PREFIX}/bin/
    cp -r ${VEBA_REPOSITORY_DIRECTORY}/bin/scripts/ ${ENV_PREFIX}/bin/
    ln -sf ${ENV_PREFIX}/bin/scripts/* ${ENV_PREFIX}/bin/

    # Version
    cp -rf ${VEBA_REPOSITORY_DIRECTORY}/VERSION ${ENV_PREFIX}/bin/VEBA_VERSION
    cp -r ${VEBA_REPOSITORY_DIRECTORY}/bin/VEBA_SCRIPT_VERSIONS ${ENV_PREFIX}/bin/
done

