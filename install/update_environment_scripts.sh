#!/usr/bin/env bash
# __version__ = "2024.3.5"

# Usage: git clone https://github.com/jolespin/veba && update_environment_scripts.sh /path/to/veba_repository
echo "-----------------------------------------------------------------------------------------------------"
echo " Updating module and script files for VEBA"
echo "-----------------------------------------------------------------------------------------------------"


SCRIPT_DIRECTORY=$(dirname $0)
VEBA_REPOSITORY_DIRECTORY=$(realpath ${SCRIPT_DIRECTORY}/../)
VEBA_REPOSITORY_DIRECTORY=${1:-"$VEBA_REPOSITORY_DIRECTORY"}

# if [ $# -eq 0 ]; then
#     SCRIPT_DIRECTORY=$(dirname $0)
#     echo " * No VEBA repository directory provided.  Using ${VEBA_REPOSITORY_DIRECTORY}"
#     # Update permissions
#     echo " * Updating permissions: ${VEBA_REPOSITORY_DIRECTORY}"
#     chmod 755 ${VEBA_REPOSITORY_DIRECTORY}/bin/*
#     chmod 755 ${VEBA_REPOSITORY_DIRECTORY}/bin/scripts/*
#     else
#     VEBA_REPOSITORY_DIRECTORY=$1
# fi

# CONDA_BASE=$(conda info --base)
CONDA_ENVS_PATH=${2:-"$(conda info --base)/envs/"}

echo "-----------------------------------------------------------------------------------------------------"
echo " * Source VEBA: ${VEBA_REPOSITORY_DIRECTORY}"
echo " * Destination VEBA environments CONDA_ENVS_PATH: ${CONDA_ENVS_PATH}"
echo "-----------------------------------------------------------------------------------------------------"

# Environemnts
for ENV_PREFIX in ${CONDA_ENVS_PATH}/VEBA ${CONDA_ENVS_PATH}/VEBA-*; 
do
    echo $ENV_PREFIX
    cp ${VEBA_REPOSITORY_DIRECTORY}/bin/*.py ${ENV_PREFIX}/bin/
    cp -r ${VEBA_REPOSITORY_DIRECTORY}/bin/scripts/ ${ENV_PREFIX}/bin/

    DST=${ENV_PREFIX}/bin/
    for SRC in ${ENV_PREFIX}/bin/scripts/*;
        do
        SRC=$(realpath --relative-to ${DST} ${SRC})
        ln -sf ${SRC} ${DST}
        done

    # Version
    cp -rf ${VEBA_REPOSITORY_DIRECTORY}/VERSION ${ENV_PREFIX}/bin/VEBA_VERSION
    cp -r ${VEBA_REPOSITORY_DIRECTORY}/bin/VEBA_SCRIPT_VERSIONS ${ENV_PREFIX}/bin/
done

