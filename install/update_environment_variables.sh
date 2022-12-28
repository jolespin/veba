#!/bin/bash
# __VERSION__ = "2022.12.27"

# Create database
DATABASE_DIRECTORY=${1:-"."}
REALPATH_DATABASE_DIRECTORY=$(realpath $DATABASE_DIRECTORY)
# CONDA_BASE=$(which conda | python -c "import sys; print('/'.join(sys.stdin.read().split('/')[:-2]))")
CONDA_BASE=$(conda run -n base bash -c "echo \${CONDA_PREFIX}")

echo ". .. ... ..... ........ ............."
echo "i * Adding the following environment variable to VEBA environments: export VEBA_DATABASE=${REALPATH_DATABASE_DIRECTORY}"


# VEBA
for ENV_PREFIX in ${CONDA_BASE}/envs/VEBA-*; do 
    echo $ENV_PREFIX;
    mkdir -v -p ${ENV_PREFIX}/etc/conda/activate.d/
    mkdir -v -p ${ENV_PREFIX}/etc/conda/deactivate.d/
    echo "export VEBA_DATABASE=${REALPATH_DATABASE_DIRECTORY}" > ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset VEBA_DATABASE" > ${ENV_PREFIX}/etc/conda/deactivate.d/veba.sh
    done

#GTDB-Tk/CheckM
echo ". .. ... ..... ........ ............."
echo "ii * Adding the following environment variable to VEBA environments: export GTDBTK_DATA_PATH=${REALPATH_DATABASE_DIRECTORY}/Classify/GTDBTk/"
for ENV_NAME in VEBA-binning-prokaryotic_env VEBA-classify_env; do 
    ENV_PREFIX=${CONDA_BASE}/envs/${ENV_NAME}
    # GTDB-Tk
    # GTDBTK_DATABASE_VERSION=$(ls ${REALPATH_DATABASE_DIRECTORY}/Classify/GTDBTk)
    echo "export GTDBTK_DATA_PATH=${REALPATH_DATABASE_DIRECTORY}/Classify/GTDBTk/" >> ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset GTDBTK_DATA_PATH" >> ${ENV_PREFIX}/etc/conda/deactivate.d/veba.sh
    # CheckM
    echo "export CHECKM_DATA_PATH=${REALPATH_DATABASE_DIRECTORY}/Classify/CheckM/" >> ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset CHECKM_DATA_PATH" >> ${ENV_PREFIX}/etc/conda/deactivate.d/veba.sh    
    done 

# CheckV
echo ". .. ... ..... ........ ............."
echo "iii * Adding the following environment variable to VEBA environments: export CHECKVDB=${REALPATH_DATABASE_DIRECTORY}/Classify/CheckV/"
for ENV_NAME in VEBA-binning-viral_env; do 
    ENV_PREFIX=${CONDA_BASE}/envs/${ENV_NAME}
    echo "export CHECKVDB=${REALPATH_DATABASE_DIRECTORY}/Classify/CheckV" >> ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset CHECKVDB" >> ${ENV_PREFIX}/etc/conda/deactivate.d/veba.sh
    done

echo -e " _    _ _______ ______  _______\n  \  /  |______ |_____] |_____|\n   \/   |______ |_____] |     |"
echo -e "........................................."
echo -e "  Environment Variable Update Complete   "
echo -e "........................................."
echo -e "The VEBA database environment variable is set in your VEBA conda environments: \n\tVEBA_DATABASE=${REALPATH_DATABASE_DIRECTORY}"