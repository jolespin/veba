#!/bin/bash

SCRIPT_PATH=$(realpath $0)
PREFIX=$(echo $SCRIPT_PATH | python -c "import sys; print('/'.join(sys.stdin.read().split('/')[:-1]))")
CONDA_BASE=$(conda run -n base bash -c "echo \${CONDA_PREFIX}")

# Environemnts
for ENV_YAML in ${PREFIX}/environments/VEBA*.yml; do
    ENV_NAME=$(basename $ENV_YAML .yml)
    echo "Creating ${ENV_NAME} environment"
    conda env create -n $ENV_NAME -f $ENV_YAML
    cp ${PREFIX}/../src/*.py ${CONDA_BASE}/${ENV_NAME}/bin/
    cp ${PREFIX}/../src/scripts/ ${CONDA_BASE}/${ENV_NAME}/bin/
    ln -sf ${PREFIX}/../src/scripts/* ${CONDA_BASE}/${ENV_NAME}/bin/
    done

echo "Please run `download_databases.sh` script"