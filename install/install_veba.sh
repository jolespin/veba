#!/bin/bash

PREFIX=$(cat $0 | python -c "import sys; print('/'.join(sys.stdin.read().split('/')[:-1]))")
CONDA_BASE=$(which conda | python -c "import sys; print('/'.join(sys.stdin.read().split('/')[:-2]))")

# Environemnts
for ENV_YAML in ${PREFIX}/environments/VEBA*.yml; do 
    ENV_NAME = $(basename $ENV_YAML .yml)
    echo "Creating ${ENV_NAME} environment"
    conda create -n $ENV_NAME -f $ENV_YAML -y
    cp ${PREFiX}/../src/*.py ${CONDA_BASE}/${ENV_NAME}/bin/
    cp ${PREFiX}/../src/scripts/ ${CONDA_BASE}/${ENV_NAME}/bin/
    ln -sf ${PREFiX}/../src/scripts/* ${CONDA_BASE}/${ENV_NAME}/bin/
    done 

# Database
DATABASE_DIRECTORY=$1
if [ -n "$DATABASE_DIRECTORY" ]; then
    echo "VEBA database will be in the following directory: $DATABASE_DIRECTORY"
else
    echo "Please specify a VEBA database directory"
    exit 1
fi

bash ${PREFIX}/install_database.sh $DATABASE_DIRECTORY




