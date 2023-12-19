#!/bin/bash
# __version__ = "2023.12.19"

SCRIPT_PATH=$(realpath $0)
PREFIX=$(echo $SCRIPT_PATH | python -c "import sys; print('/'.join(sys.stdin.read().split('/')[:-1]))")
# CONDA_BASE=$(conda run -n base bash -c "echo \${CONDA_PREFIX}")
CONDA_BASE=$(conda info --base)

# Update permissions
echo "Updating permissions for scripts in ${PREFIX}/../src"
chmod 755 ${PREFIX}/../src/veba
chmod 755 ${PREFIX}/../src/*.py
chmod 755 ${PREFIX}/../src/scripts/*

# Install mamba
conda install -c conda-forge mamba -y 
# conda update mamba -y  # Recommended

# Environemnts
# Main environment
echo "Creating ${VEBA} main environment"

ENV_NAME="VEBA"
mamba create -y -n $ENV_NAME -c conda-forge -c bioconda -c jolespin seqkit genopype networkx biopython biom-format anndata || (echo "Error when creating main VEBA environment" ; exit 1) &> ${PREFIX}/environments/VEBA.log

# Copy main executable
echo -e "\t*Copying main VEBA executable into ${ENV_NAME} environment path"
cp -r ${PREFIX}/../src/veba ${CONDA_BASE}/envs/${ENV_NAME}/bin/
# Copy over files to environment bin/
echo -e "\t*Copying VEBA modules into ${ENV_NAME} environment path"
cp -r ${PREFIX}/../src/*.py ${CONDA_BASE}/envs/${ENV_NAME}/bin/
echo -e "\t*Copying VEBA utility scripts into ${ENV_NAME} environment path"
cp -r ${PREFIX}/../src/scripts/ ${CONDA_BASE}/envs/${ENV_NAME}/bin/
# Symlink the utility scripts to bin/
echo -e "\t*Symlinking VEBA utility scripts into ${ENV_NAME} environment path"
ln -sf ${CONDA_BASE}/envs/${ENV_NAME}/bin/scripts/* ${CONDA_BASE}/envs/${ENV_NAME}/bin/

# Version
cp -rf ${PREFIX}/../VERSION ${CONDA_BASE}/envs/${ENV_NAME}/bin/VEBA_VERSION

# Module environments
for ENV_YAML in ${PREFIX}/environments/VEBA*.yml; do
    # Get environment name
    ENV_NAME=$(basename $ENV_YAML .yml)

    # Create conda environment
    echo "Creating ${ENV_NAME} module environment"
    mamba env create -n $ENV_NAME -f $ENV_YAML || (echo "Error when creating VEBA environment: ${ENV_YAML}" ; exit 1) &> ${ENV_YAML}.log

    # Copy over files to environment bin/
    echo -e "\t*Copying VEBA modules into ${ENV_NAME} environment path"
    cp -r ${PREFIX}/../src/*.py ${CONDA_BASE}/envs/${ENV_NAME}/bin/
    echo -e "\t*Copying VEBA utility scripts into ${ENV_NAME} environment path"
    cp -r ${PREFIX}/../src/scripts/ ${CONDA_BASE}/envs/${ENV_NAME}/bin/
    # Symlink the utility scripts to bin/
    echo -e "\t*Symlinking VEBA utility scripts into ${ENV_NAME} environment path"
    ln -sf ${CONDA_BASE}/envs/${ENV_NAME}/bin/scripts/* ${CONDA_BASE}/envs/${ENV_NAME}/bin/

    # Version
    cp -rf ${PREFIX}/../VERSION ${CONDA_BASE}/envs/${ENV_NAME}/bin/VEBA_VERSION

    done

echo -e " _    _ _______ ______  _______\n  \  /  |______ |_____] |_____|\n   \/   |______ |_____] |     |"
echo -e "..............................."
echo -e "     Installation Complete     "
echo -e "..............................."
echo -e "Please run 'download_databases.sh' script available in the installation directory. If you need to redownload:"
echo -e "\twget https://raw.githubusercontent.com/jolespin/veba/main/install/download_databases.sh"
echo -e "For help or instructions, refer to the installation walkthrough: \n\thttps://github.com/jolespin/veba/blob/main/install/README.md"