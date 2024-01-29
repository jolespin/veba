#!/bin/bash
# __version__ = "2024.1.23"

SCRIPT_PATH=$(realpath $0)
SCRIPT_DIRECTORY=$(dirname $0)
VEBA_REPOSITORY_DIRECTORY=$(realpath ${SCRIPT_DIRECTORY}/../)

LOG_DIRECTORY=${1:-"logs/veba_installation"}
mkdir -pv ${LOG_DIRECTORY}

# CONDA_BASE=$(conda run -n base bash -c "echo \${CONDA_PREFIX}")
CONDA_BASE=$(conda info --base)

# Update permissions
echo "Updating permissions for scripts in ${VEBA_REPOSITORY_DIRECTORY}/bin"
chmod 755 ${VEBA_REPOSITORY_DIRECTORY}/bin/veba
chmod 755 ${VEBA_REPOSITORY_DIRECTORY}/bin/*.py
chmod 755 ${VEBA_REPOSITORY_DIRECTORY}/bin/scripts/*

# Install mamba
conda install -c conda-forge mamba -y 
# conda update mamba -y  # Recommended

# Environemnts
# Main environment
echo "Creating ${VEBA} main environment"

ENV_NAME="VEBA"
mamba create -y -n $ENV_NAME -c conda-forge -c bioconda -c jolespin seqkit genopype networkx biopython biom-format anndata || (echo "Error when creating main VEBA environment" ; exit 1) &> ${LOG_DIRECTORY}/VEBA.log

# Copy main executable
echo -e "\t*Copying main VEBA executable into ${ENV_NAME} environment path"
cp -r ${VEBA_REPOSITORY_DIRECTORY}/bin/veba ${CONDA_BASE}/envs/${ENV_NAME}/bin/
# Copy over files to environment bin/
echo -e "\t*Copying VEBA modules into ${ENV_NAME} environment path"
cp -r ${VEBA_REPOSITORY_DIRECTORY}/bin/*.py ${CONDA_BASE}/envs/${ENV_NAME}/bin/
echo -e "\t*Copying VEBA utility scripts into ${ENV_NAME} environment path"
cp -r ${VEBA_REPOSITORY_DIRECTORY}/bin/scripts/ ${CONDA_BASE}/envs/${ENV_NAME}/bin/
# Symlink the utility scripts to bin/
echo -e "\t*Symlinking VEBA utility scripts into ${ENV_NAME} environment path"
ln -sf ${CONDA_BASE}/envs/${ENV_NAME}/bin/scripts/* ${CONDA_BASE}/envs/${ENV_NAME}/bin/

# Version
cp -rf ${VEBA_REPOSITORY_DIRECTORY}/VERSION ${CONDA_BASE}/envs/${ENV_NAME}/bin/VEBA_VERSION
cp -r ${VEBA_REPOSITORY_DIRECTORY}/bin/VEBA_SCRIPT_VERSIONS ${CONDA_BASE}/envs/${ENV_NAME}/bin/

# Module environments
for ENV_YAML in ${VEBA_REPOSITORY_DIRECTORY}/install/environments/VEBA*.yml; do
    # Get environment name
    ENV_NAME=$(basename $ENV_YAML .yml)

    # Create conda environment
    echo "Creating ${ENV_NAME} module environment"
    mamba env create -n $ENV_NAME -f $ENV_YAML || (echo "Error when creating VEBA environment: ${ENV_YAML}" ; exit 1) &> ${LOG_DIRECTORY}/${ENV_NAME}.log

    # Copy over files to environment bin/
    echo -e "\t*Copying VEBA modules into ${ENV_NAME} environment path"
    cp -r ${VEBA_REPOSITORY_DIRECTORY}//bin/*.py ${CONDA_BASE}/envs/${ENV_NAME}/bin/
    echo -e "\t*Copying VEBA utility scripts into ${ENV_NAME} environment path"
    cp -r ${VEBA_REPOSITORY_DIRECTORY}/bin/scripts/ ${CONDA_BASE}/envs/${ENV_NAME}/bin/
    # Symlink the utility scripts to bin/
    echo -e "\t*Symlinking VEBA utility scripts into ${ENV_NAME} environment path"
    ln -sf ${CONDA_BASE}/envs/${ENV_NAME}/bin/scripts/* ${CONDA_BASE}/envs/${ENV_NAME}/bin/

    # Version
    cp -rf ${VEBA_REPOSITORY_DIRECTORY}/VERSION ${CONDA_BASE}/envs/${ENV_NAME}/bin/VEBA_VERSION
    cp -r ${VEBA_REPOSITORY_DIRECTORY}/bin/VEBA_SCRIPT_VERSIONS ${CONDA_BASE}/envs/${ENV_NAME}/bin/

    done

echo -e " _    _ _______ ______  _______\n  \  /  |______ |_____] |_____|\n   \/   |______ |_____] |     |"
echo -e "..............................."
echo -e "     Installation Complete     "
echo -e "..............................."
echo -e "Please run 'download_databases.sh' script available in the installation directory. If you need to redownload:"
echo -e "\twget https://raw.githubusercontent.com/jolespin/veba/main/install/download_databases.sh"
echo -e "For help or instructions, refer to the installation walkthrough: \n\thttps://github.com/jolespin/veba/blob/main/install/README.md"