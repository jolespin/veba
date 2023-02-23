#!/bin/bash
# __VERSION__ = "2023.1.24"
SCRIPT_PATH=$(realpath $0)
PREFIX=$(echo $SCRIPT_PATH | python -c "import sys; print('/'.join(sys.stdin.read().split('/')[:-1]))")
CONDA_BASE=$(conda run -n base bash -c "echo \${CONDA_PREFIX}")

# Update permissions
echo "Updating permissions for scripts in ${PREFIX}/../src"
chmod 755 ${PREFIX}/../src/*.py
chmod 755 ${PREFIX}/../src/scripts/*

# Install mamba
conda install -c conda-forge mamba -y 
# conda update mamba -y  # Recommended

# Environemnts
for ENV_YAML in ${PREFIX}/environments/VEBA*.yml; do
    # Get environment name
    ENV_NAME=$(basename $ENV_YAML .yml)

    # Create conda environment
    echo "Creating ${ENV_NAME} environment"
    time(mamba env create -n $ENV_NAME -f $ENV_YAML || (echo "Error when creating VEBA environment: ${ENV_YAML}" ; exit 1)) &> ${ENV_YAML}.log

    # Copy over files to environment bin/
    echo -e "\t*Copying VEBA modules into ${ENV_NAME} environment path"
    cp -r ${PREFIX}/../src/*.py ${CONDA_BASE}/envs/${ENV_NAME}/bin/
    echo -e "\t*Copying VEBA utility scripts into ${ENV_NAME} environment path"
    cp -r ${PREFIX}/../src/scripts/ ${CONDA_BASE}/envs/${ENV_NAME}/bin/
    # Symlink the utility scripts to bin/
    echo -e "\t*Symlinking VEBA utility scripts into ${ENV_NAME} environment path"
    ln -sf ${CONDA_BASE}/envs/${ENV_NAME}/bin/scripts/* ${CONDA_BASE}/envs/${ENV_NAME}/bin/

    done

echo -e " _    _ _______ ______  _______\n  \  /  |______ |_____] |_____|\n   \/   |______ |_____] |     |"
echo -e "..............................."
echo -e "     Installation Complete     "
echo -e "..............................."
echo -e "Please run 'download_databases.sh' script available in the installation directory. If you need to redownload:"
echo -e "\twget https://raw.githubusercontent.com/jolespin/veba/main/install/download_databases.sh"
echo -e "For help or instructions, refer to the installation walkthrough: \n\thttps://github.com/jolespin/veba/blob/main/install/README.md"