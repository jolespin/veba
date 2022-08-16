#!/bin/bash
SCRIPT_PATH=$(realpath $0)
PREFIX=$(echo $SCRIPT_PATH | python -c "import sys; print('/'.join(sys.stdin.read().split('/')[:-1]))")
CONDA_BASE=$(conda run -n base bash -c "echo \${CONDA_PREFIX}")

# Update permissions
echo "Updating permissions for scripts in ${PREFIX}/../src"
chmod 755 ${PREFIX}/../src/*.py
chmod 755 ${PREFIX}/../src/scripts/*

# Environemnts
for ENV_YAML in ${PREFIX}/environments/VEBA*.yml; do
    # Get environment name
    ENV_NAME=$(basename $ENV_YAML .yml)

    # Create conda environment
    echo "Creating ${ENV_NAME} environment"
    conda env create -n $ENV_NAME -f $ENV_YAML

    # Copy over files to environment bin/
    cp -r ${PREFIX}/../src/*.py ${CONDA_BASE}/envs/${ENV_NAME}/bin/
    cp -r ${PREFIX}/../src/scripts/ ${CONDA_BASE}/envs/${ENV_NAME}/bin/

    # Symlink the accessory scripts to bin/
    ln -sf ${PREFIX}/../src/scripts/* ${CONDA_BASE}/envs/${ENV_NAME}/bin/

    done

echo -e " _    _ _______ ______  _______\n  \  /  |______ |_____] |_____|\n   \/   |______ |_____] |     |"
echo -e "..............................."
echo -e "     Installation Complete     "
echo -e "..............................."
echo -e "Please run 'download_databases.sh' script available in the installation directory. If you need to redownload:"
echo -e "\twget https://raw.githubusercontent.com/jolespin/veba/main/install/download_databases.sh"
echo -e "For help or instructions, refer to the installation walkthrough: \n\thttps://github.com/jolespin/veba/blob/main/install/README.md"