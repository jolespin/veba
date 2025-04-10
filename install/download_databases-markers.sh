#!/bin/bash
# __version__ = "2025.4.3"
# VEBA_DATABASE_VERSION = "VEBA-DB_v9"
# MICROEUKAYROTIC_DATABASE_VERSION = "MicroEuk_v3"
# usage: bash veba/download_databases-preprocess.sh /path/to/veba_database_destination/

# Create database
DATABASE_DIRECTORY=${1:-"."}
REALPATH_DATABASE_DIRECTORY=$(realpath $DATABASE_DIRECTORY)
SCRIPT_DIRECTORY=$(dirname "$0")

MAXIMUM_NUMBER_OF_CPU=$(python -c "from multiprocessing import cpu_count; print(cpu_count())")
N_JOBS=${3:-${MAXIMUM_NUMBER_OF_CPU}}

# Database structure
echo ". .. ... ..... ........ ............."
echo "Creating directories"
echo ". .. ... ..... ........ ............."

mkdir -vp $DATABASE_DIRECTORY
mkdir -vp ${DATABASE_DIRECTORY}/Annotate
mkdir -vp ${DATABASE_DIRECTORY}/Classify
# >${DATABASE_DIRECTORY}/log

# Versions
DATE=$(date)
echo $DATE > ${DATABASE_DIRECTORY}/ACCESS_DATE

# MarkerSets
echo ". .. ... ..... ........ ............."
echo " * Processing profile HMM marker sets"
echo ". .. ... ..... ........ ............."
wget -v -O ${DATABASE_DIRECTORY}/MarkerSets.tar.gz https://figshare.com/ndownloader/files/36201486
tar xvzf ${DATABASE_DIRECTORY}/MarkerSets.tar.gz -C ${DATABASE_DIRECTORY}
gzip ${DATABASE_DIRECTORY}/MarkerSets/*.hmm
rm -rf ${DATABASE_DIRECTORY}/MarkerSets.tar.gz

echo -e " _    _ _______ ______  _______\n  \  /  |______ |_____] |_____|\n   \/   |______ |_____] |     |"
echo -e "..................................................."
echo -e "     (markers) Database Configuration Complete     "
echo -e "..................................................."
echo -e "[partial-database] If you are only installing individual databases you will need to set the following environment variable manually or use `update_environment_variables.sh` script:"
echo -e "[partial-database] \tVEBA_DATABASE=${REALPATH_DATABASE_DIRECTORY}"
echo -e "[partial-database] For walkthroughs on different workflows, please refer to the documentation:"
echo -e "[partial-database] \thttps://github.com/jolespin/veba/blob/main/walkthroughs/README.md"