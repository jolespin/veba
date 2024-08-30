#!/bin/bash
# __version__ = "2024.8.30"
# VEBA_DATABASE_VERSION = "VDB_v7"
# MICROEUKAYROTIC_DATABASE_VERSION = "MicroEuk_v3"
# usage: bash veba/download_databases-contamination.sh /path/to/veba_database_destination/

# Create database
DATABASE_DIRECTORY=${1:-"."}
REALPATH_DATABASE_DIRECTORY=$(realpath $DATABASE_DIRECTORY)
SCRIPT_DIRECTORY=$(dirname "$0")

MAXIMUM_NUMBER_OF_CPU=$(python -c "from multiprocessing import cpu_count; print(cpu_count())")
N_JOBS=$(2:-${MAXIMUM_NUMBER_OF_CPU})

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

# Contamination
echo ". .. ... ..... ........ ............."
echo " * Processing contamination databases"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Contamination

# Ribokmers
mkdir -v -p ${DATABASE_DIRECTORY}/Contamination/kmers
wget -v -O ${DATABASE_DIRECTORY}/Contamination/kmers/ribokmers.fa.gz https://figshare.com/ndownloader/files/36220587

# T2T-CHM13v2.0
# Human (Bowtie2 Index)
wget -v -P ${DATABASE_DIRECTORY} https://genome-idx.s3.amazonaws.com/bt/chm13v2.0.zip
unzip -d ${DATABASE_DIRECTORY}/Contamination/ ${DATABASE_DIRECTORY}/chm13v2.0.zip
rm -rf ${DATABASE_DIRECTORY}/chm13v2.0.zip

# # Human (MiniMap2 Index) (Uncomment if you plan on using long reads (7.1 GB))
# wget -v -P ${DATABASE_DIRECTORY} https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
# minimap2 -d ${DATABASE_DIRECTORY}/Contamination/chm13v2.0/chm13v2.0.mmi ${DATABASE_DIRECTORY}/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
# rm -rf ${DATABASE_DIRECTORY}/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz


# # Mouse (Bowtie2 Index) (Uncomment if you plan on using mouse microbiomes)
# wget -v -P ${DATABASE_DIRECTORY} https://genome-idx.s3.amazonaws.com/bt/GRCm39.zip
# unzip -d ${DATABASE_DIRECTORY}/Contamination/ ${DATABASE_DIRECTORY}/GRCm39.zip
# rm -rf ${DATABASE_DIRECTORY}/GRCm39.zip

# Contamination
# AntiFam
mkdir -v -p ${DATABASE_DIRECTORY}/Contamination/AntiFam
wget -v -O ${DATABASE_DIRECTORY}/Antifam.tar.gz https://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz
tar xzfv ${DATABASE_DIRECTORY}/Antifam.tar.gz -C ${DATABASE_DIRECTORY}/Contamination/AntiFam 
rm ${DATABASE_DIRECTORY}/Contamination/AntiFam/AntiFam_*.hmm
gzip ${DATABASE_DIRECTORY}/Contamination/AntiFam/*.hmm
rm -rf ${DATABASE_DIRECTORY}/Antifam.tar.gz
rm -rf ${DATABASE_DIRECTORY}/Contamination/AntiFam/*.seed


echo -e " _    _ _______ ______  _______\n  \  /  |______ |_____] |_____|\n   \/   |______ |_____] |     |"
echo -e "......................................................"
echo -e "     (contamination) Database Configuration Complete     "
echo -e "......................................................"
echo -e "[partial-database] If you are only installing individual databases you will need to set the following environment variable manually or use `update_environment_variables.sh` script:"
echo -e "[partial-database] \tVEBA_DATABASE=${REALPATH_DATABASE_DIRECTORY}"
echo -e "[partial-database] For walkthroughs on different workflows, please refer to the documentation:"
echo -e "[partial-database] \thttps://github.com/jolespin/veba/blob/main/walkthroughs/README.md"