#!/bin/bash
# __version__ = "2024.6.5"
# VEBA_DATABASE_VERSION = "VDB_v7"
# MICROEUKAYROTIC_DATABASE_VERSION = "MicroEuk_v3"
# usage: bash veba/download_databases-classify.sh /path/to/veba_database_destination/

# Create database
DATABASE_DIRECTORY=${1:-"."}
REALPATH_DATABASE_DIRECTORY=$(realpath $DATABASE_DIRECTORY)
SCRIPT_DIRECTORY=$(dirname "$0")

# N_JOBS=$(2:-"1")

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

# NCBI Taxonomy
echo ". .. ... ..... ........ ............."
echo " * Processing NCBITaxonomy"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy
# wget -v -P ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz 
wget -v -P ${DATABASE_DIRECTORY} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
# python -c 'import sys; from ete3 import NCBITaxa; NCBITaxa(taxdump_file="%s/taxdump.tar.gz"%(sys.argv[1]), dbfile="%s/Classify/NCBITaxonomy/taxa.sqlite"%(sys.argv[1]))' $DATABASE_DIRECTORY
tar xzfv ${DATABASE_DIRECTORY}/taxdump.tar.gz -C ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/
rm -rf ${DATABASE_DIRECTORY}/taxdump.tar.gz

# GTDB-Tk
echo ". .. ... ..... ........ ............."
echo " * Processing GTDB-Tk"
echo ". .. ... ..... ........ ............."

# # GTDB r220 (For future VEBA ≥ 2.2.0, VDB_7)
# Download from the data.ace.uq.edu.au mirror b/c it's way faster than data.gtdb.ecogenomic.org.
GTDB_VERSION="220"
wget -v -P ${DATABASE_DIRECTORY} https://data.ace.uq.edu.au/public/gtdb/data/releases/release${GTDB_VERSION}/${GTDB_VERSION}.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r${GTDB_VERSION}_data.tar.gz
tar xvzf ${DATABASE_DIRECTORY}/gtdbtk_r${GTDB_VERSION}_data.tar.gz -C ${DATABASE_DIRECTORY}
mv ${DATABASE_DIRECTORY}/release${GTDB_VERSION} ${DATABASE_DIRECTORY}/Classify/GTDB
echo "r${GTDB_VERSION}" > ${DATABASE_DIRECTORY}/Classify/GTDB/database_version
wget -P ${DATABASE_DIRECTORY}/Classify/GTDB/ https://data.ace.uq.edu.au/public/gtdb/data/releases/release${GTDB_VERSION}/${GTDB_VERSION}.0/RELEASE_NOTES.txt 
rm -rf ${DATABASE_DIRECTORY}/gtdbtk_r${GTDB_VERSION}_data.tar.gz

# GTDB r220 mash sketch
GTDB_ZENODO_RECORD_ID="11494307"
wget -v -O ${DATABASE_DIRECTORY}/gtdb_r${GTDB_VERSION}.msh https://zenodo.org/records/${GTDB_ZENODO_RECORD_ID}/files/gtdb_r${GTDB_VERSION}.msh?download=1
mkdir -p ${DATABASE_DIRECTORY}/Classify/GTDB/mash/
mv ${DATABASE_DIRECTORY}/gtdb_r${GTDB_VERSION}.msh ${DATABASE_DIRECTORY}/Classify/GTDB/mash/gtdb.msh

# # GTDB r214.1
# # data.gtdb.ecogenomic.org is slower than the faster mirror (data.ace.uq.edu.au/)
# # wget -v -P ${DATABASE_DIRECTORY} https://data.gtdb.ecogenomic.org/releases/release214/214.1/auxillary_files/gtdbtk_r214_data.tar.gz
# wget -v -P ${DATABASE_DIRECTORY} https://data.ace.uq.edu.au/public/gtdb/data/releases/release214/214.1/auxillary_files/gtdbtk_r214_data.tar.gz
# tar xvzf ${DATABASE_DIRECTORY}/gtdbtk_r214_data.tar.gz -C ${DATABASE_DIRECTORY}
# mv ${DATABASE_DIRECTORY}/release214 ${DATABASE_DIRECTORY}/Classify/GTDB
# rm -rf ${DATABASE_DIRECTORY}/gtdbtk_r214_data.tar.gz

# # GTDB r214.1 mash sketch
# wget -v -O ${DATABASE_DIRECTORY}/gtdb_r214.msh https://zenodo.org/record/8048187/files/gtdb_r214.msh?download=1
# mkdir -p ${DATABASE_DIRECTORY}/Classify/GTDB/mash/
# mv ${DATABASE_DIRECTORY}/gtdb_r214.msh ${DATABASE_DIRECTORY}/Classify/GTDB/mash/

# CheckV
echo ". .. ... ..... ........ ............."
echo " * Processing CheckV"
echo ". .. ... ..... ........ ............."
rm -rf ${DATABASE_DIRECTORY}/Classify/CheckV
CHECKVDB_VERSION="v1.5"
wget -v -P ${DATABASE_DIRECTORY} https://portal.nersc.gov/CheckV/checkv-db-${CHECKVDB_VERSION}.tar.gz
tar xvzf ${DATABASE_DIRECTORY}/checkv-db-${CHECKVDB_VERSION}.tar.gz -C ${DATABASE_DIRECTORY}
mv ${DATABASE_DIRECTORY}/checkv-db-${CHECKVDB_VERSION} ${DATABASE_DIRECTORY}/Classify/CheckV
echo "${CHECKV_VERSION}" > ${DATABASE_DIRECTORY}/Classify/CheckV/database_version
diamond makedb --in ${DATABASE_DIRECTORY}/Classify/CheckV/genome_db/checkv_reps.faa --db ${DATABASE_DIRECTORY}/Classify/CheckV/genome_db/checkv_reps.dmnd
rm -rf ${DATABASE_DIRECTORY}/checkv-db-${CHECKVDB_VERSION}.tar.gz

# geNomad
mkdir -p ${DATABASE_DIRECTORY}/Classify/geNomad
GENOMADDB_VERSION="v1.2"
wget -v -O ${DATABASE_DIRECTORY}/genomad_db_${GENOMADDB_VERSION}}.tar.gz https://zenodo.org/record/7586412/files/genomad_db_${GENOMADDB_VERSION}.tar.gz?download=1
tar xvzf ${DATABASE_DIRECTORY}/genomad_db_${GENOMADDB_VERSION}.tar.gz -C ${DATABASE_DIRECTORY}/Classify/geNomad --strip-components=1
rm -rf ${DATABASE_DIRECTORY}/genomad_db_${GENOMADDB_VERSION}.tar.gz

# CheckM2
echo ". .. ... ..... ........ ............."
echo " * Processing CheckM2"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Classify/CheckM2
wget -v -P ${DATABASE_DIRECTORY} https://zenodo.org/api/files/fd3bc532-cd84-4907-b078-2e05a1e46803/checkm2_database.tar.gz
tar xzfv ${DATABASE_DIRECTORY}/checkm2_database.tar.gz -C ${DATABASE_DIRECTORY}/Classify/CheckM2 --strip-components=1
rm -rf ${DATABASE_DIRECTORY}/checkm2_database.tar.gz

# Microeukaryotic 
echo ". .. ... ..... ........ ............."
echo " * Processing Microeukaryotic MMSEQS2 database"
echo ". .. ... ..... ........ ............."

## Download v2.1 from Zenodo
# wget -v -O ${DATABASE_DIRECTORY}/Microeukaryotic.tar.gz https://zenodo.org/record/7485114/files/VDB-Microeukaryotic_v2.tar.gz?download=1
# mkdir -p ${DATABASE_DIRECTORY}/Classify/Microeukaryotic && tar -xvzf ${DATABASE_DIRECTORY}/Microeukaryotic.tar.gz -C ${DATABASE_DIRECTORY}/Classify/Microeukaryotic --strip-components=1
# mmseqs createdb --compressed 1 ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.faa.gz ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/microeukaryotic
# rm -rf ${DATABASE_DIRECTORY}/Microeukaryotic.tar.gz

# # eukaryota_odb10 subset of Microeukaryotic Protein Database
# wget -v -O ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.eukaryota_odb10.list https://zenodo.org/record/7485114/files/reference.eukaryota_odb10.list?download=1
# seqkit grep -f ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.eukaryota_odb10.list ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.faa.gz > ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.eukaryota_odb10.faa
# mmseqs createdb --compressed 1 ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.eukaryota_odb10.faa ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/microeukaryotic.eukaryota_odb10
# rm -rf ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.eukaryota_odb10.faa
# rm -rf ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.faa.gz # Comment this out if you want to keep the actual protein sequences

# Download MicroEuk_v3 from Zenodo
wget -v -O ${DATABASE_DIRECTORY}/MicroEuk_v3.tar.gz https://zenodo.org/records/10139451/files/MicroEuk_v3.tar.gz?download=1 
tar xvzf ${DATABASE_DIRECTORY}/MicroEuk_v3.tar.gz -C ${DATABASE_DIRECTORY}
mkdir -p ${DATABASE_DIRECTORY}/Classify/MicroEuk

# Source Taxonomy
cp -rf ${DATABASE_DIRECTORY}/MicroEuk_v3/source_taxonomy.tsv.gz ${DATABASE_DIRECTORY}/Classify/MicroEuk

# MicroEuk100
gzip -d ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk100.faa.gz
mmseqs createdb --compressed 1 ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk100.faa ${DATABASE_DIRECTORY}/Classify/MicroEuk/MicroEuk100

# MicroEuk100.eukaryota_odb10
gzip -d ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk100.eukaryota_odb10.list.gz
seqkit grep -f ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk100.eukaryota_odb10.list ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk100.faa | mmseqs createdb --compressed 1 stdin ${DATABASE_DIRECTORY}/Classify/MicroEuk/MicroEuk100.eukaryota_odb10

# MicroEuk90
gzip -d -c ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk90_clusters.tsv.gz | cut -f1 | sort -u > ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk90.list
seqkit grep -f ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk90.list ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk100.faa | mmseqs createdb --compressed 1 stdin ${DATABASE_DIRECTORY}/Classify/MicroEuk/MicroEuk90

# MicroEuk50
gzip -d -c ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk50_clusters.tsv.gz | cut -f1 | sort -u > ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk50.list
seqkit grep -f ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk50.list ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk100.faa | mmseqs createdb --compressed 1 stdin ${DATABASE_DIRECTORY}/Classify/MicroEuk/MicroEuk50

# source_to_lineage.dict.pkl.gz
build_source_to_lineage_dictionary.py -i ${DATABASE_DIRECTORY}/MicroEuk_v3/source_taxonomy.tsv.gz -o ${DATABASE_DIRECTORY}/Classify/MicroEuk/source_to_lineage.dict.pkl.gz

# target_to_source.dict.pkl.gz
build_target_to_source_dictionary.py -i ${DATABASE_DIRECTORY}/MicroEuk_v3/identifier_mapping.proteins.tsv.gz -o ${DATABASE_DIRECTORY}/Classify/MicroEuk/target_to_source.dict.pkl.gz

# Remove intermediate files
rm -rf ${DATABASE_DIRECTORY}/MicroEuk_v3/
rm -rf ${DATABASE_DIRECTORY}/MicroEuk_v3.tar.gz

echo -e " _    _ _______ ______  _______\n  \  /  |______ |_____] |_____|\n   \/   |______ |_____] |     |"
echo -e "...................................................."
echo -e "     (classify) Database Configuration Complete     "
echo -e "...................................................."
echo -e "[partial-database] If you are only installing individual databases you will need to set the following environment variable manually or use `update_environment_variables.sh` script:"
echo -e "[partial-database] \tVEBA_DATABASE=${REALPATH_DATABASE_DIRECTORY}"
echo -e "[partial-database] For walkthroughs on different workflows, please refer to the documentation:"
echo -e "[partial-database] \thttps://github.com/jolespin/veba/blob/main/walkthroughs/README.md"