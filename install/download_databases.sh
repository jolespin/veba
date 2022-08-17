#!/bin/bash

# Create database
DATABASE_DIRECTORY=${1:-"."}

# Database structure
echo ". .. ... ..... ........ ............."
echo "Creating directories"
echo ". .. ... ..... ........ ............."

mkdir -vp $DATABASE_DIRECTORY
mkdir -vp ${DATABASE_DIRECTORY}/Annotate
mkdir -vp ${DATABASE_DIRECTORY}/Classify
# mkdir -vp ${DATABASE_DIRECTORY}/MarkerSets
# >${DATABASE_DIRECTORY}/log

# Versions
DATE=$(date)
echo $DATE > ${DATABASE_DIRECTORY}/ACCESS_DATE

# NCBI Taxonomy
echo ". .. ... ..... ........ ............."
echo "i * Processing NCBITaxonomy"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy
wget -v -P ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz 
wget -v -P ${DATABASE_DIRECTORY} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
python -c 'import sys; from ete3 import NCBITaxa; NCBITaxa(taxdump_file="%s/taxdump.tar.gz"%(sys.argv[1]), dbfile="%s/Classify/NCBITaxonomy/taxa.sqlite"%(sys.argv[1]))' $DATABASE_DIRECTORY
tar xzfv ${DATABASE_DIRECTORY}/taxdump.tar.gz -C ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/
rm -rf ${DATABASE_DIRECTORY}/taxdump.tar.gz

# GTDB-Tk
echo ". .. ... ..... ........ ............."
echo "ii * Processing GTDB-Tk" # Need to update VEBA to use GTDB-Tk v2.1.0
echo ". .. ... ..... ........ ............."
# mkdir -v -p ${DATABASE_DIRECTORY}/Classify/gtdbtk
# wget -v -P ${DATABASE_DIRECTORY} https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
wget -v -P ${DATABASE_DIRECTORY} https://data.gtdb.ecogenomic.org/releases/release202/202.0/auxillary_files/gtdbtk_r202_data.tar.gz
tar xvzf ${DATABASE_DIRECTORY}/gtdbtk_r202_data.tar.gz -C ${DATABASE_DIRECTORY}
mv ${DATABASE_DIRECTORY}/release202 ${DATABASE_DIRECTORY}/Classify/GTDBTk
rm -rf ${DATABASE_DIRECTORY}/gtdbtk_r202_data.tar.gz

# CheckV
echo ". .. ... ..... ........ ............."
echo "iii * Processing CheckV"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Classify/CheckV
wget -v -P ${DATABASE_DIRECTORY} https://portal.nersc.gov/CheckV/checkv-db-v1.0.tar.gz
tar xvzf ${DATABASE_DIRECTORY}/checkv-db-v1.0.tar.gz -C ${DATABASE_DIRECTORY}
mv ${DATABASE_DIRECTORY}/checkv-db-v1.0 ${DATABASE_DIRECTORY}/Classify/CheckV
rm -rf ${DATABASE_DIRECTORY}/checkv-db-v1.0.tar.gz

# CheckM
echo ". .. ... ..... ........ ............."
echo "iv * Processing CheckM"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Classify/CheckM
wget -v -P ${DATABASE_DIRECTORY} https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar xvzf ${DATABASE_DIRECTORY}/checkm_data_2015_01_16.tar.gz -C ${DATABASE_DIRECTORY}/Classify/CheckM
rm -rf ${DATABASE_DIRECTORY}/checkm_data_2015_01_16.tar.gz

# Microeukaryotic 
echo ". .. ... ..... ........ ............."
echo "v * Processing Microeukaryotic MMSEQS2 database"
echo ". .. ... ..... ........ ............."
# Download from FigShare
wget -v -O Microeukaryotic.tar.gz https://figshare.com/ndownloader/files/34929255
tar xvzf Microeukaryotic.tar.gz -C ${DATABASE_DIRECTORY}/Classify
mmseqs createdb ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.rmdup.iupac.relabeled.no_deprecated.complete_lineage.faa.gz ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/microeukaryotic
rm -rf ${DATABASE_DIRECTORY}/Microeukaryotic.tar.gz

# MarkerSets
echo ". .. ... ..... ........ ............."
echo "vi * Processing profile HMM marker sets"
echo ". .. ... ..... ........ ............."
wget -v -O MarkerSets.tar.gz https://figshare.com/ndownloader/files/36201486
tar xvzf MarkerSets.tar.gz -C ${DATABASE_DIRECTORY}
rm -rf MarkerSets.tar.gz

# KOFAMSCAN
echo ". .. ... ..... ........ ............."
echo "vii * Processing KOFAMSCAN profile HMM marker sets"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Annotate/KOFAM
wget -v -O - ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz | gzip -d > ${DATABASE_DIRECTORY}/Annotate/KOFAM/ko_list
wget -v -c ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz -O - |  tar -xz
mv profiles ${DATABASE_DIRECTORY}/Annotate/KOFAM/

# Pfam
echo ". .. ... ..... ........ ............."
echo "viii * Processing Pfam profile HMM marker sets"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Annotate/Pfam
wget -v -P ${DATABASE_DIRECTORY}/Annotate/Pfam http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

# NCBI non-redundant
echo ". .. ... ..... ........ ............."
echo "ix * Processing NCBI non-redundant diamond database"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Annotate/nr
wget -v -P ${DATABASE_DIRECTORY} https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
diamond makedb --in ${DATABASE_DIRECTORY}/nr.gz --db ${DATABASE_DIRECTORY}/Annotate/nr.dmnd --taxonmap ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/prot.accession2taxid.FULL.gz --taxonnodes ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/nodes.dmp --taxonnames ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/names.dmp
rm -rf ${DATABASE_DIRECTORY}/nr.gz

# Contamination
echo ". .. ... ..... ........ ............."
echo "ix * Processing contamination databases"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Contamination
mkdir -v -p ${DATABASE_DIRECTORY}/Contamination/grch38
# mkdir -v -p ${DATABASE_DIRECTORY}/Contamination/antifam
mkdir -v -p ${DATABASE_DIRECTORY}/Contamination/kmers
wget -v -O ${DATABASE_DIRECTORY}/Contamination/kmers/ribokmers.fa.gz https://figshare.com/ndownloader/files/36220587
wget -v -P ${DATABASE_DIRECTORY} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
tar xvzf ${DATABASE_DIRECTORY}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz -C ${DATABASE_DIRECTORY}/Contamination/grch38

echo ". .. ... ..... ........ ............."
echo "xi * Adding the following environment variable to VEBA environments: export VEBA_DATABASE=${1}"
# CONDA_BASE=$(which conda | python -c "import sys; print('/'.join(sys.stdin.read().split('/')[:-2]))")
CONDA_BASE=$(conda run -n base bash -c "echo \${CONDA_PREFIX}")

# VEBA (Note, this will add extra lines if this is run more than once.  Might need to grep -v to another file then overwrite)
REALPATH_DATABASE_DIRECTORY=$(realpath $DATABASE_DIRECTORY)
for ENV_PREFIX in ${CONDA_BASE}/envs/VEBA-*; do 
    echo $ENV_PREFIX;
    mkdir -v -p ${ENV_PREFIX}/etc/conda/activate.d/
    mkdir -v -p ${ENV_PREFIX}/etc/conda/deactivate.d/
    echo "export VEBA_DATABASE=${REALPATH_DATABASE_DIRECTORY}" > ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset VEBA_DATABASE" > ${ENV_PREFIX}/etc/conda/deactivate.d/veba.sh
    done

#GTDB-Tk
echo ". .. ... ..... ........ ............."
echo "xii * Adding the following environment variable to VEBA environments: export GTDBTK_DATA_PATH=${REALPATH_DATABASE_DIRECTORY}/Classify/GTDBTk/"
for ENV_NAME in VEBA-binning-prokaryotic_env VEBA-classify_env; do 
    ENV_PREFIX=${CONDA_BASE}/envs/${ENV_NAME}
    # GTDBTK_DATABASE_VERSION=$(ls ${REALPATH_DATABASE_DIRECTORY}/Classify/GTDBTk)
    echo "export GTDBTK_DATA_PATH=${REALPATH_DATABASE_DIRECTORY}/Classify/GTDBTk/" >> ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset GTDBTK_DATA_PATH" >> ${ENV_PREFIX}/etc/conda/deactivate.d/veba.sh
    done 

# CheckV
echo ". .. ... ..... ........ ............."
echo "xiii * Adding the following environment variable to VEBA environments: export CHECKVDB=${REALPATH_DATABASE_DIRECTORY}/Classify/CheckV}"
for ENV_NAME in VEBA-binning-viral_env; do 
    ENV_PREFIX=${CONDA_BASE}/envs/${ENV_NAME}
    echo "export CHECKVDB=${REALPATH_DATABASE_DIRECTORY}/Classify/CheckV" >> ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset CHECKVDB" >> ${ENV_PREFIX}/etc/conda/deactivate.d/veba.sh
    done