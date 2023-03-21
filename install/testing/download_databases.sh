#!/bin/bash
# __VERSION__ = "2023.3.6"
# VEBA_DATABASE_VERSION = "VDB_v4.1"
# MICROEUKAYROTIC_DATABASE_VERSION = "VDB-Microeukaryotic_v2.1"

# Create database
DATABASE_DIRECTORY=${1:-"."}
REALPATH_DATABASE_DIRECTORY=$(realpath $DATABASE_DIRECTORY)

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
echo "i * Processing NCBITaxonomy"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy
wget -v -P ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz 
wget -v -P ${DATABASE_DIRECTORY} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
# python -c 'import sys; from ete3 import NCBITaxa; NCBITaxa(taxdump_file="%s/taxdump.tar.gz"%(sys.argv[1]), dbfile="%s/Classify/NCBITaxonomy/taxa.sqlite"%(sys.argv[1]))' $DATABASE_DIRECTORY
tar xzfv ${DATABASE_DIRECTORY}/taxdump.tar.gz -C ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/
rm -rf ${DATABASE_DIRECTORY}/taxdump.tar.gz

# GTDB-Tk
echo ". .. ... ..... ........ ............."
echo "ii * Processing GTDB-Tk"
echo ". .. ... ..... ........ ............."

# For GTDBTk v2
wget -v -P ${DATABASE_DIRECTORY} https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz
tar xvzf ${DATABASE_DIRECTORY}/gtdbtk_r207_v2_data.tar.gz -C ${DATABASE_DIRECTORY}
mv ${DATABASE_DIRECTORY}/release207_v2 ${DATABASE_DIRECTORY}/Classify/GTDBTk
rm -rf ${DATABASE_DIRECTORY}/gtdbtk_r207_v2_data.tar.gz

# CheckV
echo ". .. ... ..... ........ ............."
echo "iii * Processing CheckV"
echo ". .. ... ..... ........ ............."
rm -rf ${DATABASE_DIRECTORY}/Classify/CheckV
wget -v -P ${DATABASE_DIRECTORY} https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz
tar xvzf ${DATABASE_DIRECTORY}/checkv-db-v1.5.tar.gz -C ${DATABASE_DIRECTORY}
mv ${DATABASE_DIRECTORY}/checkv-db-v1.5 ${DATABASE_DIRECTORY}/Classify/CheckV
diamond makedb --in ${DATABASE_DIRECTORY}/Classify/CheckV/genome_db/checkv_reps.faa --db ${DATABASE_DIRECTORY}/Classify/CheckV/genome_db/checkv_reps.dmnd
rm -rf ${DATABASE_DIRECTORY}/checkv-db-v1.5.tar.gz

# geNomad
mkdir -p ${DATABASE_DIRECTORY}/Classify/geNomad
wget -v -O ${DATABASE_DIRECTORY}/genomad_db_v1.2.tar.gz https://zenodo.org/record/7586412/files/genomad_db_v1.2.tar.gz?download=1
tar xvzf ${DATABASE_DIRECTORY}/genomad_db_v1.2.tar.gz -C ${DATABASE_DIRECTORY}/Classify/geNomad --strip-components=1
rm -rf ${DATABASE_DIRECTORY}/genomad_db_v1.2.tar.gz

# CheckM2
echo ". .. ... ..... ........ ............."
echo "iv * Processing CheckM2"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Classify/CheckM2
wget -v -P ${DATABASE_DIRECTORY} https://zenodo.org/api/files/fd3bc532-cd84-4907-b078-2e05a1e46803/checkm2_database.tar.gz
tar xzfv ${DATABASE_DIRECTORY}/checkm2_database.tar.gz -C ${DATABASE_DIRECTORY}/Classify/CheckM2 --strip-components=1
rm -rf ${DATABASE_DIRECTORY}/checkm2_database.tar.gz

# Microeukaryotic 
echo ". .. ... ..... ........ ............."
echo "v * Processing Microeukaryotic MMSEQS2 database"
echo ". .. ... ..... ........ ............."

# Download v2.1 from Zenodo
wget -v -O ${DATABASE_DIRECTORY}/Microeukaryotic.tar.gz https://zenodo.org/record/7485114/files/VDB-Microeukaryotic_v2.tar.gz?download=1
mkdir -p ${DATABASE_DIRECTORY}/Classify/Microeukaryotic && tar -xvzf ${DATABASE_DIRECTORY}/Microeukaryotic.tar.gz -C ${DATABASE_DIRECTORY}/Classify/Microeukaryotic --strip-components=1
mmseqs createdb ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.faa.gz ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/microeukaryotic
rm -rf ${DATABASE_DIRECTORY}/Microeukaryotic.tar.gz

# eukaryota_odb10 subset of Microeukaryotic Protein Database
wget -v -O ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.eukaryota_odb10.list https://zenodo.org/record/7485114/files/reference.eukaryota_odb10.list?download=1
seqkit grep -f ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.eukaryota_odb10.list ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.faa.gz > ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.eukaryota_odb10.faa
mmseqs createdb ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.eukaryota_odb10.faa ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/microeukaryotic.eukaryota_odb10
rm -rf ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.eukaryota_odb10.faa
rm -rf ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.faa.gz # Comment this out if you want to keep the actual protein sequences

# MarkerSets
echo ". .. ... ..... ........ ............."
echo "vi * Processing profile HMM marker sets"
echo ". .. ... ..... ........ ............."
wget -v -O ${DATABASE_DIRECTORY}/MarkerSets.tar.gz https://figshare.com/ndownloader/files/36201486
tar xvzf ${DATABASE_DIRECTORY}/MarkerSets.tar.gz -C ${DATABASE_DIRECTORY}
gzip ${DATABASE_DIRECTORY}/MarkerSets/*.hmm
rm -rf ${DATABASE_DIRECTORY}/MarkerSets.tar.gz

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
wget -v -P ${DATABASE_DIRECTORY}/Annotate/Pfam http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt

# NCBIfam-AMRFinder
echo ". .. ... ..... ........ ............."
echo "ix * Processing NCBIfam-AMRFinder profile HMM marker sets"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Annotate/NCBIfam-AMRFinder
wget -v -P ${DATABASE_DIRECTORY}/Annotate/NCBIfam-AMRFinder https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/latest/NCBIfam-AMRFinder.changelog.txt
wget -v -P ${DATABASE_DIRECTORY}/Annotate/NCBIfam-AMRFinder https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/latest/NCBIfam-AMRFinder.tsv
wget -v -O ${DATABASE_DIRECTORY}/NCBIfam-AMRFinder.HMM.tar.gz https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/latest/NCBIfam-AMRFinder.HMM.tar.gz
tar xzfv ${DATABASE_DIRECTORY}/NCBIfam-AMRFinder.HMM.tar.gz -C ${DATABASE_DIRECTORY}/Annotate/NCBIfam-AMRFinder --strip-components=1
cat ${DATABASE_DIRECTORY}/Annotate/NCBIfam-AMRFinder/*.HMM | gzip > ${DATABASE_DIRECTORY}/Annotate/NCBIfam-AMRFinder/NCBIfam-AMRFinder.hmm.gz
rm -rf ${DATABASE_DIRECTORY}/Annotate/NCBIfam-AMRFinder/*.HMM
rm -rf ${DATABASE_DIRECTORY}/NCBIfam-AMRFinder.HMM.tar.gz

# NCBI non-redundant
echo ". .. ... ..... ........ ............."
echo "x * Processing NCBI non-redundant diamond database"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Annotate/nr
wget -v -P ${DATABASE_DIRECTORY} https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
diamond makedb --in ${DATABASE_DIRECTORY}/nr.gz --db ${DATABASE_DIRECTORY}/Annotate/nr/nr.dmnd --taxonmap ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/prot.accession2taxid.FULL.gz --taxonnodes ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/nodes.dmp --taxonnames ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/names.dmp
rm -rf ${DATABASE_DIRECTORY}/nr.gz

#MiBIG
mkdir -v -p ${DATABASE_DIRECTORY}/Annotate/MIBiG
wget -v -P ${DATABASE_DIRECTORY} https://dl.secondarymetabolites.org/mibig/mibig_prot_seqs_3.1.fasta
seqkit rmdup -s ${DATABASE_DIRECTORY}/mibig_prot_seqs_3.1.fasta > ${DATABASE_DIRECTORY}/mibig_prot_seqs_3.1.rmdup.fasta
diamond makedb --in ${DATABASE_DIRECTORY}/mibig_prot_seqs_3.1.rmdup.fasta --db ${DATABASE_DIRECTORY}/Annotate/MIBiG/mibig_v3.1.dmnd
rm -rf ${DATABASE_DIRECTORY}/mibig_prot_seqs_3.1.fasta
rm -rf ${DATABASE_DIRECTORY}/mibig_prot_seqs_3.1.rmdup.fasta

# Contamination
echo ". .. ... ..... ........ ............."
echo "xi * Processing contamination databases"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Contamination

# AntiFam
mkdir -v -p ${DATABASE_DIRECTORY}/Contamination/AntiFam
wget -v -O ${DATABASE_DIRECTORY}/Antifam.tar.gz https://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz
tar xzfv ${DATABASE_DIRECTORY}/Antifam.tar.gz -C ${DATABASE_DIRECTORY}/Contamination/AntiFam 
rm ${DATABASE_DIRECTORY}/Contamination/AntiFam/AntiFam_*.hmm
gzip ${DATABASE_DIRECTORY}/Contamination/AntiFam/*.hmm
rm -rf ${DATABASE_DIRECTORY}/Antifam.tar.gz
rm -rf ${DATABASE_DIRECTORY}/Contamination/AntiFam/*.seed

# Ribokmers
mkdir -v -p ${DATABASE_DIRECTORY}/Contamination/kmers
wget -v -O ${DATABASE_DIRECTORY}/Contamination/kmers/ribokmers.fa.gz https://figshare.com/ndownloader/files/36220587

# Replacing GRCh38 with CHM13v2.0 in v2022.10.18
wget -v -P ${DATABASE_DIRECTORY} https://genome-idx.s3.amazonaws.com/bt/chm13v2.0.zip
unzip -d ${DATABASE_DIRECTORY}/Contamination/ ${DATABASE_DIRECTORY}/chm13v2.0.zip
rm -rf ${DATABASE_DIRECTORY}/chm13v2.0.zip

echo ". .. ... ..... ........ ............."
echo "xii * Adding the following environment variable to VEBA environments: export VEBA_DATABASE=${REALPATH_DATABASE_DIRECTORY}"
# CONDA_BASE=$(which conda | python -c "import sys; print('/'.join(sys.stdin.read().split('/')[:-2]))")
CONDA_BASE=$(conda run -n base bash -c "echo \${CONDA_PREFIX}")

# VEBA
for ENV_PREFIX in ${CONDA_BASE}/envs/test-VEBA-*; do 
    echo $ENV_PREFIX;
    mkdir -v -p ${ENV_PREFIX}/etc/conda/activate.d/
    mkdir -v -p ${ENV_PREFIX}/etc/conda/deactivate.d/
    echo "export VEBA_DATABASE=${REALPATH_DATABASE_DIRECTORY}" > ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset VEBA_DATABASE" > ${ENV_PREFIX}/etc/conda/deactivate.d/veba.sh
    done

#CheckM2
echo ". .. ... ..... ........ ............."
echo "xiii * Adding the following environment variable to VEBA environments: export CHECKM2DB=${REALPATH_DATABASE_DIRECTORY}/Classify/CheckM2/uniref100.KO.1.dmnd"
for ENV_NAME in VEBA-binning-prokaryotic_env; do 
    ENV_NAME=test-${ENV_NAME}
    ENV_PREFIX=${CONDA_BASE}/envs/${ENV_NAME}
    # CheckM2
    echo "export CHECKM2DB=${REALPATH_DATABASE_DIRECTORY}/Classify/CheckM2/uniref100.KO.1.dmnd" >> ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset CHECKM2DB" >> ${ENV_PREFIX}/etc/conda/deactivate.d/veba.sh    
    done 

#GTDB-Tk
echo ". .. ... ..... ........ ............."
echo "xiv * Adding the following environment variable to VEBA environments: export GTDBTK_DATA_PATH=${REALPATH_DATABASE_DIRECTORY}/Classify/GTDBTk/"
for ENV_NAME in VEBA-classify_env; do 
    ENV_NAME=test-${ENV_NAME}
    ENV_PREFIX=${CONDA_BASE}/envs/${ENV_NAME}
    # GTDB-Tk
    echo "export GTDBTK_DATA_PATH=${REALPATH_DATABASE_DIRECTORY}/Classify/GTDBTk/" >> ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset GTDBTK_DATA_PATH" >> ${ENV_PREFIX}/etc/conda/deactivate.d/veba.sh
    done 

# CheckV
echo ". .. ... ..... ........ ............."
echo "xv * Adding the following environment variable to VEBA environments: export CHECKVDB=${REALPATH_DATABASE_DIRECTORY}/Classify/CheckV/"
for ENV_NAME in VEBA-binning-viral_env; do 
    ENV_NAME=test-${ENV_NAME}
    ENV_PREFIX=${CONDA_BASE}/envs/${ENV_NAME}
    echo "export CHECKVDB=${REALPATH_DATABASE_DIRECTORY}/Classify/CheckV" >> ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset CHECKVDB" >> ${ENV_PREFIX}/etc/conda/deactivate.d/veba.sh
    done

echo -e " _    _ _______ ______  _______\n  \  /  |______ |_____] |_____|\n   \/   |______ |_____] |     |"
echo -e "........................................."
echo -e "     Database Configuration Complete     "
echo -e ".........................................."
echo -e "The VEBA database environment variable is set in your VEBA conda environments: \n\tVEBA_DATABASE=${REALPATH_DATABASE_DIRECTORY}"
echo -e "For walkthroughs on different workflows, please refer to the documentation: \n\thttps://github.com/jolespin/veba/blob/main/walkthroughs/README.md"
