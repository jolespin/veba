# Create database
DATABASE_DIRECTORY=${1:-"."}

# Database structure
echo "Creating directories"
mkdir -vp $DATABASE_DIRECTORY
mkdir -vp ${DATABASE_DIRECTORY}/Annotate
mkdir -vp ${DATABASE_DIRECTORY}/Classify
mkdir -vp ${DATABASE_DIRECTORY}/MarkerSets

# Versions
DATE=$(date)
echo $DATE > ${DATABASE_DIRECTORY}/ACCESS_DATE

# NCBI Taxonomy
echo "Processing NCBITaxonomy"
mkdir -vp ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy
wget -vP ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz 
wget -vP ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
python -c 'import sys; from ete3 import NCBITaxa; NCBITaxa(taxdump_file="%s/Classify/NCBITaxonomy/taxdump.tar.gz"%(sys.argv[1]), dbfile="%s/Classify/NCBITaxonomy/taxa.sqlite"%(sys.argv[1]))' $DATABASE_DIRECTORY
tar -xzfv ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/taxdump.tar.gz

# GTDB-Tk
echo "Processing GTDB-Tk"
wget -vP ${DATABASE_DIRECTORY}/Classify https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar xvzf ${DATABASE_DIRECTORY}/Classify/gtdbtk_data.tar.gz
mv ${DATABASE_DIRECTORY}/Classify/gtdbtk_data ${DATABASE_DIRECTORY}/Classify/gtdbtk
rm -rf ${DATABASE_DIRECTORY}/Classify/gtdbtk_data.tar.gz

# CheckV
echo "Processing CheckV"
wget -vP ${DATABASE_DIRECTORY}/Classify https://portal.nersc.gov/CheckV/checkv-db-v1.0.tar.gz
tar xvzf ${DATABASE_DIRECTORY}/Classify/checkv-db-v1.0.tar.gz
mv ${DATABASE_DIRECTORY}/Classify/checkv-db-v1.0 ${DATABASE_DIRECTORY}/Classify/checkv
rm -rf ${DATABASE_DIRECTORY}/Classify/checkv-db-v1.0.tar.gz

# CheckM
echo "Processing CheckM"
wget -vP ${DATABASE_DIRECTORY}/Classify https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar xvzf ${DATABASE_DIRECTORY}/Classify/checkm_data_2015_01_16.tar.gz
mv ${DATABASE_DIRECTORY}/Classify/checkm_data_2015_01_16 ${DATABASE_DIRECTORY}/Classify/checkm
rm -rf ${DATABASE_DIRECTORY}/Classify/checkm_data_2015_01_16.tar.gz

# Microeukaryotic 
echo "Processing Microeukaryotic MMSEQS2 database"
# Download from FigShare
tar -xzf ${DATABASE_DIRECTORY}/Classify/Microeukaryotic.tar.gz
mmseqs createdb ${DATABASE_DIRECTORY}/Classify/Microeukaryotic/reference.rmdup.iupac.relabeled.no_deprecated.complete_lineage.faa.gz microeukaryotic
rm -rf ${DATABASE_DIRECTORY}/Classify/Microeukaryotic.tar.gz

# MarkerSets
echo "Processing profile HMM marker sets"
wget -O MarkerSets.tar.gz https://ndownloader.figshare.com/files/34844712
tar -xzf MarkerSets.tar.gz
rm -rf MarkerSets.tar.gz

# KOFAMSCAN
echo "Processing KOFAMSCAN profile HMM marker sets"
mkdir -vp ${DATABASE_DIRECTORY}/Annotate/kofam
wget -O - ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz | gzip -d > ${DATABASE_DIRECTORY}/Annotate/kofam/ko_list
wget -c ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz -O - |  tar -xz
mv profiles ${DATABASE_DIRECTORY}/Annotate/kofam/

# Pfam
echo "Processing Pfam profile HMM marker sets"
mkdir -vp ${DATABASE_DIRECTORY}/Annotate/pfam
wget -P ${DATABASE_DIRECTORY}/Annotate/nr/pfam http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz

# NCBI non-redundant
echo "Processing NCBI non-redundant diamond database"
mkdir -vp ${DATABASE_DIRECTORY}/Annotate/nr
wget -P ${DATABASE_DIRECTORY}/Annotate/nr https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
diamond makedb --in ${DATABASE_DIRECTORY}/Annotate/nr.gz --db ${DATABASE_DIRECTORY}/Annotate/nr.dmnd --taxonmap ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/prot.accession2taxid.FULL.gz --taxonnodes ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/nodes.dmp --taxonnames ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/names.dmp

echo "Adding the following environment variable to VEBA environments: export VEBA_DATABASE=${1}"
CONDA_BASE=$(which conda | python -c "import sys; print('/'.join(sys.stdin.read().split('/')[:-2]))")

# VEBA
for ENV_PREFIX in ${CONDA_BASE}/envs/VEBA-*; do 
    echo $NAME;
    mkdir -vp ${ENV_PREFIX}/etc/conda/activate.d/
    mkdir -vp ${ENV_PREFIX}/etc/conda/deactivate.d/
    echo "export VEBA_DATABASE=${1}" > ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset VEBA_DATABASE" > ${ENV_PREFIX}/etc/conda/deactivate.d/veba.sh
    done

#GTDB-Tk
for ENV_PREFIX in VEBA-binning-prokaryotic_env VEBA-classify-prokaryotic_env; do 
    GTDBTK_DATABASE_VERSION=$(ls ${DATABASE_DIRECTORY}/Classify/gtdbtk)
    echo "export GTDBTK_DATA_PATH=${DATABASE_DIRECTORY}/Classify/gtdbtk/${GTDBTK_DATABASE_VERSION}/}" >> ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset GTDBTK_DATA_PATH" >> ${ENV_PREFIX}/etc/conda/activate.d/veba.sh

# CheckV
for ENV_PREFIX in VEBA-binning-viral_env; do 
    echo "export CHECKVDB=${DATABASE_DIRECTORY}/Classify/checkv}" >> ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
    echo "unset CHECKVDB" >> ${ENV_PREFIX}/etc/conda/activate.d/veba.sh
