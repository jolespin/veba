#!/bin/bash
# __version__ = "2024.8.30"
# VEBA_DATABASE_VERSION = "VDB_v7"
# MICROEUKAYROTIC_DATABASE_VERSION = "MicroEuk_v3"
# usage: bash veba/download_databases-annotate.sh /path/to/veba_database_destination/

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

# KOfam
echo ". .. ... ..... ........ ............."
echo " * Processing KEGG profile HMM marker sets"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Annotate/KOfam/
wget -v -O - ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz | gzip -d > ${DATABASE_DIRECTORY}/Annotate/KOfam/ko_list
wget -v -c ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz -O - |  tar -xz
mv profiles ${DATABASE_DIRECTORY}/Annotate/KOfam/

## KEGG MicrobeAnnotator
#wget -v -O ${DATABASE_DIRECTORY}/MicrobeAnnotator-KEGG.tar.gz https://zenodo.org/records/10020074/files/MicrobeAnnotator-KEGG.tar.gz?download=1
#tar xvzf ${DATABASE_DIRECTORY}/MicrobeAnnotator-KEGG.tar.gz -C ${DATABASE_DIRECTORY}/Annotate --no-xattrs
#rm -rf ${DATABASE_DIRECTORY}/Annotate/._MicrobeAnnotator-KEGG
#rm -rf ${DATABASE_DIRECTORY}/MicrobeAnnotator-KEGG.tar.gz

# KEGG Pathway Profiler
mkdir -p ${DATABASE_DIRECTORY}/Annotate/KEGG-Pathway-Profiler/
build-pathway-database.py --force --download --intermediate_directory ${DATABASE_DIRECTORY}/Annotate/KEGG-Pathway-Profiler/ --database ${DATABASE_DIRECTORY}/Annotate/KEGG-Pathway-Profiler/database.pkl.gz

# Pfam
echo ". .. ... ..... ........ ............."
echo " * Processing Pfam profile HMM marker sets"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Annotate/Pfam
wget -v -P ${DATABASE_DIRECTORY}/Annotate/Pfam http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget -v -P ${DATABASE_DIRECTORY}/Annotate/Pfam http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt
wget -v -P ${DATABASE_DIRECTORY}/Annotate/Pfam http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz

# NCBIfam-AMRFinder
echo ". .. ... ..... ........ ............."
echo " * Processing NCBIfam-AMRFinder profile HMM marker sets"
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
# echo ". .. ... ..... ........ ............."
# echo "x * Processing NCBI non-redundant diamond database"
# echo ". .. ... ..... ........ ............."
# mkdir -v -p ${DATABASE_DIRECTORY}/Annotate/nr
# wget -v -P ${DATABASE_DIRECTORY} https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
# diamond makedb --in ${DATABASE_DIRECTORY}/nr.gz --db ${DATABASE_DIRECTORY}/Annotate/nr/nr.dmnd --taxonmap ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/prot.accession2taxid.FULL.gz --taxonnodes ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/nodes.dmp --taxonnames ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/names.dmp
# rm -rf ${DATABASE_DIRECTORY}/nr.gz

# UniRef
echo ". .. ... ..... ........ ............."
echo " * Processing UniRef diamond database"
echo ". .. ... ..... ........ ............."
mkdir -v -p ${DATABASE_DIRECTORY}/Annotate/UniRef

wget -v -P ${DATABASE_DIRECTORY}/Annotate/UniRef/ https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.release_note
wget -v -P ${DATABASE_DIRECTORY} https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
diamond makedb --in ${DATABASE_DIRECTORY}/uniref90.fasta.gz --db ${DATABASE_DIRECTORY}/Annotate/UniRef/uniref90.dmnd --threads ${N_JOBS}
rm -rf ${DATABASE_DIRECTORY}/uniref90.fasta.gz

wget -v -P ${DATABASE_DIRECTORY}/Annotate/UniRef/ https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref50/uniref50.release_note
wget -v -P ${DATABASE_DIRECTORY} https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
diamond makedb --in ${DATABASE_DIRECTORY}/uniref50.fasta.gz --db ${DATABASE_DIRECTORY}/Annotate/UniRef/uniref50.dmnd --threads ${N_JOBS}
rm -rf ${DATABASE_DIRECTORY}/uniref50.fasta.gz

#MiBIG
mkdir -v -p ${DATABASE_DIRECTORY}/Annotate/MIBiG
wget -v -P ${DATABASE_DIRECTORY} https://dl.secondarymetabolites.org/mibig/mibig_prot_seqs_3.1.fasta
seqkit rmdup -s ${DATABASE_DIRECTORY}/mibig_prot_seqs_3.1.fasta > ${DATABASE_DIRECTORY}/mibig_prot_seqs_3.1.rmdup.fasta
diamond makedb --in ${DATABASE_DIRECTORY}/mibig_prot_seqs_3.1.rmdup.fasta --db ${DATABASE_DIRECTORY}/Annotate/MIBiG/mibig_v3.1.dmnd --threads ${N_JOBS}
rm -rf ${DATABASE_DIRECTORY}/mibig_prot_seqs_3.1.fasta
rm -rf ${DATABASE_DIRECTORY}/mibig_prot_seqs_3.1.rmdup.fasta

# #BiG-SLiCE
# mkdir -v -p ${DATABASE_DIRECTORY}/Annotate/BiG-SLiCE
# wget -v -P ${DATABASE_DIRECTORY} https://github.com/medema-group/bigslice/releases/download/v1.0.0/bigslice-models.2020-04-27.tar.gz
# tar xzfv ${DATABASE_DIRECTORY}/bigslice-models.2020-04-27.tar.gz -C ${DATABASE_DIRECTORY}/Annotate/BiG-SLiCE
# rm -rf ${DATABASE_DIRECTORY}/bigslice-models.2020-04-27.tar.gz

# VFDB
mkdir -v -p ${DATABASE_DIRECTORY}/Annotate/VFDB
wget -v -P ${DATABASE_DIRECTORY} http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
wget -v -P ${DATABASE_DIRECTORY}/Annotate/VFDB/ http://www.mgc.ac.cn/VFs/Down/VFs.xls.gz
diamond makedb --in ${DATABASE_DIRECTORY}/VFDB_setA_pro.fas.gz --db ${DATABASE_DIRECTORY}/Annotate/VFDB/VFDB_setA_pro.dmnd --threads ${N_JOBS}
rm -rf ${DATABASE_DIRECTORY}/VFDB_setA_pro.fas.gz

# CAZy
mkdir -v -p ${DATABASE_DIRECTORY}/Annotate/CAZy
wget -v -P ${DATABASE_DIRECTORY} https://bcb.unl.edu/dbCAN2/download/CAZyDB.07262023.fa
diamond makedb --in ${DATABASE_DIRECTORY}/CAZyDB.07262023.fa --db ${DATABASE_DIRECTORY}/Annotate/CAZy/CAZyDB.07262023.dmnd --threads ${N_JOBS}
rm -rf ${DATABASE_DIRECTORY}/CAZyDB.07262023.fa


echo -e " _    _ _______ ______  _______\n  \  /  |______ |_____] |_____|\n   \/   |______ |_____] |     |"
echo -e "...................................................."
echo -e "     (annotate) Database Configuration Complete     "
echo -e "...................................................."
echo -e "[partial-database] If you are only installing individual databases you will need to set the following environment variable manually or use `update_environment_variables.sh` script:"
echo -e "[partial-database] \tVEBA_DATABASE=${REALPATH_DATABASE_DIRECTORY}"
echo -e "[partial-database] For walkthroughs on different workflows, please refer to the documentation:"
echo -e "[partial-database] \thttps://github.com/jolespin/veba/blob/main/walkthroughs/README.md"