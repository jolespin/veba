# Setup
module="annotate"
job_name="${module}"
output_directory="../Analysis/veba_output/annotations/"
veba_database="/home/ubuntu/jolespin-volume/s3/newatlantis-raw-veba-db-prod/VEBA/VEBA-DB_v9/"

# Make the directory and clear existing files
mkdir -p ${output_directory}/
mkdir -p logs/
proteins="../Analysis/veba_output/misc/all_genomes.all_proteins.lt100k.faa"
identifier_mapping="../Analysis/veba_output/cluster/output/global/identifier_mapping.proteins.tsv.gz"
cat ../Analysis/veba_output/binning/*/*/output/genomes/*.faa | seqkit seq -M 99999 > ${proteins}

params="-i ${identifier_mapping} -a ${proteins} -o ${output_directory} --n_jobs=-1 --veba_database ${veba_database}"
/usr/bin/time -v veba --module ${module} --params="${params}" 2> logs/${job_name}.e 1> logs/${job_name}.o
