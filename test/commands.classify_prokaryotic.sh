# Setup
organism_type="prokaryotic"
module="classify-${organism_type}"
job_name="${module}"
output_directory="../Analysis/veba_output/classify/${organism_type}"
veba_database="/home/ubuntu/jolespin-volume/s3/newatlantis-raw-veba-db-prod/VEBA/VEBA-DB_v9/"

# Make the directory and clear existing files
mkdir -p ${output_directory}/
mkdir -p logs/

params="-i ../Analysis/veba_output/binning/${organism_type}/ -c ../Analysis/veba_output/cluster/output/global/mags_to_slcs.tsv -o ${output_directory} -p=-1 --veba_database ${veba_database}"
/usr/bin/time -v veba --module ${module} --params="${params}" 2> logs/${job_name}.e 1> logs/${job_name}.o

