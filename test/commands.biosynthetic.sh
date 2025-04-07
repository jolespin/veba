# Setup
module="biosynthetic"
job_name="${module}"
output_directory="../Analysis/veba_output/biosynthetic/"
veba_database="/home/ubuntu/jolespin-volume/s3/newatlantis-raw-veba-db-prod/VEBA/VEBA-DB_v9/"

# Make the directory and clear existing files
mkdir -p ${output_directory}/
mkdir -p logs/

params="-i ../Analysis/veba_output/misc/biosynthetic_table.tsv -o ../Analysis/veba_output/biosynthetic/ --n_jobs=-1 --veba_database ${veba_database}"
/usr/bin/time -v veba --module ${module} --params="${params}" 2> logs/${job_name}.e 1> logs/${job_name}.o
