#!/usr/bin/env bash

# Setup
module="assembly"
output_directory="../Analysis/veba_output/${module}/"
commands_file="commands.${module}.list"

# Make the directory and clear existing files
mkdir -p ${output_directory}/
mkdir -p logs/
rm -f ${commands_file}
# Iterate through all samples in the identifier list
for id in $(cat ../identifiers.list);
do
        job_name="${module}__${id}"
        echo ${job_name}
        r1="../Analysis/veba_output/preprocess/${id}/output/trimmed_1.fastq.gz"
        r2="../Analysis/veba_output/preprocess/${id}/output/trimmed_2.fastq.gz"
        params="-1 ${r1} -2 ${r2} -n ${id} -o ${output_directory} -p=16"
        cmd="/usr/bin/time -v veba --module ${module} --params \"${params}\""
        echo "${cmd} 2> logs/${job_name}.e 1> logs/${job_name}.o" >> ${commands_file}

done
