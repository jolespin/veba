CODE=0714
N_JOBS=8
OUT_DIR=veba_output/assembly

mkdir -p logs
mkdir -p ${OUT_DIR}

for ID in $(cat identifiers.list);
	do N="assembly__${ID}";
	rm -f logs/${N}.*
	R1=veba_output/preprocess/${ID}/output/cleaned_1.fastq.gz
	R2=veba_output/preprocess/${ID}/output/cleaned_2.fastq.gz
	CMD="source activate VEBA-assembly_env && python ~/Algorithms/Pipelines/veba_pipeline/src/assembly.py -1 ${R1} -2 ${R2} -n ${ID} -o ${OUT_DIR} -p ${N_JOBS} -m 1024"
	qsub  -l himem7 -N ${N} -cwd -P ${CODE} -o logs/${N}.o -e logs/${N}.e -pe threaded ${N_JOBS} "${CMD}"
	done
