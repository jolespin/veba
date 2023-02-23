CODE=0714
N_JOBS=4

mkdir -p logs
mkdir -p ${OUT_DIR}

for ID in $(cat identifiers.list);
	do N="mapping-global__${ID}";
	rm -f logs/${N}.*
	R1=veba_output/preprocess/${ID}/output/cleaned_1.fastq.gz
	R2=veba_output/preprocess/${ID}/output/cleaned_2.fastq.gz
	INDEX=veba_output/index/global/output/
	OUT_DIR=veba_output/mapping/global/${ID}
	CMD="source activate VEBA-mapping_env && python ~/Algorithms/Pipelines/veba_pipeline/src/mapping.py -1 ${R1} -2 ${R2} -n ${ID} -o ${OUT_DIR} -p ${N_JOBS} -x ${INDEX}"
	qsub  -l centos7 -N ${N} -cwd -P ${CODE} -o logs/${N}.o -e logs/${N}.e -pe threaded ${N_JOBS} "${CMD}"
	done
