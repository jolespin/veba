CODE=0714
N_JOBS=4

for ID in $(cat identifiers.list);
	do N="binning-viral__${ID}";
	rm -f logs/${N}.*
	FASTA=veba_output/binning/eukaryotic/${ID}/output/unbinned.fasta
	CMD="source activate VEBA-binning_env && python ~/Algorithms/Pipelines/veba_pipeline/src/binning-viral.py -f ${FASTA} -n ${ID} -p ${N_JOBS} -m 1500 -o veba_output/binning/viral"
	qsub -V -l centos7 -N ${N} -cwd -P ${CODE} -o logs/${N}.o -e logs/${N}.e -pe threaded ${N_JOBS} "${CMD}"
	done
