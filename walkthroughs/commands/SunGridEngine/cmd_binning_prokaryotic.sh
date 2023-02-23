CODE=0714
N_JOBS=16

for ID in $(cat identifiers.list);
	do N="binning-prokaryotic__${ID}";
	rm -f logs/${N}.*
	#FASTA=veba_output/assembly/${ID}/output/scaffolds.fasta
	#BAM=veba_output/assembly/${ID}/output/mapped.sorted.bam
	FASTA=veba_output/assembly/${ID}/intermediate/1__assembly/scaffolds.fasta
	BAM=veba_output/assembly/${ID}/intermediate/2__alignment/mapped.sorted.bam
	CMD="source activate VEBA-binning-prokaryotic_env && python ~/Algorithms/Pipelines/veba_pipeline/src/binning-prokaryotic.py -f ${FASTA} -b ${BAM} -n ${ID} -p ${N_JOBS} -m 1500 -I 10"
	qsub -V -l himem7 -N ${N} -cwd -P ${CODE} -o logs/${N}.o -e logs/${N}.e -pe threaded ${N_JOBS} "${CMD}"
	done
