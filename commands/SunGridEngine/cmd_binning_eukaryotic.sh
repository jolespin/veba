CODE=0714
N_JOBS=4
#DB=/usr/local/scratch/CORE/jespinoz/db/veba/eukaryotic_classification/eukaryotic
DB=/usr/local/scratch/CORE/jespinoz/db/veba/v1.0/Classify/Eukaryotic/eukaryotic

for ID in $(cat identifiers.list);
#for ID in $(cat identifiers.rerun.list);
#for ID in SRR17458614 SRR17458615 SRR17458630 SRR17458638;
#for ID in SRR17458613;
	do N="binning-eukaryotic-metabat2__${ID}";
	rm -f logs/${N}.*
	#rm -rf veba_output/binning/eukaryotic/${ID}
	#FASTA=veba_output/assembly/${ID}/output/scaffolds.fasta
	FASTA=veba_output/binning/prokaryotic/${ID}/output/unbinned.fasta
	BAM=veba_output/assembly/${ID}/output/mapped.sorted.bam
	CMD="source activate VEBA-binning-eukaryotic_env && python ~/Algorithms/Pipelines/veba_pipeline/src/binning-eukaryotic.py -f ${FASTA} -b ${BAM} -n ${ID} -p ${N_JOBS} -m 1500 -a metabat2 --metaeuk_database $DB -o veba_output/binning/eukaryotic"
	qsub -V -l himem7 -N ${N} -cwd -P ${CODE} -o logs/${N}.o -e logs/${N}.e -pe threaded ${N_JOBS} "${CMD}"
	done
