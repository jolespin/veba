
N_JOBS=4
DB_HUMAN=/usr/local/scratch/CORE/jespinoz/db/genomes/human/GRCh38.p13/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index
DB_RIBO=/usr/local/scratch/CORE/jespinoz/db/bbtools/ribokmers.fa.gz
for ID in  $(cat identifiers.list);
	do echo $ID;
	R1=Fastq/${ID}_1.fastq.gz
	R2=Fastq/${ID}_2.fastq.gz
	N=preprocessing__${ID}
	rm -f logs/${N}.*
	qsub -l himem7 -o logs/${N}.o -e logs/${N}.e -cwd -P 0714 -N ${N} -j y -pe threaded ${N_JOBS} "source activate VEBA-preprocess_env && python ~/Algorithms/Pipelines/veba_pipeline/src/preprocess.py -n ${ID} -1 ${R1} -2 ${R2} -p ${N_JOBS} -x ${DB_HUMAN} -k ${DB_RIBO}"
	done

