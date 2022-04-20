N_JOBS=8
CODE=0714
#for NAME in prokaryotic eukaryotic viral; do 
for NAME in eukaryotic; do
	N=cluster-50_${NAME}
	rm -f logs/${N}.*
	PREFIX_CHAR=$(echo $NAME | cut -c1)
	qsub -l himem7 -N ${N} -o logs/${N}.o -e logs/${N}.e -cwd -P ${CODE} -pe threaded ${N_JOBS} -V -M jespinoz@jcvi.org -m e "source activate VEBA-cluster_env && python ~/Algorithms/Pipelines/veba_pipeline/src/cluster.py -i completeness_50/${NAME}_scaffolds_to_bins.tsv -m completeness_50/${NAME}_genomes.list -a completeness_50/${NAME}_proteins.list -o veba_output/cluster/completeness_50/${NAME} -p ${N_JOBS} --cluster_prefix ${PREFIX_CHAR^^}SLC"
	done


