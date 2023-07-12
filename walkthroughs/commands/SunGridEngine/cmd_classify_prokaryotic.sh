N_JOBS=4
CODE=0714
N=classify-prokaryotic
rm -f logs/${N}.*
#GENOMES=Output/Prokaryotic
INPUT=veba_output/binning/prokaryotic
#CLUSTERS=veba_output/cluster/output/clusters.tsv
qsub -l centos7 -N ${N} -o logs/${N}.o -e logs/${N}.e -cwd -P ${CODE} -pe threaded ${N_JOBS} -V -M jespinoz@jcvi.org -m e "source activate VEBA-classify_env && python ~/Algorithms/Pipelines/veba_pipeline/src/classify-prokaryotic.py -i $INPUT -p ${N_JOBS} -o veba_output/classify/prokaryotic"



