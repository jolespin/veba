N_JOBS=16
N=index-global
rm -f logs/${N}.*
#qsub -N ${N} -o logs/${N}.o -e logs/${N}.e -l centos7 -cwd -P 0714 -pe threaded ${N_JOBS} -M jespinoz@jcvi.org -m e "source activate VEBA-mapping_env && python ~/Algorithms/Pipelines/veba_pipeline/src/index.py -r genomes.list -g gene_models.list -o veba_output/index/global/ -p ${N_JOBS}"

N=index-local
rm -f logs/${N}.*
qsub -N ${N} -o logs/${N}.o -e logs/${N}.e -l centos7 -cwd -P 0714 -pe threaded ${N_JOBS} -M jespinoz@jcvi.org -m e "source activate VEBA-mapping_env && python ~/Algorithms/Pipelines/veba_pipeline/src/index.py -r sample_to_genome.tsv -g sample_to_gff.tsv -o veba_output/index/local/ -p ${N_JOBS}"

