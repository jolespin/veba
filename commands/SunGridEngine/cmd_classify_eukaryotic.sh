CODE=0714
N=classify-eukaryotic
rm -rf logs/${N}.*
qsub -M jespinoz@jcvi.org -m e -l centos7 -N ${N} -o logs/${N}.o -e logs/${N}.e -cwd -P ${CODE} "source activate VEBA-classify_env && python ~/Algorithms/Pipelines/veba_pipeline/src/classify-eukaryotic.py -i veba_output/binning/eukaryotic/ -o veba_output/classify/eukaryotic -l 1.0"


