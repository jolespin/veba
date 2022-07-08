CODE=0714
N_JOBS=1

#for ID in $(cat identifiers.list);
#	do N="annotate__${ID}";
#	rm -f logs/${N}.*
#	OUT_DIR="veba_output/annotation/${ID}"
#	mkdir -p ${OUT_DIR}
#	#PROTEIN=veba_output/binning/prokaryotic/${ID}/intermediate/2__prodigal/gene_models.faa
#	PROTEIN=veba_output/binning/prokaryotic/${ID}/output/genomes/
#	MAPPING=veba_output/binning/prokaryotic/${ID}/output/genomes/identifier_mapping.tsv
#	CMD="source activate VEBA-annotate_env && python ~/Algorithms/Pipelines/veba_pipeline/src/annotate.py -a ${PROTEIN} -o ${OUT_DIR} -p ${N_JOBS} -i ${MAPPING}"
#	qsub  -l centos7 -N ${N} -cwd -P ${CODE} -o logs/${N}.o -e logs/${N}.e -pe threaded ${N_JOBS} "${CMD}"
#	done

N="annotate"
rm -f logs/${N}.*
OUT_DIR="veba_output/annotation/"
mkdir -p ${OUT_DIR}
PROTEIN=all_proteins.faa
MAPPING=identifier_mapping.tsv
CMD="source activate VEBA-annotate_env && python ~/Algorithms/Pipelines/veba_pipeline/src/annotate.py -a ${PROTEIN} -o ${OUT_DIR} -p ${N_JOBS} -i ${MAPPING}"
qsub  -l centos7 -N ${N} -cwd -P ${CODE} -o logs/${N}.o -e logs/${N}.e -pe threaded ${N_JOBS} "${CMD}"
