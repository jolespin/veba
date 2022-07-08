ID=$1
N_JOBS=$2

FASTA=veba_output/assembly/${ID}/intermediate/1__assembly/scaffolds.fasta
BAM=veba_output/assembly/${ID}/intermediate/2__alignment/mapped.sorted.bam
source activate VEBA-binning-prokaryotic_env && python ~/Algorithms/Pipelines/veba_pipeline/src/binning-prokaryotic.py -f ${FASTA} -b ${BAM} -n ${ID} -p ${N_JOBS} -m 1500 -I 10

UNBINNED=veba_output/binning/prokaryotic/${ID}/output/unbinned.fasta
source activate VEBA-binning-eukaryotic_env && python ~/Algorithms/Pipelines/veba_pipeline/src/binning-eukaryotic.py -f ${UNBINNED} -b ${BAM} -n ${ID} -p ${N_JOBS}

UNBINNED=veba_output/binning/eukaryotic/${ID}/output/unbinned.fasta
source activate VEBA-binning-viral_env && python ~/Algorithms/Pipelines/veba_pipeline/src/binning-eukaryotic.py -f ${UNBINNED} -n ${ID} -p ${N_JOBS}
