#!/usr/bin/env bash
# __version__ = "2023.7.12"

TAG="1.2.0"


# Assembly
NAME="VEBA-assembly"
DOCKER_IMAGE="jolespin/veba_assembly:${TAG}"
CMD="export VEBA_DATABASE=/volumes/database/ && (assembly.py -h >/dev/null && echo '${NAME} Passed') || (echo '! ${NAME}' Failed)"
docker run --name ${NAME} --rm -it ${DOCKER_IMAGE}  -c "${CMD}"

# Annotate
NAME="VEBA-annotate"
DOCKER_IMAGE="jolespin/veba_annotate:${TAG}"
CMD="export VEBA_DATABASE=/volumes/database/ && (annotate.py -h >/dev/null && echo '${NAME} Passed') || (echo '! ${NAME}' Failed)"
docker run --name ${NAME} --rm -it ${DOCKER_IMAGE}  -c "${CMD}"

# Binning Eukaryotic
NAME="VEBA-binning-eukaryotic"
DOCKER_IMAGE="jolespin/veba_binning-eukaryotic:${TAG}"
CMD="export VEBA_DATABASE=/volumes/database/ && (binning-eukaryotic.py -h >/dev/null && echo '${NAME} Passed') || (echo '! ${NAME}' Failed)"
docker run --name ${NAME} --rm -it ${DOCKER_IMAGE}  -c "${CMD}"

# Binning Prokaryotic
NAME="VEBA-binning-prokaryotic"
DOCKER_IMAGE="jolespin/veba_binning-prokaryotic:${TAG}"
CMD="export VEBA_DATABASE=/volumes/database/ && (binning-prokaryotic.py -h >/dev/null && echo '${NAME} Passed') || (echo '! ${NAME}' Failed)"
docker run --name ${NAME} --rm -it ${DOCKER_IMAGE}  -c "${CMD}"

# Binning Viral
NAME="VEBA-binning-viral"
DOCKER_IMAGE="jolespin/veba_binning-viral:${TAG}"
CMD="export VEBA_DATABASE=/volumes/database/ && (binning-prokaryotic.py -h >/dev/null && echo '${NAME} Passed') || (echo '! ${NAME}' Failed)"
docker run --name ${NAME} --rm -it ${DOCKER_IMAGE}  -c "${CMD}"

# Classify
NAME="VEBA-classify"
DOCKER_IMAGE="jolespin/veba_classify:${TAG}"
CMD="export VEBA_DATABASE=/volumes/database/ && ((classify-eukaryotic.py -h && classify-prokaryotic.py -h && classify-viral.py) >/dev/null && echo '${NAME} Passed') || (echo '! ${NAME}' Failed)"
docker run --name ${NAME} --rm -it ${DOCKER_IMAGE}  -c "${CMD}"

# Cluster
NAME="VEBA-cluster"
DOCKER_IMAGE="jolespin/veba_cluster:${TAG}"
CMD="export VEBA_DATABASE=/volumes/database/ && (cluster.py -h >/dev/null && echo '${NAME} Passed') || (echo '! ${NAME}' Failed)"
docker run --name ${NAME} --rm -it ${DOCKER_IMAGE}  -c "${CMD}"

# Mapping
NAME="VEBA-mapping"
DOCKER_IMAGE="jolespin/veba_mapping:${TAG}"
CMD="export VEBA_DATABASE=/volumes/database/ && ((index.py -h && mapping.py -h) >/dev/null && echo '${NAME} Passed') || (echo '! ${NAME}' Failed)"
docker run --name ${NAME} --rm -it ${DOCKER_IMAGE}  -c "${CMD}"

# Phylogeny
NAME="VEBA-phylogeny"
DOCKER_IMAGE="jolespin/veba_phylogeny:${TAG}"
CMD="export VEBA_DATABASE=/volumes/database/ && (phylogeny.py -h >/dev/null && echo '${NAME} Passed') || (echo '! ${NAME}' Failed)"
docker run --name ${NAME} --rm -it ${DOCKER_IMAGE}  -c "${CMD}"

# Preprocess
NAME="VEBA-preprocess"
DOCKER_IMAGE="jolespin/veba_preprocess:${TAG}"
CMD="export VEBA_DATABASE=/volumes/database/ && (preprocess.py -h >/dev/null && echo '${NAME} Passed') || (echo '! ${NAME}' Failed)"
docker run --name ${NAME} --rm -it ${DOCKER_IMAGE}  -c "${CMD}"
