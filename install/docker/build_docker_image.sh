#!/usr/bin/env bash
# __version__ = "2023.7.11"

ENV_NAME=$1;
VERSION=$(head -n 1 ../../VERSION)

echo -e "================"
echo -e "VEBA v${VERSION}"
echo -e "================"

# Check if environment file exists
FILE=../environments/${ENV_NAME}.yml
if [ -f "$FILE" ]; then
    NAME=$(echo $ENV_NAME | cut -f1 -d "_" | cut -c6-)
    TAG="jolespin/veba_${NAME}:${VERSION}"
    echo -e "Creating Docker image ${TAG} for ${ENV_NAME}"
    docker build --build-arg ENV_NAME=${ENV_NAME} -t ${TAG} -f Dockerfile ../../
else 
    echo -e "$FILE does not exist."
fi
