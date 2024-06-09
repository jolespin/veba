#!/usr/bin/env bash
# __version__ = "2023.7.11"

for ENV_NAME in $(ls ../environments/ | grep ".yml"); 
do 
    ENV_NAME=$(echo $ENV_NAME | cut -f1 -d ".")
    echo $ENV_NAME
    time(bash build_docker_image.sh ${ENV_NAME})
done
