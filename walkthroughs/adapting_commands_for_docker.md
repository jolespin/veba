### Adapting commands for use with Docker
Containerization is a solution to using VEBA on any system and portability for resources such as AWS or Google Cloud.  To address this, I've been containerizing all of the modules.  Here is the guide for using these containers.

_____________________________________________________

#### Steps:

1. Install Docker Engine
2. Pull the image for the module
3. Run Docker container
4. Get the results
5. Use Singularity (Coming Soon)

_____________________________________________________


#### 1. Install Docker Engine

Refer to the [Docker documentation](https://docs.docker.com/engine/install/).  


#### 2. Pull Docker image for the module

Let's say you wanted to use the `preprocess` module.  Download the Docker image as so: 

```
VERSION=1.4.0
docker image pull jolespin/veba_preprocess:${VERSION}
``` 

#### 3. Run Docker container

One key difference with running a Docker container is that you need to specify the path IN the Docker container but it's pretty simple.  Basically, we link a local directory to a container directory using the `--volume` argument.   

For example, here's how we would run the `preprocess.py` module.  First let's just look at the options:

```bash
# Version
VERSION=1.4.0

# Image
DOCKER_IMAGE="jolespin/veba_preprocess:${VERSION}"

docker run --name VEBA-preprocess --rm -it ${DOCKER_IMAGE}  preprocess.py -h
```

If we wanted to run it interactively, start the container with `bash` (it automatically loads the appropriate `conda` environment):

```
docker run --name VEBA-preprocess --rm -it ${DOCKER_IMAGE}  bash
```

Though, it's the `preprocess.py` module so it you're running anything other than a toy dataset, then you probably want to run it on the grid so you can go to do something else. 

Below, we specify the `LOCAL_WORKING_DIRECTORY` and in this case it's the local current working directory.  We need to link the `CONTAINER_INPUT_DIRECTORY` to the `/volumes/input/` on the container.  

For the output directory, it requires an additional step.  We first need to specify the `LOCAL_OUTPUT_PARENT_DIRECTORY` and link this to the `CONTAINER_OUTPUT_DIRECTORY` which will be `/volumes/output/` on the container.  After that we need to specify the `RELATIVE_OUTPUT_DIRECTORY` which is `veba_output/preprocess/` in this case. 

We link these with the `--volume` argument so any file in the `LOCAL_WORKING_DIRECTORY` will be mirrored in the `CONTAINER_INPUT_DIRECTORY` and any files created in the `CONTAINER_OUTPUT_DIRECTORY` (i.e., the `RELATIVE_OUTPUT_DIRECTORY`) will be mirrored in the `LOCAL_OUTPUT_PARENT_DIRECTORY`. 

For some modules, you will need the VEBA Database.  To set the `LOCAL_DATABASE_DIRECTORY` to the `${VEBA_DATABASE}` environment variable or alternatively the path to the VEBA Database (refer to the [*VEBA Database Documentation*](https://github.com/jolespin/veba/blob/main/install/DATABASE.md#database-structure)). We mount this to `CONTAINER_DATABASE_DIRECTORY` which is `/volumes/database/` volume in the container. 


**Note:**

If we don't link the local and container output directories then the output files will be stranded in the container.

```bash
# Directories
LOCAL_WORKING_DIRECTORY=$(pwd)
LOCAL_WORKING_DIRECTORY=$(realpath -m ${LOCAL_WORKING_DIRECTORY})
LOCAL_OUTPUT_PARENT_DIRECTORY=../
LOCAL_OUTPUT_PARENT_DIRECTORY=$(realpath -m ${LOCAL_OUTPUT_PARENT_DIRECTORY})
LOCAL_DATABASE_DIRECTORY=${VEBA_DATABASE} # /path/to/VEBA_DATABASE/
LOCAL_DATABASE_DIRECTORY=$(realpath -m ${LOCAL_DATABASE_DIRECTORY})

CONTAINER_INPUT_DIRECTORY=/volumes/input/
CONTAINER_OUTPUT_DIRECTORY=/volumes/output/
CONTAINER_DATABASE_DIRECTORY=/volumes/database/

# Parameters
ID=S1
R1=Fastq/${ID}_1.fastq.gz
R2=Fastq/${ID}_2.fastq.gz
NAME=VEBA-preprocess__${ID}
RELATIVE_OUTPUT_DIRECTORY=veba_output/preprocess/

# Command
CMD="preprocess.py -1 ${CONTAINER_INPUT_DIRECTORY}/${R1} -2 ${CONTAINER_INPUT_DIRECTORY}/${R2} -n ${ID} -o ${CONTAINER_OUTPUT_DIRECTORY}/${RELATIVE_OUTPUT_DIRECTORY} -x ${CONTAINER_DATABASE_DIRECTORY}/Contamination/chm13v2.0/chm13v2.0"

# Docker
# Version
VERSION=1.4.0

# Image
DOCKER_IMAGE="jolespin/veba_preprocess:${VERSION}"

# Run
docker run \
    --name ${NAME} \
    --rm \
    --volume ${LOCAL_WORKING_DIRECTORY}:${CONTAINER_INPUT_DIRECTORY}:ro \
    --volume ${LOCAL_OUTPUT_PARENT_DIRECTORY}:${CONTAINER_OUTPUT_DIRECTORY}:rw \
    --volume ${LOCAL_DATABASE_DIRECTORY}:${CONTAINER_DATABASE_DIRECTORY}:ro \
    ${DOCKER_IMAGE} \
    ${CMD}

```


#### 4. Get the results

Now that the container has finished running the commands. Let's view the results.  If you don't have `tree` in your environment, you should download it because it's useful. `mamba install -c conda-forge tree`

To view it all:

```
tree ${LOCAL_OUTPUT_PARENT_DIRECTORY}/veba_output/preprocess/${ID}/
```

or just the output: 

```
ls -lhS ${LOCAL_OUTPUT_PARENT_DIRECTORY}/veba_output/preprocess/${ID}/output
```
