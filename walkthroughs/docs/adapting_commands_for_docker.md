### Adapting commands for use with Docker
Containerization is a solution to using VEBA on any system and portability for resources such as AWS or Google Cloud.  To address this, I've been containerizing all of the modules.  Here is the guide for using these containers with Docker.

_____________________________________________________

#### Steps:

1. Install Docker Engine
2. Pull the image for the module
3. Run Docker container
4. Get the results
5. Use the results in a subsequent step

_____________________________________________________

**Developer Note:** 

`VEBA ≥ v2.0.0` Docker images are built-on Apple Silicon emulating Linux AMD64 architecture.  If you experience any issues with these images, please post an issue on GitHub.
_____________________________________________________

#### 1. Install Docker Engine

Refer to the [Docker documentation](https://docs.docker.com/engine/install/).  


#### 2. Pull Docker image for the module

Let's say you wanted to use the `preprocess` module.  Download the Docker image you want to use: 

```bash
# Version
VERSION=2.3.0

# Image
DOCKER_IMAGE="jolespin/veba_preprocess:${VERSION}"

# Pull image
docker image pull ${DOCKER_IMAGE}
``` 

#### 3. Run Docker container

There are a few intricacies in using `VEBA` via containers:

* Must explicitly call module scripts (i.e., old syntax).  For example, if you're using the `veba_preprocess` container, you will need to use `preprocess.py` to call the preprocessing module (not `veba --module preprocess`).
* Paths must be specified in relation to internal mounting points WITHIN the Docker container.  The way we do this is by mounting a local directory to a container directory using the `--volume` argument.   **If we don't link the local and container output directories then the output files will be marooned in the container.**

For example, here's how we would run the `preprocess.py` module.  First let's just view the help menu:

```bash
docker run --name VEBA-preprocess --rm -it ${DOCKER_IMAGE}  preprocess.py -h
```

If we wanted to run it interactively, start the container with `bash` (it automatically loads the appropriate `conda` environment):

```bash
docker run --name VEBA-preprocess --rm -it ${DOCKER_IMAGE}  bash
```

Though, it's the `preprocess.py` module so if you're running anything other than a toy dataset, then you probably want to run it on the grid so you can go to do something else while it runs. 


#### 4. Running VEBA modules from within the specified container

For some modules, you will need the *VEBA Database*.  To set the `LOCAL_DATABASE_DIRECTORY` to the `${VEBA_DATABASE}` environment variable or alternatively the path to the VEBA Database (refer to the [*VEBA Database Documentation*](https://github.com/jolespin/veba/blob/main/install/DATABASE.md#database-structure)). We mount this to `CONTAINER_DATABASE_DIRECTORY` in the container (i.e., `/volumes/database/`) with read-only permissions (`ro`).  We will provide this to `VEBA` with the `--veba_database` argument (See examples).

#### 4a. Mounting your working directory to an internal volume within the Docker container as a workspace

Below, we specify the `LOCAL_WORKING_DIRECTORY` as the current directory so we can mount to the `CONTAINER_WORKING_DIRECTORY` in the container (i.e., `/volumes/workspace/`) with read/write permissions (`rw`).

We mount these with the `--volume` argument so any file in the `LOCAL_WORKING_DIRECTORY` will be mirrored in the `CONTAINER_WORKING_DIRECTORY` including files that are generated in the container by `VEBA`.

```bash
# Directories
LOCAL_WORKING_DIRECTORY=$(pwd)
LOCAL_WORKING_DIRECTORY=$(realpath -m ${LOCAL_WORKING_DIRECTORY})
LOCAL_DATABASE_DIRECTORY=${VEBA_DATABASE} # /path/to/VEBA_DATABASE/
LOCAL_DATABASE_DIRECTORY=$(realpath -m ${LOCAL_DATABASE_DIRECTORY})

CONTAINER_WORKING_DIRECTORY=/volumes/workspace/
CONTAINER_DATABASE_DIRECTORY=/volumes/database/

# Parameters
ID=S1
R1=Fastq/${ID}_1.fastq.gz
R2=Fastq/${ID}_2.fastq.gz
NAME=VEBA-preprocess__${ID}
RELATIVE_OUTPUT_DIRECTORY=veba_output/preprocess/

# Command
CMD="preprocess.py -1 ${CONTAINER_WORKING_DIRECTORY}/${R1} -2 ${CONTAINER_WORKING_DIRECTORY}/${R2} -n ${ID} -o ${CONTAINER_WORKING_DIRECTORY}/${RELATIVE_OUTPUT_DIRECTORY} -x ${CONTAINER_DATABASE_DIRECTORY}/Contamination/chm13v2.0/chm13v2.0"

# Docker
# Version
VERSION=2.3.0

# Image
DOCKER_IMAGE="jolespin/veba_preprocess:${VERSION}"

# Run
docker run \
    --name ${NAME} \
    --rm \
    --volume ${LOCAL_WORKING_DIRECTORY}:${CONTAINER_WORKING_DIRECTORY}:rw \
    --volume ${LOCAL_DATABASE_DIRECTORY}:${CONTAINER_DATABASE_DIRECTORY}:ro \
    ${DOCKER_IMAGE} \
    ${CMD}
```

#### 4b. Mounting separate input and output directories to internal volumes within the Docker container

This usage is slightly more advanced but useful in systems where your input files are in a different directory than your analysis files.  One point to note about this usage is that you will need to update the local input volume directory if you are using results from VEBA as input...in which case, you should then switch over to the workspace usage in `4a` above.

For the input directory, we specify the `LOCAL_INPUT_DIRECTORY` and in this case it's the local current working directory that we mount to the `CONTAINER_INPUT_DIRECTORY` in the container (i.e., `/volumes/input/`) with read-only permissions (`ro`).

For the output directory, we specify the `LOCAL_OUTPUT_PARENT_DIRECTORY` and mount this to the `CONTAINER_OUTPUT_DIRECTORY` in the container (i.e., `/volumes/output/`) with read-write permissions.  After that we need to specify the `RELATIVE_OUTPUT_DIRECTORY` which is `veba_output/preprocess/` in this case. 

We mount these local directories with the `--volume` argument so any file in the `LOCAL_WORKING_DIRECTORY` will be mirrored in the `CONTAINER_INPUT_DIRECTORY` and any files created in the `CONTAINER_OUTPUT_DIRECTORY` will be mirrored in the `${LOCAL_OUTPUT_PARENT_DIRECTORY}/${RELATIVE_OUTPUT_DIRECTORY}`. 

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
VERSION=2.3.0

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
tree ${RELATIVE_OUTPUT_DIRECTORY}/${ID}/
```

or just the output: 

```
ls -lhS ${RELATIVE_OUTPUT_DIRECTORY}/${ID}/output
```

#### 5. Next step(s) of pipeline

Below is assuming you ran the `preprocess` module of VEBA and you want to follow it up with the assembly and binning. At this point, we should use the `workspace` mounting method.  We are also going to use the maximum amount of threads available on the system by settings `-p/--n_jobs=-1` similar to `scikit-learn` API.

##### 5a. Assembly

We are going to use the preprocessed reads from the previous step to perform metagenomic assembly:

```bash
# Directories
LOCAL_WORKING_DIRECTORY=$(pwd)
LOCAL_WORKING_DIRECTORY=$(realpath -m ${LOCAL_WORKING_DIRECTORY})

CONTAINER_WORKING_DIRECTORY=/volumes/workspace/

# Parameters
ID=S1
NAME=VEBA-assembly__${ID}
RELATIVE_OUTPUT_DIRECTORY=veba_output/assembly/

# If you ended up providing a contamination database: 
R1=veba_output/preprocess/${ID}/output/cleaned_1.fastq.gz
R2=veba_output/preprocess/${ID}/output/cleaned_2.fastq.gz

## If not, then you have trimmed fastq instead:
# R1=veba_output/preprocess/${ID}/output/trimmed_1.fastq.gz
# R2=veba_output/preprocess/${ID}/output/trimmed_2.fastq.gz

# Command
CMD="assembly.py -1 ${CONTAINER_WORKING_DIRECTORY}/${R1} -2 ${CONTAINER_WORKING_DIRECTORY}/${R2} -n ${ID} -o ${CONTAINER_WORKING_DIRECTORY}/${RELATIVE_OUTPUT_DIRECTORY} -p=-1"

# Docker
# Version
VERSION=2.3.0

# Image
DOCKER_IMAGE="jolespin/veba_assembly:${VERSION}"

# Run
docker run \
    --name ${NAME} \
    --rm \
    --volume ${LOCAL_OUTPUT_PARENT_DIRECTORY}:${CONTAINER_WORKING_DIRECTORY}:rw \
    ${DOCKER_IMAGE} \
    ${CMD}
```

##### 5b. Prokaryotic binning

We are going to metagenomic assembly and mapped sorted BAM file to bin out prokaryotic genomes.  Notice that we have to specify the database here:

```bash
# Directories
LOCAL_WORKING_DIRECTORY=$(pwd)
LOCAL_WORKING_DIRECTORY=$(realpath -m ${LOCAL_WORKING_DIRECTORY})
LOCAL_DATABASE_DIRECTORY=${VEBA_DATABASE} # /path/to/VEBA_DATABASE/
LOCAL_DATABASE_DIRECTORY=$(realpath -m ${LOCAL_DATABASE_DIRECTORY})

CONTAINER_WORKING_DIRECTORY=/volumes/workspace/
CONTAINER_DATABASE_DIRECTORY=/volumes/database/

# Parameters
ID=S1
NAME=VEBA-binning-prokaryotic__${ID}
FASTA=veba_output/assembly/${ID}/output/scaffolds.fasta
BAM=veba_output/assembly/${ID}/output/mapped.sorted.bam
RELATIVE_OUTPUT_DIRECTORY=veba_output/binning/prokaryotic/

# Command
CMD="binning-prokaryotic.py -f ${CONTAINER_WORKING_DIRECTORY}/${FASTA} -b ${CONTAINER_WORKING_DIRECTORY}/${BAM} -n ${ID} -o ${CONTAINER_WORKING_DIRECTORY}/${RELATIVE_OUTPUT_DIRECTORY} -p=-1 --veba_database ${CONTAINER_DATABASE_DIRECTORY}"

# Docker
# Version
VERSION=2.3.0

# Image
DOCKER_IMAGE="jolespin/veba_binning-prokaryotic:${VERSION}"

# Run
docker run \
    --name ${NAME} \
    --rm \
    --volume ${LOCAL_OUTPUT_PARENT_DIRECTORY}:${CONTAINER_WORKING_DIRECTORY}:rw \
    --volume ${LOCAL_DATABASE_DIRECTORY}:${CONTAINER_DATABASE_DIRECTORY}:ro \
    ${DOCKER_IMAGE} \
    ${CMD}
```

You get idea...

Some trickier situations might be if you need to run an utility script but most of those have very basic requirements like `pandas`, `tqdm`, `biopython`, etc.

If you have questions, feel free to post a GitHub issue with `[Question]` as the prefix.