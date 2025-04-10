### Adapting commands for use with Singularity
Containerization is a solution to using VEBA on any system and portability for resources such as AWS or Google Cloud.  To address this, I've been containerizing all of the modules.  Here is the guide for using these containers with `Singularity`.  For more detailed support on containers, please refer to the [Docker Walkthrough](adapting_commands_for_docker.md) for more details.

**Disclaimer:** I am the *least* familiar with `Singularity` and can only provide limited support based on the default version available on our HPC.  This workflow is based on the following [conversation](https://github.com/jolespin/veba/issues/45).

_____________________________________________________

#### Steps:

1. Load the Singularity module
2. Configure the image for the module
3. Run Singularity container
4. Get the results


_____________________________________________________

#### 1. Load the Singularity module

I'm not sure how this is installed but if you're using a HPC then contact IT to see if it's available and their recommended usage. On the our HPC, this is how it is loaded but your version or activation might be different:

```bash
declare -xr SINGULARITY_MODULE='singularitypro/3.9'

module purge
module load "${SINGULARITY_MODULE}"
```

#### 2. Configure the image for the module

Let's say you wanted to use the `binning-prokaryotic` module.  Pull the Docker image you want to use and configure it for Singularity: 

```bash
# Version
VERSION=2.5.0

# Image
MODULE="veba_binning-prokaryotic"
DOCKER_IMAGE="jolespin/${MODULE}:${VERSION}"

# Configure image
mkdir -p containers/
SINGULARITY_IMAGE=containers/${MODULE}__${VERSION}.sif
singularity pull ${SINGULARITY_IMAGE} docker://jolespin/${DOCKER_IMAGE}
``` 

#### 3. Run Singuarlity container

There are a few intricacies in using `VEBA` via containers:

* Must explicitly call module scripts (i.e., old syntax).  For example, if you're using the `veba_binning-prokaryotic` container, you will need to use `binning-prokaryotic.py` to call the prokaryotic binning module (not `veba --module binning-prokaryotic`).
* Paths must be specified in relation to internal mounting points WITHIN the Docker container.  The way we do this is by mounting a local directory to a container directory using the `--volume` argument.   **If we don't link the local and container output directories then the output files will be marooned in the container.**

These are described in detail in the [Docker Walkthrough](adapting_commands_for_docker.md).

⚠️ With Singularity, you need to export your path (this is not needed with Docker):

```bash
export PATH=/usr/bin/:/opt/conda/bin/; export CONDA_PREFIX=/opt/conda/;
```

#### 4. Running VEBA modules from within the specified container

For some modules, you will need the *VEBA Database*.  To set the `LOCAL_DATABASE_DIRECTORY` to the `${VEBA_DATABASE}` environment variable or alternatively the path to the VEBA Database (refer to the [*VEBA Database Documentation*](https://github.com/jolespin/veba/blob/main/install/DATABASE.md#database-structure)). We mount this to `CONTAINER_DATABASE_DIRECTORY` in the container (i.e., `/volumes/database/`) with read-only permissions (`ro`).  We will provide this to `VEBA` with the `--veba_database` argument (See examples).

#### Mounting your working directory to an internal volume within the Singularity container as a workspace

Below, we specify the `LOCAL_WORKING_DIRECTORY` as the current directory so we can mount to the `CONTAINER_WORKING_DIRECTORY` in the container (i.e., `/volumes/workspace/`) with read/write permissions (`rw`).

We mount these with the `--volume` argument so any file in the `LOCAL_WORKING_DIRECTORY` will be mirrored in the `CONTAINER_WORKING_DIRECTORY` including files that are generated in the container by `VEBA`.

```bash
# Local directories
VEBA_DATABASE=path/to/veba_database/
LOCAL_WORKING_DIRECTORY=$(pwd)
LOCAL_WORKING_DIRECTORY=$(realpath -m ${LOCAL_WORKING_DIRECTORY})
LOCAL_DATABASE_DIRECTORY=${VEBA_DATABASE}
LOCAL_DATABASE_DIRECTORY=$(realpath -m ${LOCAL_DATABASE_DIRECTORY})

# Container directories
CONTAINER_WORKING_DIRECTORY=/volumes/workspace/
CONTAINER_DATABASE_DIRECTORY=/volumes/database/

FASTA=${CONTAINER_WORKING_DIRECTORY}/veba_output/assembly/S1/output/scaffolds.fasta
BAM=${CONTAINER_WORKING_DIRECTORY}/veba_output/assembly/S1/output/mapped.sorted.bam
OUTPUT_DIRECTORY=${CONTAINER_WORKING_DIRECTORY}/veba_output/binning/prokaryotic/
NAME="S1"

SINGULARITY_IMAGE="containers/veba_binning-prokaryotic__${VERSION}.sif"
singularity exec \
    --bind ${LOCAL_WORKING_DIRECTORY}:${CONTAINER_WORKSPACE_DIRECTORY},${LOCAL_DATABASE_DIRECTORY}:${CONTAINER_DATABASE_DIRECTORY} \
     ${SINGULARITY_IMAGE} \
    bash -c \
"export PATH=/opt/conda/bin:/usr/bin; export CONDA_PREFIX=/opt/conda; binning-prokaryotic.py -f ${FASTA} -b ${BAM} -n ${NAME} -o ${OUTPUT_DIRECTORY} --veba_database ${CONTAINER_DATABASE_DIRECTORY}"
```

You get idea...

Some trickier situations might be if you need to run an utility script but most of those have very basic requirements like `pandas`, `tqdm`, `biopython`, etc.

If you have questions, feel free to post a GitHub issue with `[Question]` as the prefix.