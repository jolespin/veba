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

Let's say you wanted to use the `assembly.py` module.  Download the Docker image as so: 

```
VERSION=1.1.2
docker image pull jolespin/veba_assembly:${VERSION}
``` 

#### 3. Run Docker container

One key difference with running a Docker container is that you need to specify the path IN the Docker container but it's pretty simple.  Basically, we link a local directory to a container directory using the `--volume` argument.   

For example, here's how we would run the `assembly.py` module.  First let's just look at the options:

```bash
# Version
VERSION=1.1.2

# Image
DOCKER_IMAGE="jolespin/veba_assembly:${VERSION}"

docker run --name VEBA-assembly --rm -it ${DOCKER_IMAGE}  -c "assembly.py -h"
```

If we wanted to run it interactively, start the container with `bash` (it automatically loads the appropriate `conda` environment):

```
docker run --name VEBA-assembly --rm -it ${DOCKER_IMAGE}  -c "bash"
```

Though, it's the `assembly.py` module so it you're running anything other than a toy dataset, then you probably want to run it on the grid so you can go to do something else. 

Below, we specify the `LOCAL_WORKING_DIRECTORY` which is just the current local directory.  We also need to specify the `CONTAINER_WORKING_DIRECTORY` which will be `/data/` on the volume.  We link these with the `--volume` argument so anything created in the `CONTAINER_WORKING_DIRECTORY` will get mirrored into the `LOCAL_WORKING_DIRECTORY`.  That is where we want to put the output files.

Note: If we don't specify the `${CONTAINER_WORKING_DIRECTORY}` prefix for the output then the output will be stranded in the container.

```bash

# Directories
LOCAL_WORKING_DIRECTORY=.
CONTAINER_WORKING_DIRECTORY=/data/

# Inputs
ID=S1
R1=Fastq/${ID}_1.fastq.gz
R2=Fastq/${ID}_2.fastq.gz

# Output
OUTPUT_DIRECTORY=veba_output/assembly

# Command
CMD="assembly.py -1 ${CONTAINER_WORKING_DIRECTORY}/${R1} -2 ${CONTAINER_WORKING_DIRECTORY}/${R2} -n ${ID} -o ${CONTAINER_WORKING_DIRECTORY}/${OUTPUT_DIRECTORY}"

docker run \
	--name VEBA-assembly__${ID} \
	--rm \
	--volume ${LOCAL_WORKING_DIRECTORY}:${CONTAINER_WORKING_DIRECTORY} \
	${DOCKER_IMAGE} \
	-c "${CMD}"
```

#### 4. Get the results

Now that the container has finished running the commands. Let's view the results.  If you don't have `tree` in your environment, you should download it because it's useful. `mamba install -c conda-forge tree`

To view it all:

```
tree veba_output/assembly/${ID}/
```

or just the output: 

```
ls -lhS veba_output/assembly/${ID}/output
```


#### Next steps:

Whatever you want to do.