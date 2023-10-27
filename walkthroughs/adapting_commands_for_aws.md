### Adapting commands for use with Docker and AWS
Containerization is a solution to using VEBA on any system and portability for resources such as AWS or Google Cloud.  Here is the guide for using these containers specifically with AWS.

_____________________________________________________

#### Steps:

1. Set up AWS infrastructure
2. Create and register a job definition
3. Submit job definition

_____________________________________________________


#### 1. Set up AWS infrastructure

Out of scope for this tutorial but essentially you need to do the following: 

* Set up AWS EFS (Elastic File System) via Terraform to read/write/mount data
* Compile database in EFS
* Create compute environment
* Create job queue linked to compute environment


#### 2. Create and register a job definition

Once the job queue is properly set up, next is to create a job definition and then submit the job definition to the queue.

The preferred way to submit jobs with AWS Batch is using json files for the job definition through Fargate.  

Here is a template you can use for a job definition.  

This job definition pulls the [jolespin/veba_preprocess](https://hub.docker.com/r/jolespin/veba_preprocess/tags) Docker image and mounts EFS directories to volumes within the Docker container.  The actual job runs the [preprocess.py module](https://github.com/jolespin/veba/tree/main/src#preprocesspy) of VEBA for a sample called [S1](https://zenodo.org/record/7946802). 


```json
{
  "jobDefinitionName": "preprocess__S1",
  "type": "container",
  "containerProperties": {
    "image": "jolespin/veba_preprocess:1.3.0",
    "command": [
      "preprocess.py",
      "-1",
      "/volumes/input/Fastq/S1_1.fastq.gz",
      "-2",
      "/volumes/input/Fastq/S1_2.fastq.gz",
      "-n",
      "1",
      "-o",
      "/volumes/output/veba_output/preprocess",
      "-p",
      "16"
      "-x",
      "/volumes/database/Contamination/chm13v2.0/chm13v2.0"
    ],
    "jobRoleArn": "arn:aws:iam::xxx:role/ecsTaskExecutionRole",
    "executionRoleArn": "arn:aws:iam::xxx:role/ecsTaskExecutionRole",
    "volumes": [
      {
        "name": "efs-volume-database",
        "efsVolumeConfiguration": {
          "fileSystemId": "fs-xxx",
          "transitEncryption": "ENABLED",
          "rootDirectory": "databases/veba/VDB_v5.1/"
        }
      },
      {
        "name": "efs-volume-input",
        "efsVolumeConfiguration": {
          "fileSystemId": "fs-xxx",
          "transitEncryption": "ENABLED",
          "rootDirectory": "path/to/efs/input/"
        }
      },
      {
        "name": "efs-volume-output",
        "efsVolumeConfiguration": {
          "fileSystemId": "fs-xxx",
          "transitEncryption": "ENABLED",
          "rootDirectory": "path/to/efs/output/"
        }
      }
    ],
    "mountPoints": [
    {
        "sourceVolume": "efs-volume-database",
        "containerPath": "/volumes/database",
        "readOnly": true
      },
      {
        "sourceVolume": "efs-volume-input",
        "containerPath": "/volumes/input",
        "readOnly": true
      },
      {
        "sourceVolume": "efs-volume-output",
        "containerPath": "/volumes/output",
        "readOnly": false
      }
    ],
    "environment": [],
    "ulimits": [],
    "resourceRequirements": [
      {
        "value": "16.0",
        "type": "VCPU"
      },
      {
        "value": "8000",
        "type": "MEMORY"
      }
    ],
    "networkConfiguration": {
      "assignPublicIp": "ENABLED"
    },
    "fargatePlatformConfiguration": {
      "platformVersion": "LATEST"
    },
    "ephemeralStorage": {
      "sizeInGiB": 40
    }
  },
  "tags": {
    "Name": "preprocess__S1"
  },
  "platformCapabilities": [
    "FARGATE"
  ]
}
```

Now register the job definition: 

```
FILE=/path/to/preprocess/S1.json
aws batch register-job-definition --cli-input-json file://${FILE}
```

#### 3. Run Docker container

Next step is to submit the job to the queue.

```
aws batch submit-job --job-definition ${JOB} --job-name ${JOB} --job-queue ${QUEUE}
```