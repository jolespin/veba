## Installation and Database Configuration Guide
____________________________________________________________
### Software installation
One issue with having large-scale pipeline suites with open-source software is the issue of dependencies.  One solution for this is to have a modular software structure where each module has its own `conda` environment (or `Docker` container).  This allows for minimizing dependency constraints as this software suite uses an array of diverse packages from different developers. 

The basis for these environments is creating a separate environment for each module with the `VEBA-` prefix and `_env` as the suffix.  For example `VEBA-assembly_env` or `VEBA-binning-prokaryotic_env`.  Because of this, `VEBA` is currently not available as a `conda` package but each module will be in the near future.  In the meantime, please use the `veba/install/install.sh` script which installs each environment from the yaml files in `veba/install/environments/`. After installing the environments, use the `veba/install/download_databases.sh` script to download and configure the databases while also adding the environment variables to the activate/deactivate scripts in each environment.  

**For advanced users installing custom builds or using preexisting databases, please refer to the FAQ and read the well-documented scripts.**

Each major version will be packaged as a [release](https://github.com/jolespin/veba/releases) which will include a log of module and script versions. 

____________________________________________________________

### VEBA Database: 

Please refer to the [VEBA Database](DATABASE.md) documentation.
____________________________________________________________

### Install via Conda (Recommended):

Please download [miniconda distribution](https://docs.conda.io/projects/miniconda/en/latest/) (what I use) or the full [anaconda distribution](https://www.anaconda.com/products/distribution). 

Currently, **Conda environments for VEBA are ONLY configured for Linux** and, due to the large databases, this software was designed to be used via HPC. 

**There are 3 steps to install *VEBA*:**

* Download repository from GitHub

* Install conda environments

* Download/configure databases

**0. Clean up your conda installation [Optional, but highly recommended]**

The `VEBA` installation is going to configure some `conda` environments for you and some of them have quite a bit of packages.  To minimize the likelihood of [weird errors](https://forum.qiime2.org/t/valueerror-unsupported-format-character-t-0x54-at-index-3312-when-creating-environment-from-environment-file/25237), it's recommended to do the following:


* Use this as your [`~/.condarc`](https://conda.io/projects/conda/en/latest/user-guide/configuration/use-condarc.html).  If you're not familiar with the `.condarc` file, then you probably don't have one configured.  You can use an editor like [nano](https://anaconda.org/conda-forge/nano) (which is what I use), [vim](https://anaconda.org/conda-forge/vim), or [emacs](https://anaconda.org/conda-forge/emacs) to copy/paste the following into [`~/.condarc`](condarc).
	
	```
	channel_priority: flexible
	channels:
	  - conda-forge
	  - bioconda
	  - jolespin
	  - defaults
	  - qiime2
		
	report_errors: true
	```
	
* Make sure your `conda` is [initialized](https://docs.conda.io/projects/conda/en/latest/commands/init.html). I use `bash` for my initialization.  
	
	```
	conda init bash
	```
	
	
* Clean up your `conda` environment with the following command.
	
	```
	conda clean --all -y
	```
	
* Update your `conda`.
	
	```
	conda update -n base --all -y
	```
	
* Install and update [`mamba`](https://mamba.readthedocs.io/en/latest/installation.html).
	
	```
	conda install -c conda-forge mamba -y
	conda update mamba -y
	```

**1. Download repository**

```
# For stable version, download and decompress the tarball:

VERSION="2.2.0"
wget https://github.com/jolespin/veba/archive/refs/tags/v${VERSION}.tar.gz
tar -xvf v${VERSION}.tar.gz && mv veba-${VERSION} veba

# Alternative download
# wget https://github.com/jolespin/veba/releases/download/v${VERSION}/v${VERSION}.zip
# unzip -d veba v${VERSION}.zip

# For developmental version, clone the repository:
# Note: This is not recommended because between v2.1.0 and v2.2.0, case changes were introduced (KOFAM -> KOfam)
# and these changes are not updating on GitHub.  Please use official releases instead of pulling the repo:
# git clone --branch devel https://github.com/jolespin/veba.git

# Update the permissions
chmod 775 veba/bin/*.py
chmod 775 veba/bin/scripts/*
chmod 775 veba/install/*.sh

# Go into the install directory
cd veba/install
``` 

**2. Install VEBA environments**

**Recommended resource allocatation:** 4 hours with 16 GB memory (include extra time for variable I/O speed for various hosts)

The update from `CheckM1` -> `CheckM2` and installation of `antiSMASH` require more memory and may require grid access if head node is limited.

```
bash install.sh
```

If you want to specify a certain log directory: 
```
bash install.sh path/to/log_directory/ # Default: logs/veba_installation/
```

If you need custom locations for you `conda` environments: # Yes, you need to use the log positional argument too
```
bash install.sh path/to/log_directory/ path/to/conda_environments_directory/
```

**3. Activate the database conda environment, download, and configure databases**

**Recommended resource allocatation:**  48 GB memory (time is dependent on I/O of database repositories)

The databases are installed with a series scripts to allow for custom installations and database builds.  The easiest way to download and configure the databases is to just run main: `veba/install/download_databases.sh` which calls the 4 database scripts and the environment variable script. See FAQ for more details.  

⚠️ **This step should use ~48 GB memory** and should be run using a compute grid via `SLURM` or `SunGridEngine`.  **If this command is run on the head node it will likely fail or timeout if a connection is interrupted.** The most computationally intensive steps are creating a `Diamond` database of `UniRef` and a `MMSEQS2` database of the `MicroEuk100/90/50`.  The longest step is downloading and decompressing `GTDB`.  

The duration will depend on several factors including your internet connection speed and the I/O of public repositories.

Total size is `~300 GB` but if you have certain databases installed already then you can just symlink them so the `VEBA_DATABASE` path has the correct structure.  The exact size may vary as `Pfam` and `UniRef` are updated regularly. 

If issues arise, please [submit a GitHub issue](https://github.com/jolespin/veba/issues) prefixed with `[Database]`. We are here to help :)

**Running an interactive queue:**

```
conda activate VEBA-database_env

bash download_databases.sh /path/to/veba_database/
```

**Running with a job scheduler (e.g., sbatch or qsub):**

[If you need a refresher on SLURM or SunGridEngine, we got you](https://github.com/jolespin/veba/blob/main/walkthroughs/README.md#basics).  If you've never heard of [SLURM](https://slurm.schedmd.com) or [SunGridEngine (e.g., qsub)](https://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html) then you may need to contact your IT for consulting on how to use job scheduling on your system.

Running `conda activate` on a compute server might prompt you to run `conda init` even if you've already initilized on the head node.  To get around this you can use `source activate ${ENV}` where `${ENV}` is the name or path of the `conda` environment.  

Using the `source activate` command requires you to be in `base` conda environment.  You can do this via `conda deactivate` or `conda activate base` before you submit your job.  

Below is an example of how to do this: 

```
# Activate your base environment
conda activate base

# Set the number of threads you want to use.  
# Keep in mind that not all steps are parallelized (e.g., wget)
N_JOBS=1

# Create a log directory
mkdir -p logs/

# Set name for log files when running on the grid
N="database_config"

# Adapt your command to a one-liner
CMD="source activate VEBA-database_env && bash download_databases.sh /path/to/veba_database/"
	
# Note: You should either use SunGridEngine or SLURM not both. 
# You might need to adapt slightly but these worked on our systems.

# SunGridEngine:
qsub -o logs/${N}.o -e logs/${N}.e -cwd -N ${N} -j y -pe threaded ${N_JOBS} "${CMD}"
	
# SLURM:
# For SLURM you might need to specify which account and partition you are associated with for grid jobs
PARTITION=[partition name]
ACCOUNT=[account name]

sbatch -A ${ACCOUNT} -p ${PARTITION} -J ${N} -N 1 -c ${N_JOBS} --ntasks-per-node=1 -o logs/${N}.o -e logs/${N}.e --export=ALL -t 16:00:00 --mem=24G --wrap="${CMD}"
```

Now, you should have the following environments:

```
VEBA-annotate_env
VEBA-assembly_env
VEBA-binning-eukaryotic_env
VEBA-binning-prokaryotic_env
VEBA-binning-viral_env
VEBA-biosynthetic_env
VEBA-classify-eukaryotic_env
VEBA-classify-prokaryotic_env
VEBA-classify-viral_env
VEBA-cluster_env
VEBA-database_env
VEBA-mapping_env
VEBA-phylogeny_env
VEBA-preprocess_env
VEBA-profile_env
```

All the environments should have the `VEBA_DATABASE` environment variable set. If not, then add it manually to ~/.bash_profile: `export VEBA_DATABASE=/path/to/veba_database`.

You can check to make sure the `conda` environments were created and all of the environment variables were created using the following command:

```
bash check_installation.sh
```

Future versions will have `bioconda` installation available.

#### Updating environment executables/scripts:


Let's say you have a separate VEBA repository you want to install scripts from (e.g., you pulled the developmental from GitHub)

```
bash update_environment_scripts.sh path/to/veba_repository/
```

or if you have your `conda` environments in a separate directory than your base `conda`:

```
bash update_environment_scripts.sh path/to/veba_repository path/to/conda_environments/
```


#### Updating environment variables with pre-configured database:

Let's say you moved your database directory or configured a new one somewhere else.  You can do the following:

```
bash update_environment_variables.sh path/to/veba_database/
```

or if you have your environments in a custom path:

```
bash update_environment_variables.sh path/to/database/ path/to/conda_environments/
```

#### Custom install of a single module and a subset of the database:

See [FAQ for custom installation tutorial](https://github.com/jolespin/veba/blob/main/FAQ.md#how-can-i-install-just-a-single-module-and-a-subset-of-the-database-required-for-that-module). 

____________________________________________________________

### Uninstall:

**There are 2 steps to uninstall *VEBA*:**

* Remove conda environments

* Remove database directory

```
# Remove conda enivronments
bash uninstall.sh

# Remove VEBA database
rm -rfv /path/to/veba_database/
```
____________________________________________________________

### Updating VEBA: 

There are currently 2 ways to update veba:

1. Basic uninstall/reinstall - You can uninstall and reinstall using the scripts in `veba/install/` directory.  It's recommended to do a fresh reinstall when do a major update.
2. Patching existing installation - Please post an issue if you have questions updating a previously installed version with a newer or developmental version. This is more advanced and should be done on a case-by-case basis.  

____________________________________________________________

### Install via Docker or Singularity:

If you are using containerized version via `Docker` you will need to make sure it's installed on your system.  If you are using a shared HPC then you probably have access to `Singularity` since `Docker` needs root access which you probably don't have on your user account.  If you have your own machine, then you will probably be using `Docker`. 

Please refer to the [adapting commands for Docker walkthrough](https://github.com/jolespin/veba/blob/main/walkthroughs/adapting_commands_for_docker.md).

Docker containers are now available (starting with `v1.1.2`) for all modules via [DockerHub](https://hub.docker.com/repositories/jolespin)

____________________________________________________________


**Developer Note:** `VEBA ≥ v2.0.0` Docker images is built-on Apple Silicon emulating Linux AMD64 architecture.  If you experience any issues with these images, please post an issue on GitHub. 




