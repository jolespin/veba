<a name="readme-top"></a>

![Maintainer](https://img.shields.io/badge/Maintainer-@jolespin-blue) ![License](https://img.shields.io/badge/License-AGPLv3-blue) ![DOI:10.1186/s12859-022-04973-8](https://zenodo.org/badge/DOI/10.1186/s12859-022-04973-8.svg)

[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]



[forks-shield]: https://img.shields.io/github/forks/jolespin/veba.svg?style=for-the-badge
[forks-url]: https://github.com/jolespin/veba/members
[stars-shield]: https://img.shields.io/github/stars/jolespin/veba.svg?style=for-the-badge
[stars-url]: https://github.com/jolespin/veba/stargazers
[issues-shield]: https://img.shields.io/github/issues/jolespin/veba.svg?style=for-the-badge
[issues-url]: https://github.com/jolespin/veba/issues

```
 _    _ _______ ______  _______
  \  /  |______ |_____] |_____|
   \/   |______ |_____] |     |
```

### What is VEBA? 
The *Viral Eukaryotic Bacterial Archaeal* (VEBA) is an open-source software suite developed with all domains of microorganisms as the primary objective (not post hoc adjustments) including prokaryotic, eukaryotic, and viral organisms.  VEBA is an end-to-end metagenomics and bioprospecting software suite that can directly recover and analyze eukaryotic and viral genomes in addition to prokaryotic genomes with native support for candidate phyla radiation (CPR). VEBA implements a novel iterative binning procedure and an optional hybrid sample-specific/multi-sample framework that recovers more genomes than non-iterative methods.  To optimize the microeukaryotic gene calling and taxonomic classifications, VEBA includes a consensus microeukaryotic database containing protists and fungi compiled from several existing databases. VEBA also provides a unique clustering-based dereplication strategy allowing for sample-specific genomes and proteins to be directly compared across non-overlapping biological samples. VEBA also automates biosynthetic gene cluster identification and novelty scores for bioprospecting.

VEBA's mission is to make robust (meta-)genomics/transcriptomics analysis effortless.  The philosophy of VEBA is that workflows should be modular, generalizable, and easy-to-use with minimal intermediate steps.  The approach implemented in VEBA is to (try and) think 2 steps ahead of what you may need to do and automate the task for you.

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### Citation

Espinoza JL, Dupont CL. VEBA: a modular end-to-end suite for in silico recovery, clustering, and analysis of prokaryotic, microeukaryotic, and viral genomes from metagenomes. BMC Bioinformatics. 2022 Oct 12;23(1):419. [doi: 10.1186/s12859-022-04973-8](https://doi.org/10.1186/s12859-022-04973-8). PMID: 36224545.

Please cite the software dependencies described under the [*Dependency Citation Table*](CITATIONS.md).

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### Announcements

* **What's new in `VEBA v1.4.0`?**

* **`VEBA` Modules:**

	* Added `profile-taxonomic.py` module which uses `sylph` to build a sketch database for genomes and queries the genome database for taxonomic abundance.
	* Added long read support for `fastq_preprocessor`, `preprocess.py`, `assembly-long.py`, `coverage-long`, and all binning modules.
	* Redesign `binning-eukaryotic` module to handle custom `MetaEuk` databases
	* Added new usage syntax `veba --module preprocess --params “${PARAMS}”` where the Conda environment is abstracted and determined automatically in the backend.  Changed all the walkthroughs to reflect this change.
	* Added `skani` which is the new default for genome-level clustering based on ANI.
	* Added `Diamond DeepClust` as an alternative to `MMSEQS2` for protein clustering.

* **`VEBA` Database (`VDB_v6`)**:

	* Completely rebuilt `VEBA's Microeukaryotic Protein Database` to produce a clustered database `MicroEuk100/90/50` similar to `UniRef100/90/50`. Available on [doi:10.5281/zenodo.10139450](https://zenodo.org/records/10139451).


Check out the [*VEBA* Change Log](CHANGELOG.md) for insight into what is being implemented in the upcoming version.

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________


### Installation and databases

**Current Stable Version:** [`v1.4.0`](https://github.com/jolespin/veba/releases/tag/v1.4.0)

**Current Database Version:** `VDB_v6`

Please refer to the [*Installation and Database Configuration Guide*](install/README.md) for software installation and database configuration.

Docker containers are now available (starting with `v1.1.2`) for all modules via [DockerHub](https://hub.docker.com/repositories/jolespin)

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### Getting started with *VEBA*

[*Usage and Resource Requirements Guide*](src/README.md) for parameters and module descriptions

[*Walkthrough Guides*](walkthroughs/README.md) for tutorials and workflows on how to get started

**Usage Example:**

Running `preprocess` module. 

1) Available with `v1.4.0+`:

```
source activate VEBA
veba --module preprocess --params "{PARAMS}" 
```

2) Available with `v1.0.0 - v1.4.0+`:

```
source activate VEBA-preprocess_env
preprocess.py "{PARAMS}"
```

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### What does *VEBA* do?

Please refer to the [*Modules*](src/README.md) for a description of all *VEBA* modules and their functionality.

If you wish *VEBA* did something that isn't implemented, please submit a [`[Feature Request Issue]`](https://github.com/jolespin/veba/issues/new/choose).

[![Schematic](images/Schematic.png)](images/Schematic.pdf)

<p align="right"><a href="#readme-top">^__^</a></p>


___________________________________________________________________


### Output structure
*VEBA*'s is built on the [*GenoPype*](https://github.com/jolespin/genopype) archituecture which creates a reproducible and easy-to-navigate directory structure.  *GenoPype*'s philosophy is to use the same names for all files but to have sample names as subdirectories.  This makes it easier to glob files for grepping, concatenating, etc. *NextFlow* support is in the works...

Example of *GenoPype*'s layout:

```
# Project directory
project_directory/

# Temporary directory
project_directory/tmp/

# Log directory
project_directory/logs/
project_directory/logs/[step]__[program-name].e
project_directory/logs/[step]__[program-name].o
project_directory/logs/[step]__[program-name].returncode

# Checkpoint directory
project_directory/checkpoints/
project_directory/checkpoints/

# Intermediate directories for each step
project_directory/intermediate/
project_directory/intermediate/[step]__[program-name]/

# Output directory
project_directory/output/

# Commands
project_directory/commands.sh
```



For *VEBA*, it has all the directories created by `GenoPype` above but is built for having multiple samples under the same project. 

Example of *VEBA*'s default directory layout:

```
ID="sample_1"

# Main output directory
veba_output/

# Assembly directory
veba_output/assembly

# Assembly output for ${ID} sample
veba_output/assembly/${ID}/output/

# Prokaryotic binning for ${ID} sample
veba_output/binning/prokaryotic/${ID}/output/ 

# Eukaryotic binning
veba_output/binning/eukaryotic/${ID}/output/

# Viral binning
veba_output/binning/viral/${ID}/output/
```

The above are default output locations but they can be customized.


<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### Frequently Asked Questions

If perusing the [*Frequently Asked Questions*](FAQ.md) doesn't address your question, feel free to submit a [[`Question Issue`]](https://github.com/jolespin/veba/issues/new) 

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

