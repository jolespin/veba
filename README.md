```
 _    _ _______ ______  _______
  \  /  |______ |_____] |_____|
   \/   |______ |_____] |     |

```
### Description
The *Viral Eukaryotic Bacterial Archaeal* (VEBA) is an open-source software suite developed with all domains of microorganisms as the primary objective (not post hoc adjustments) including prokaryotic, eukaryotic, and viral organisms.  To our knowledge, VEBA is the first end-to-end metagenomics software suite that can directly recover and analyze eukaryotic and viral genomes in addition to prokaryotic genomes with integrated support for CPR. VEBA implements a novel iterative binning procedure and an optional hybrid sample-specific/multi-sample framework that recovers more genomes than non-iterative methods.  To optimize the microeukaryotic gene calling and taxonomic classifications, VEBA includes a consensus microeukaryotic database containing protists and fungi compiled from several existing databases. VEBA also provides a unique clustering-based dereplication strategy allowing for sample-specific genomes and proteins to be directly compared across non-overlapping biological samples.  In addition, VEBA is the only pipeline that automates the detection of CPR bacteria and implements the appropriate genome quality assessments for said organisms.   

___________________________________________________________________
### Citation
Espinoza et al. (Accepted)

___________________________________________________________________

### Installation and databases
Please refer to the [*Installation and Database Configuration Guide*](install/README.md) for software installation and database configuration.

**Current Version:** [`v1.0.0`](https://github.com/jolespin/veba/releases/tag/v1.0.0)

___________________________________________________________________
### Getting started with *VEBA*

Please refer to the [*Walkthrough Guides*](walkthroughs/README.md) for tutorials and workflows on how to get started.

___________________________________________________________________

### Frequently Asked Questions

If perusing the [*Frequently Asked Questions*](FAQ.md) doesn't address your question, feel free to submit a [GitHub issue](https://github.com/jolespin/veba/issues/new) with the `[Question]` prefix in the title.


___________________________________________________________________

### *VEBA* Modules

Please refer to the [*Modules*](src/README.md) for a description of all *VEBA* modules and their functionality.

[![Schematic](images/Schematic.png)](images/Schematic.pdf)


* **preprocess** – Fastq quality trimming, adapter removal, decontamination, and read statistics calculations

* **assembly** – Assemble reads, align reads to assembly, and count mapped reads

* **coverage** – Align reads to (concatenated) reference and counts mapped reads

* **binning-prokaryotic** – Iterative consensus binning for recovering prokaryotic genomes with lineage-specific quality assessment

* **binning-eukaryotic** – Binning for recovering eukaryotic genomes with exon-aware gene modeling and lineage-specific quality assessment

* **binning-viral** – Detection of viral genomes and quality assessment

* **classify-prokaryotic** – Taxonomic classification and candidate phyla radiation adjusted quality 

* **classify-eukaryotic** – Taxonomic classification of eukaryotic genomes

* **classify-viral** – Taxonomic classification and isolation source of viral genomes

* **cluster** – Species-level clustering of genomes and lineage-specific orthogroup detection

* **annotate** – Annotates translated gene calls against NR, Pfam, and KOFAM

* **phylogeny** – Constructs phylogenetic trees given a marker set

* **index** – Builds local or global index for alignment to genomes
 
* **mapping** – Aligns reads to local or global index of genomes


___________________________________________________________________

### Output structure
*VEBA*'s is built on the [*GenoPype*](https://github.com/jolespin/genopype) archituecture which creates a reproducible and easy-to-navigate directory structure.  *GenoPype*'s philosophy is to use the same names for all files but to have sample names as subdirectories.  This makes it easier to glob files for grepping, concatenating, etc. 

Here is an example of *GenoPype*'s layout:

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

For *VEBA*, it has all the directories created by `GenoPype` above but is built for having multiple samples under the same project:

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
