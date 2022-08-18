#### Frequently Asked Questions

**1. *VEBA* has so many modules and capabilities, how can I get a feel for how to use it for my dataset(s)?**

Check out the [walkthroughs](https://github.com/jolespin/veba/tree/main/walkthroughs) where there are step-by-step workflows for different types of data.

**2. I already have some of the databases downloaded, can I use these with *VEBA* instead of have redundant databases?**

Yes! Just symlink them so it fits the database structure detailed out [here](https://github.com/jolespin/veba/tree/main/install#database-structure). 

**3. I already have genomes binned and/or genes modeled from another program or downloaded from a repository (e.g., NCBI), can I use them with *VEBA*?**

Yes! *VEBA* isn't restrictive with the source of the data for most of the modules.  If you have genomes or gene models derived from another source, you can still use the following modules: `coverage.py`, `cluster.py`, `annotate.py`, `phylogeny.py`, `index.py`, `mapping.py`, and any of the [utility scripts](https://github.com/jolespin/veba/tree/main/src/scripts) that apply. 

**4. How can I speed up the installation?**

You should be able to replace all of the `conda` references to `mamba` but this hasn't been tested yet.  It's on my list of action items.  Right now, it takes ~1.5 hours to install all the environments and ~4.5 to download/configure the databases.  **Please refer to the [documentation](https://github.com/jolespin/veba/blob/main/install/README.md)** to make sure you allocate enough resources to run `Diamond` and `MMSEQS2` in the backend of the database config.

**5. Are there any database versions that are mandatory?**

Yes, there a few and they are detailed out in the [version notes](https://github.com/jolespin/veba/blob/main/install/README.md#version-notes).

**6. Why are there different conda environments for different modules and how do I know which one to use?**

This is because there are SO MANY packages throughout all the workflows that it's literally impossible to install all of them in one environment.  I tried to make the environments as straight forward as possible but I understand this could be confusing so I'm actively working on this.  The environment names are pretty straight forward (e.g., use `VEBA-annotate_env` for the `annotate.py` module) but if you have questions, they are listed out [here](https://github.com/jolespin/veba/blob/main/src/README.md).

**7. What is [GenoPype](https://github.com/jolespin/) and why does *VEBA* use it instead of Snakemake or NextFlow?** 

`GenoPype` is a solution I developed to meet the needs for my personal pipelines.  It creates checkpoints, log files, intermediate directories, validates i/o, and everything under the sun.  Future versions may use `Snakemake` but right now `GenoPype` was designed specifically for *VEBA* and since its inception I've used it as the framework for many production-level pipelines.

**8. Can I install this via bioconda?**

Currently, not directly but the install scripts are all built around conda so you are essentially doing the same thing.  However, I will work on getting these up on bioconda soon.


**9. Why does `preprocess.py` use [fastq_preprocessor](https://github.com/jolespin/fastq_preprocessor) instead of KneadData?**

`KneadData` is great and I've used it for years but it's a bit dated at the moment.  There are better tools available for the backend and I basically reimplemented the `KneadData` workflow using the following: `fastp` instead of `Trimmomatic`; automatic repairing with bbsuite's `repair.sh` (necessary for SPAdes assemblers); still uses `bowtie2`; bbsuite's `bbduk.sh` to quantify reads that match k-mers (e.g., ribosomal reads); and runs `seqkit stats` on all the steps for a full accounting of reads at the end.



**10. How do I report an issue or contribute?** 

*VEBA* is currently under active development. If you are interested in requesting features or wish to report a bug, please post a GitHub issue prefixed with the tag `[Feature Request]` and `[Bug]`, respectively.  If you want to contribute or have any other inquiries, contact me at `jespinoz[A|T]jcvi[DOT]org`
