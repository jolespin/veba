#### Frequently Asked Questions

**1. *VEBA* has so many modules and capabilities, how can I get a feel for how to use it for my dataset(s)?**

Check out the [walkthroughs](https://github.com/jolespin/veba/tree/main/walkthroughs) where there are step-by-step workflows for different types of data.

**2. It says the total database size is 369G but I already have some of the databases downloaded. Can I use these preexisting databases with *VEBA* instead of having redundant databases?**

Yes! Just symlink them so it fits the database structure detailed out [here](https://github.com/jolespin/veba/tree/main/install#database-structure). The bulk of the database is the `Diamond` database of NCBI's NR proteins which your institute might already have on their servers.  Just make sure they compiled it with taxonomy information like *VEBA* does [here](https://github.com/jolespin/veba/blob/1755c762f3ea5626fb4cbd327b2d24e05dfc0a2f/install/download_databases.sh#L102).  Other large-ish databases you might be able to symlink are [GTDB-Tk](https://github.com/jolespin/veba/blob/1755c762f3ea5626fb4cbd327b2d24e05dfc0a2f/install/download_databases.sh#L32), [CheckV](https://github.com/jolespin/veba/blob/1755c762f3ea5626fb4cbd327b2d24e05dfc0a2f/install/download_databases.sh#L43), [CheckM](https://github.com/jolespin/veba/blob/1755c762f3ea5626fb4cbd327b2d24e05dfc0a2f/install/download_databases.sh#L53), [KOFAM](https://github.com/jolespin/veba/blob/1755c762f3ea5626fb4cbd327b2d24e05dfc0a2f/install/download_databases.sh#L80), or [Pfam](https://github.com/jolespin/veba/blob/1755c762f3ea5626fb4cbd327b2d24e05dfc0a2f/install/download_databases.sh#L89).  Make sure your databases fall in line with the specifications in the [version notes](https://github.com/jolespin/veba/blob/main/install/README.md#version-notes). *However, if you do this option it will be difficult to diagnose errors so this should only be for advanced users.* 

**3. I already have genomes binned and/or genes modeled from another program or downloaded from a repository (e.g., NCBI), can I use them with *VEBA*?**

Yes! *VEBA* isn't restrictive with the source of the data for most of the modules.  If you have genomes or gene models derived from another source, you can still use the following modules: `coverage.py`, `cluster.py`, `annotate.py`, `phylogeny.py`, `index.py`, `mapping.py`, and any of the [utility scripts](https://github.com/jolespin/veba/tree/main/src/scripts) that apply. 

**4. How can I speed up the installation?**

You can replace all of the `conda` references to `mamba` but this hasn't been tested yet.   With `conda`, it takes ~1.5 hours to install all the environments and `mamba` only drops the time by ~10 minutes so it's not recommended. The databases take ~4.5 to download/configure.  **Please refer to the [documentation](https://github.com/jolespin/veba/blob/main/install/README.md)** to make sure you allocate enough resources to run `Diamond` and `MMSEQS2` in the backend of the database config.  If you have connection issues to your remote server, you can always use a screen so it doesn't lose your progress when you are installing VEBA (I tend to do this for large copy jobs).  Here is an example command to ssh into your remote server and launching a screen: `ssh -t [username]@[domain] 'screen -DR'`

**5. Are there any database versions that are mandatory?**

Yes, there a few and they are detailed out in the [version notes](https://github.com/jolespin/veba/blob/main/install/README.md#version-notes).

**6. Why are there different conda environments for different modules and how do I know which one to use?**

This is because there are SO MANY packages throughout all the workflows that it's literally impossible to install all of them in one environment.  I tried to make the environments as straight forward as possible but I understand this could be confusing so I'm actively working on this.  The environment names are pretty straight forward (e.g., use `VEBA-annotate_env` for the `annotate.py` module) but if you have questions, they are listed out [here](https://github.com/jolespin/veba/blob/main/src/README.md).

**7. What is [GenoPype](https://github.com/jolespin/genopype) and why does *VEBA* use it instead of Snakemake or NextFlow?** 

`GenoPype` is a solution I developed to meet the needs for my personal pipelines.  It creates checkpoints, log files, intermediate directories, validates i/o, and everything under the sun.  Future versions may use `Snakemake` but right now `GenoPype` was designed specifically for *VEBA* and since its inception I've used it as the framework for many production-level pipelines.

**8. Can I install this via Bioconda?**

Currently, not directly but the install scripts are all built around conda so you are essentially doing the same thing.  However, I will work on getting these up on bioconda soon.


**9. Why does `preprocess.py` use [fastq_preprocessor](https://github.com/jolespin/fastq_preprocessor) instead of *KneadData*?**

`KneadData` is great and I've used it for years but it's a bit dated at the moment.  There are better tools available for the backend and I basically reimplemented the `KneadData` workflow using the following: `fastp` instead of `Trimmomatic`; automatic repairing with bbsuite's `repair.sh` (necessary for SPAdes assemblers); still uses `bowtie2`; bbsuite's `bbduk.sh` to quantify reads that match k-mers (e.g., ribosomal reads); and runs `seqkit stats` on all the steps for a full accounting of reads at the end.



**10. How do I report an issue or contribute?** 

*VEBA* is currently under active development. If you are interested in requesting features or wish to report a bug, please post a GitHub issue prefixed with the tag `[Feature Request]` and `[Bug]`, respectively.  If you want to contribute or have any other inquiries, contact me at `jespinoz[A|T]jcvi[DOT]org`

**11. During installation, I got *SafetyErrors* and *ClobberErrors*.  Does my *VEBA* installation work?** 

These are known errors that have to do with `CheckM` and `Perl` dependencies, respectively. In short, these are non-fatal errors and will not affect your installation.  For more details, check this section of the [installation manual](https://github.com/jolespin/veba/tree/main/install#common-installation-errors-that-do-not-affect-veba-functionality). 

**12. Why did I get a `KeyError: 'TMPDIR'`?**

This is because CheckM can't handle long directory paths.  By default, the temporary directory is set to the TMPDIR environment variable.  If you don't have a TMPDIR environment variable for some reason, add a TMPDIR environment variable to your path either in the script or your ~/.bash_profile.  For example, `export TMPDIR=/path/to/temporary/directory/with/read/write/access

**13. Why did I get an `AssertionError` The following path does not exist /path/to/scaffolds\_to\_bins.tsv?**

This means that you don't have any MAGs that meet the quality threshold. This is typically an empty file that throws the error.  You could always lower the completeness or completion thresholds but this may yield lower quality results.


