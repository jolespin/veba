## Frequently Asked Questions
<a name="faq-top"></a>

______________________
#### *VEBA* has so many modules and capabilities, how can I get a feel for how to use it for my dataset(s)?

Check out the [walkthroughs](https://github.com/jolespin/veba/tree/main/walkthroughs) where there are step-by-step workflows for different types of data.  For a visual walkthrough of the modules, watch the [Getting started with VEBA](https://www.youtube.com/watch?v=pqrIffWNuug) YouTube video. There are several video tutorials on our [YouTube Channel @VEBA-Multiomics](https://www.youtube.com/@VEBA-Multiomics) covering topics such as how to get started, how to install/configure databases, custom installations/databases, Docker usage on local machines, and the end-to-end walkthrough in real-ish time.

<p align="right"><a href="#faq-top">^__^</a></p>

______________________
#### I already have some of the databases downloaded. Can I use these preexisting databases with *VEBA* instead of having redundant databases?

Yes! Just symlink them so it fits the database structure detailed out [here](https://github.com/jolespin/veba/tree/main/install#database-structure). Large-ish databases you might be able to symlink are GTDB-Tk, CheckV, CheckM2, KOFAM, or Pfam.  If you do this, make sure you have proper read permissions and your databases fall in line with the specifications in the [version notes](https://github.com/jolespin/veba/blob/main/install/README.md#version-notes). *However, if you do this option it will be difficult to diagnose errors so this should only be for advanced users.* 

For example, let's say you have already downloaded Pfam for another tool or analysis.  Simply symlink it as follows: 

```bash
SOURCE="/path/to/source/Pfam-A.hmm.gz"
DATABASE_DIRECTORY="/path/to/veba_database/"
ln -sf ${SOURCE} ${DATABASE_DIRECTORY}/Annotate/Pfam/
```

Since this is more advanced usage, you'll have to go through and comment out the databases you are symlinking in the download scripts.

<p align="right"><a href="#faq-top">^__^</a></p>

______________________


#### How can I install just a single module and a subset of the database required for that module?

This can be done easily with a custom installation.  For example, let's say you want to only use the `annotate.py` module.  You would go to the [module table](https://github.com/jolespin/veba/blob/main/bin/README.md) to see that `annotate.py` module uses the `VEBA-annotate_env` and the `Annotate` database.  Then you would install the custom build as follows:

```bash
# Download the release
# Follow instructions here: https://github.com/jolespin/veba/tree/main/install

# Specify environment
ENV_NAME="VEBA-annotate_env"

# Create a new environment with the required dependencies
mamba env create -n ${ENV_NAME} -f veba/install/environments/${ENV_NAME}.yml

# Update the scripts in the environments $PATH
bash veba/install/update_environment_scripts.sh veba/

# Configure the annotation database (this is going to run Diamond so you need at least 48GB of memory here)
bash veba/install/download_databases-annotate.sh /path/to/veba_database/
```

You should end up with a directory with the annotation database files and placeholders for the other directories:

```
|-- ACCESS_DATE
|-- Annotate
|   |-- CAZy
|   |   `-- CAZyDB.07262023.dmnd
|   |-- KOFAM
|   |   |-- ko_list
|   |   `-- profiles
|   |-- MIBiG
|   |   `-- mibig_v3.1.dmnd
|   |-- MicrobeAnnotator-KEGG
|   |   |-- KEGG_Bifurcating_Module_Information.pkl
|   |   |-- KEGG_Bifurcating_Module_Information.pkl.md5
|   |   |-- KEGG_Module_Information.txt
|   |   |-- KEGG_Module_Information.txt.md5
|   |   |-- KEGG_Regular_Module_Information.pkl
|   |   |-- KEGG_Regular_Module_Information.pkl.md5
|   |   |-- KEGG_Structural_Module_Information.pkl
|   |   `-- KEGG_Structural_Module_Information.pkl.md5
|   |-- NCBIfam-AMRFinder
|   |   |-- NCBIfam-AMRFinder.changelog.txt
|   |   |-- NCBIfam-AMRFinder.hmm.gz
|   |   `-- NCBIfam-AMRFinder.tsv
|   |-- Pfam
|   |   |-- Pfam-A.hmm.gz
|   |   `-- relnotes.txt
|   |-- UniRef
|   |   |-- uniref50.dmnd
|   |   |-- uniref50.release_note
|   |   |-- uniref90.dmnd
|   |   `-- uniref90.release_note
|   `-- VFDB
|       |-- VFDB_setA_pro.dmnd
|       `-- VFs.xls.gz
|-- Classify
|-- Contamination

```

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### How can I use Docker or Singularity to run VEBA?

Check out the [*VEBA* walkthroughs for Docker, Singularity, and AWS](https://github.com/jolespin/veba/tree/main/walkthroughs#containerization-and-aws).

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### I already have genomes binned and/or genes modeled from another program or downloaded from a repository (e.g., NCBI), can I use them with *VEBA*?

Yes! *VEBA* isn't restrictive with the source of the data for most of the modules.  If you have genomes or gene models derived from another source, you can still use the following modules: `coverage.py`, `cluster.py`, `annotate.py`, `phylogeny.py`, `index.py`, `mapping.py`, and any of the [utility scripts](https://github.com/jolespin/veba/tree/main/src/scripts) that apply. 

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### How can I speed up the installation?

The databases take ~4.5 to download/configure.  **Please refer to the [documentation](https://github.com/jolespin/veba/blob/main/install/README.md)** to make sure you allocate enough resources to run `Diamond` and `MMSEQS2` in the backend of the database config.  If you have connection issues to your remote server, you can always use a screen so it doesn't lose your progress when you are installing VEBA (I tend to do this for large copy jobs).  Here is an example command to ssh into your remote server and launching a screen: `ssh -t [username]@[domain] 'screen -DR'`

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### Are there any database versions that are mandatory?

Yes, there a few and they are detailed out in the [version notes](https://github.com/jolespin/veba/blob/main/install/README.md#version-notes).  The most notable would be `GTDB` which is specific to different `GTDB-Tk` versions and prebuilt mash screens *VEBA* provides.

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### Why are there different conda environments for different modules and how do I know which one to use?

This is because there are SO MANY packages throughout all the workflows that it's literally impossible to install all of them in one environment.  I tried to make the environments as straight forward as possible but I understand this could be confusing so I'm actively working on this.  The environment names are pretty straight forward (e.g., use `VEBA-annotate_env` for the `annotate.py` module) but if you have questions, they are listed out [here](https://github.com/jolespin/veba/blob/main/src/README.md).

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### What is [GenoPype](https://github.com/jolespin/genopype) and why does *VEBA* use it instead of Snakemake or NextFlow?

`GenoPype` is a solution I developed to meet the needs for my personal pipelines.  It creates checkpoints, log files, intermediate directories, validates i/o, and everything under the sun.  Future versions may use `Nextflow` but right now `GenoPype` was designed specifically for *VEBA* and since its inception I've used it as the framework for many production-level pipelines.

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### Can I install this via Bioconda?

Currently, not directly but the install scripts are all built around conda so you are essentially doing the same thing.  However, I will work on getting these up on bioconda soon.

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### Why does `preprocess.py` use [fastq_preprocessor](https://github.com/jolespin/fastq_preprocessor) instead of *KneadData*?

`KneadData` is great and I've used it for years but it's a bit dated at the moment.  There are better tools available for the backend and I basically reimplemented the `KneadData` workflow using the following: `fastp` instead of `Trimmomatic`; automatic repairing with bbsuite's `repair.sh` (necessary for SPAdes assemblers); still uses `bowtie2`; bbsuite's `bbduk.sh` to quantify reads that match k-mers (e.g., ribosomal reads); and runs `seqkit stats` on all the steps for a full accounting of reads at the end.

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### How do I report an issue or contribute?

*VEBA* is currently under active development. If you are interested in requesting features, have questions, or wish to report a bug, please post a GitHub issue prefixed with the tag `[Feature Request]`, `[Question]`, and `[Bug]`, respectively.  If you want to contribute or have any other inquiries, contact me at `jol.espinoz[A|T]gmail[DOT]com`

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### Why did I get a `KeyError: 'TMPDIR'`? ⚠️

This is because some programs can't handle long directory paths.  By default, the temporary directory is set to the TMPDIR environment variable.  If you don't have a TMPDIR environment variable for some reason, add a TMPDIR environment variable to your path either in the script or your ~/.bash_profile.  For example, `export TMPDIR=/path/to/temporary/directory/with/read/write/access.  

Here's information about the canonical `TMPDIR` environment variable:

>`TMPDIR` is the canonical environment variable in Unix and POSIX that should be used to specify a temporary directory for scratch space. Most Unix programs will honor this setting and use its value to denote the scratch area for temporary files instead of the common default of /tmp or /var/tmp.
>
>Source - https://en.wikipedia.org/wiki/TMPDIR

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### Why did I get an `AssertionError` The following path does not exist /path/to/scaffolds\_to\_bins.tsv?

This means that you don't have any MAGs that meet the quality threshold. This is typically an empty file that throws the error.  You could always lower the completeness or completion thresholds but this may yield lower quality results.

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### Why am I getting errors for `Kingfisher`?
 
Common `kingfisher` errors: 

* If it says your missing a command, this is likely because `sra-tools` was not installed properly or at all.  Installing `sra-tools` should fix the problem: `conda install -c bioconda sra-tools --force-reinstall` from the `VEBA-preprocess_env` environment. 

* If you get an error related to `prefetch` try changing the `-m` argument (e.g., `-m aws-http`): 

```
(VEBA-preprocess_env) [jespinoz@exp-15-01 Fastq]$ for ID in $(cat ../identifiers.list); do kingfisher get -r $ID -m prefetch; done
09/18/2022 04:04:25 PM INFO: Attempting download method prefetch ..
09/18/2022 04:04:25 PM WARNING: Method prefetch failed: Error was: Command prefetch -o SRR4114636.sra SRR4114636 returned non-zero exit status 127.
STDERR was: b'bash: prefetch: command not found\n'STDOUT was: b''
09/18/2022 04:04:25 PM WARNING: Method prefetch failed
Traceback (most recent call last):
  File "/expanse/projects/jcl110/anaconda3/envs/VEBA-preprocess_env/bin/kingfisher", line 261, in <module>
    main()
  File "/expanse/projects/jcl110/anaconda3/envs/VEBA-preprocess_env/bin/kingfisher", line 241, in main
    extraction_threads = args.extraction_threads,
  File "/expanse/projects/jcl110/anaconda3/envs/VEBA-preprocess_env/lib/python3.7/site-packages/kingfisher/__init__.py", line 234, in download_and_extract
    raise Exception("No more specified download methods, cannot continue")
Exception: No more specified download methods, cannot continue
```

* If SRA-Tools didn't install correctly, you may get this error when converting .sra to .fastq[.gz] files.  If so, just reinstall `sra-tools` via `conda install -c bioconda sra-tools --force-reinstall` in your `VEBA-preprocess_env` environment: 

```
STDERR was: b'bash: fasterq-dump: command not found\n'STDOUT was: b''
```

<p align="right"><a href="#faq-top">^__^</a></p>

______________________


#### Why am I getting an error at the last step of `binning-prokaryotic.py`? 

You might not actually have any high quality bins.  To check this, manually inspect the `CheckM2` results from the intermediate results.  For example: 

```bash
ID="*" # Change this to the ID you are curious about
for FP in veba_output/binning/prokaryotic/${ID}/intermediate/*__checkm2/filtered/checkm2_output.filtered.tsv; do 
	echo ${FP}
	cat ${FP}
	done
```

Some of these files might (and should) be empty and that's expected but what you're looking for are the actual results.  

If you don't have any `checkm2_output.filtered.tsv` then likely no MAGs passed `DAS Tool` then you probably have very low quality genomes if any. Next step is to check the `CheckM2` output before filtered.

```bash
ID="*" # Change this to the ID you are curious about
for FP in veba_output/binning/prokaryotic/${ID}/intermediate/*__checkm2/quality_report.tsv; do 
	echo ${FP}
	cat ${FP}
	done
```
Are there any MAGs here? If so, how are the completeness values? What about the contamination values? Are they meeting the thresholds? If so, then submit a GitHub issue because they should pass.  If not, then you probably just have poor quality data.  If all of your samples are like this then consider doing a *bona fide* coassembly (not pseudo-coassembly).  [Here is a walkthrough to do that with VEBA.](https://github.com/jolespin/veba/blob/main/walkthroughs/setting_up_coassemblies.md). If that still doesn't yield results then assembly-centric metagenomics is likely not the way forward with your dataset and you should consider using a read-based profiling tool like [Kraken2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) or [MetaPhlAn 4](https://huttenhower.sph.harvard.edu/metaphlan/), both of which are not implemented in *VEBA* but very easily be  used with *VEBA* intermediate files.  If you use an external profiling tool, [you can still use the preprocessed reads from *VEBA*](https://github.com/jolespin/veba/blob/main/walkthroughs/download_and_preprocess_reads.md#4-perform-qualityadapter-trimming-remove-human-contamination-and-count-the-ribosomal-reads-but-dont-remove-them).
	
<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### Why am I getting compatibility issues when creating my environments?

This is likely the result of your channel priority and channel order. [Here's an explanation on the .condarc config file.](https://conda.io/projects/conda/en/latest/user-guide/configuration/use-condarc.html).  Alternatively, just copy the below to the following path: `~/.condarc`

```
channel_priority: flexible
channels:
  - conda-forge
  - defaults
  - bioconda
  - jolespin
report_errors: true
```
<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### How can I use options for backend programs that are not arguments in *VEBA*?

While *VEBA* accounts for the most important parameters in the backend, it doesn't hard-code direct access to every parameter (that would be crazy!).  However, many of these options are still accessible for key programs in the backend.  For instance, if you wanted to adjust the *DAS Tool* score threshold, which is not an argument that is hard-coded in *VEBA* like `--dastool_searchengine`, you can simple add the following argument in your *VEBA* command: `--dastool_options='--score_threshold 0.6 --megabin_penalty 0.7'`.  This can be any number of additional arguments, just note that certain ones can break *VEBA* (e.g., changing basenames) so be mindful. 


**^ Please note the usage of quotes here and the equal sign^** 

When using this functionality, just make sure that the argument doesn't overlap with the specified arguments for *VEBA*.  For instance, in the case of *DAS Tool* we already hard-coded access to the `--search_engine` argument via the `--dastool_searchengine` so don't use `--dastool_options '--search_engine <value>'`. Again, be mindful when using this advanced usage.

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### I got an error, how can I diagnose the issue with the internal logs from each step?

*VEBA* is set up to produce log files for each step.  Most steps cannot proceed with an error but to allow for convergence on iterative binning errors are allowed in some steps of the `binning-prokaryotic.py` module.  If you are experiencing an error, look at the log files in the project directory.  For instance, if you are recieving an error for `binning-prokaryotic.py` look under the following directory: `veba_output/binning/prokaryotic/${ID}/logs/` where stderr and stdout are denoted by `.e` and `.o` extensions.  Also, check out the files in the corresponding intermediate directory.  

For instance, if you received an error during `binning-prokaryotic.py` then look at these files to diagnose your issues. 

  * Did any MAGs that made it pass the filters?**

  `cat veba_output/binning/prokaryotic/${ID}/intermediate/*__checkm2/filtered/checkm2_results.filtered.tsv`

  If so, then you should check the last step.  If you have 10 iterations then it will be step 63.  If you have fewer iterations, then it will be a different step that is lower.

  If not, then manually inspect the `CheckM2` results before filtering. 

  `cat veba_output/binning/prokaryotic/${ID}/intermediate/*__checkm2/quality_report.tsv`

  * Do you have MAGs there? Do any of them look legit or are they poor quality? **If your MAGs are `≥ the --checkm_completeness` and `< the --checkm_contamination` thresholds but are not making it through the step**, then please submit a GitHub issue with your log files, scaffolds, and BAM file so I can reproduce and diagnose.

  Work backwards, do you see anything in `7__dastool`? If not, were there are any bins in steps 3-6? 


If you can't figure it out, then submit a GitHub issue ticket and provide a zipped directory of the log files. 

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### How can I restart a module from a specific step?

You can do this by using the `--restart_from_checkpoint <int>` argument which is available on all of the modules.  This goes through and removes all of the checkpoints and intermediate files from that step onwards. 

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### What can I do if `MaxBin2` or `CONCOCT` is taking magnitudes longer to run than `Metabat2` in `binning-prokaryotic.py` module?

If you have a lot of samples and a lot of contigs then `MaxBin2` is likely taking forever to run.  If this is the case, you can use the `--skip_maxbin2` flag because it takes MUCH longer to run. For the Plastisphere it was going to take 40 hours per `MaxBin2` run (there are 2 `MaxBin2` runs) per iteration. `Metabat2` and `CONCOCT` can do the heavy lifting much faster and often with better results so it's recommended to skip `MaxBin2` for larger datasets.  In cases where you have a lot of sample, `CONCOCT` can take a long time so you may want to use `--skip_concoct`. 

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### The host for the microbiome I'm studying isn't human, but instead *[organism X]*.  How can I remove host contamination?

You can either make your own database or, if you are studying a model organism, you can download the precompiled index files for your organism on the [Bowtie2 website](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml). 

Here are a few shortcuts: 

* [*H. sapiens* CHM13 v2 (T2T)](https://genome-idx.s3.amazonaws.com/bt/chm13v2.0.zip)
* [*M. musculus* GRCm39](https://genome-idx.s3.amazonaws.com/bt/GRCm39.zip)
* [*A. thaliana* TAIR10](https://genome-idx.s3.amazonaws.com/bt/TAIR10.zip)

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### I'm using long reads, how can I remove host contamination?

If you are using human microbiomes, you can uncomment the minimap2 human T2T build when configuring your databases.  If you're already run the database configuration, simply copy those commands and run them manually to configure the database.  

```bash
DATABASE_DIRECTORY="path/to/veba_database/"
wget -v -P ${DATABASE_DIRECTORY} https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
minimap2 -d ${DATABASE_DIRECTORY}/Contamination/chm13v2.0/chm13v2.0.mmi ${DATABASE_DIRECTORY}/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
rm -rf ${DATABASE_DIRECTORY}/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
```

If you have a different organism or genome build, simply use the minimap2 command with your genome build.


<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### What's the difference between a co-assembly and a pseudo-coassembly?

Coassembly is when a user concatenates all forward reads into one file (e.g., `cat *_1.fastq.gz > concat_1.fastq.gz`) and all reverse reads into another file (e.g., `cat *_2.fastq.gz > concat_2.fastq.gz`) which is then input into an assembly algorithm (e.g., `metaSPAdes`) to perform "coassembly".  This is often performed when the samples are similar enough to contain similar strains of bacteria and the samples are not deep enough to yield high quality sample-specific assemblies. 

For pseudo-coassembly binning, the user first assembles all of the samples individually (i.e., sample-specific assembly) and then bins out MAGs; preferably using iterative prokaryotic binning followed by eukaryotic and viral binning if applicable.  In most pipelines, the unbinned contigs are discarded but in certain cases (e.g., when the samples are similar enough in origin such as different samples from the same location in the same study) these unbinned contigs can repurposed in a "pseudo-coassembly", a concept introduced in the *VEBA* publication, where unbinned contigs are concatenated together to produce a pseudo-coassembly (e.g., `cat */unbinned.fasta > pseudo-coassembly.fasta`).  **Note that an additional round of assembly is not performed here.** The logic for this procedure is that genomes left over after binning in each individual sample are incomplete fragments which is why they were not recovered during the sample-specific binning and pseudo-coassembly binning has the potential to combine said fragments into a complete genome with reduced likelihood of contaminated genomes than binning using the entire coassembled dataset.  

For more information on *bona fide* coassemblies and what they are, please refer to [AstrobioMike's Happy Belly Bioinformatics blogpost](https://astrobiomike.github.io/metagenomics/metagen_anvio#what-is-a-co-assembly).

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### What's the difference between a bin and a MAG?

In the *VEBA* suite, we define bins as candidate genomes output by binning algorithms that have not been quality assessed and MAGs as genomes that have been quality filtered by *CheckM*, *BUSCO*, and *CheckV* for prokaryotes, eukaryotes, and viruses, respectively.

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### Error when installing environments through `conda` (or `mamba`) saying `Encountered problems while solving` and/or `Problem: nothing provides`?

This is a common issue with `conda` (and `mamba`) and can usually be solved with 2 steps.  

1) First and foremost, make sure you have the most recent version of `conda` (or `mamba`) installed via `conda update conda` (and `conda update mamba`, respectively). [This issue has been well documented on QIIME2's forum.](https://forum.qiime2.org/t/installing-qiime2-with-mamba/21911/4)

2) The second action you can do is set your channel priorities in your `~/.condarc` (if you don't have one, then created one). 

```
(base) cat ~/.condarc
channel_priority: flexible
channels:
  - conda-forge
  - bioconda
  - jolespin
  - defaults
  - qiime2/label/r2022.2

report_errors: true
```
<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### How can I make minor updates instead of reinstalling everything?

Use the `update_environment_scripts.sh` or `update_environment_variables.sh` scripts that are in `veba/install/`.  You may want to do this if there is a minor update to a script that doesn't change the versions or required databases.

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### While running `assembly.py` my job errored (e.g., timed out or ran out of memory), how can I resume the assembly without starting over?

If you're using a `SPAdes`-based program (e.g., `metaSPAdes`) you can use one of two options: 

* `--assembler_options='--continue'` - Continue run from the last available check-point

* `--assembler_options='--restart-from last'` - Restart run with updated options and from the specified check-point

For example, the following `assembly.py` command: `source activate VEBA-assembly_env && assembly.py -1 ${R1} -2 ${R2} -n ${ID} -o ${OUT_DIR} -p ${N_JOBS} --assembler_options='--continue'`

*VEBA* handles these edge case options and removes the other arguments.

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### My job errored in the middle of a step because of a minor issue, how can I continue from the middle of a step and created a checkpoint?

Let's say you were running `assembly.py` and the module failed right after the computationally expensive assembly because of a missing Python library in `fasta_to_saf.py`.   You obviously don't want to recompute the assembly and your step is almost complete as the remaining commands will take a few seconds.

Here's the original command that failed because of a minor error: 

```
assembly.py -1 ${R1} -2 ${R2} -n ${ID} -o ${OUT_DIR} -p ${N_JOBS}
```

1. First check out the `commands.sh` file in the subdirectory for your sample. 


```
(base) head -n 2 veba_output/assembly/SRR5720219/commands.sh
# 1__assembly
( /expanse/projects/jcl110/anaconda3/envs/VEBA-assembly_env/bin/metaspades.py -o veba_output/assembly/SRR5720219/intermediate/1__assembly --tmp-dir veba_output/assembly/SRR5720219/tmp/assembly --threads 4 --memory 250 --restart-from last ) && python /expanse/projects/jcl110/anaconda3/envs/VEBA-assembly_env/bin/scripts/fasta_to_saf.py -i veba_output/assembly/SRR5720219/intermediate/1__assembly/scaffolds.fasta > veba_output/assembly/SRR5720219/intermediate/1__assembly/scaffolds.fasta.saf && rm -rf veba_output/assembly/SRR5720219/intermediate/1__assembly/before_rr.fasta && rm -rf veba_output/assembly/SRR5720219/intermediate/1__assembly/K* && rm -rf veba_output/assembly/SRR5720219/intermediate/1__assembly/misc && rm -rf veba_output/assembly/SRR5720219/intermediate/1__assembly/corrected && rm -rf veba_output/assembly/SRR5720219/intermediate/1__assembly/first_pe_contigs.fasta
```

2. Run the reminaing steps manually.  

```
python /expanse/projects/jcl110/anaconda3/envs/VEBA-assembly_env/bin/scripts/fasta_to_saf.py -i veba_output/assembly/SRR5720219/intermediate/1__assembly/scaffolds.fasta > veba_output/assembly/SRR5720219/intermediate/1__assembly/scaffolds.fasta.saf 

rm -rf veba_output/assembly/SRR5720219/intermediate/1__assembly/before_rr.fasta && rm -rf veba_output/assembly/SRR5720219/intermediate/1__assembly/K* && rm -rf veba_output/assembly/SRR5720219/intermediate/1__assembly/misc && rm -rf veba_output/assembly/SRR5720219/intermediate/1__assembly/corrected && rm -rf veba_output/assembly/SRR5720219/intermediate/1__assembly/first_pe_contigs.fasta
```

3. Lastly, create a manual checkpoint.  It can say anything in the text file but I usually add the date: 

```
echo "Manual run: $(date)" > veba_output/assembly/SRR5720219/checkpoints/1__assembly
```

4. Now rerun the original command: `assembly.py -1 ${R1} -2 ${R2} -n ${ID} -o ${OUT_DIR} -p ${N_JOBS}`

*VEBA* should register that step 1 is complete and will continue with step 2. 

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### How can I install a developmental module environment in my *VEBA* installation?

Developmental/experimental environments can be installed separately.  They are not installed automatically because they use far more compute resources and time than the other environments.  These include the following: 

* `VEBA-amplicon_env`
* `VEBA-crispr_env`

The below code shows how to install the `VEBA-crispr_env` environment as an example:

1. Specify the path to the VEBA repository directory: 

```
VEBA_REPOSITORY_DIRECTORY=path/to/veba_repository_directory (e.g., a release or from git clone https://github.com/jolespin/veba)
```

2. Create the environment: 

```
conda env create -n VEBA-crispr_env -f veba/install/environments/devel/VEBA-crispr_env.yml
```

3. Add the scripts to the environments: 

```
bash ${VEBA_REPOSITORY_DIRECTORY}/install/update_environment_scripts.sh ${VEBA_REPOSITORY_DIRECTORY}
```

4. Update the environment variables:

```
VEBA_DATABASE=/path/to/veba_database

bash ${VEBA_REPOSITORY_DIRECTORY}/install/update_environment_variables.sh ${VEBA_DATABASE}
```

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### Why are my `unbinned.fasta` files missing or empty?

If you are running `binning-prokaryotic.py` -> `binning-eukaryotic.py` -> `binning-viral.py` or similar you may encounter either a missing `unbinned.fasta` file or an empty `unbinned.fasta` file.  These go hand in hand and are most likely because no contigs met the minimum length requirements leading to an empty file.  If there is an empty `unbinned.fasta` file then trying to use it downstream (e.g., `binning-prokaryotic.py` -> `binning-eukaryotic.py`) will throw an error during file validation and will not produce any `unbinned.fasta` files downstream of this.  You can check this by using the following command: 

```bash
M=1500 # Whatever your minimum contig length is set for when you initially ran
ID= # The ID in question

# Check to see if any sequences are long enough
cat veba_output/assembly/${ID}/output/scaffolds.fasta | seqkit seq -m ${M} | grep -c "^>"

# If the output is 0 then the explanation above holds.
```

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### Does **VEBA** support Trinity for transcriptome assembly?

I've considered adding this but it adds [A LOT of dependencies](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/trinity/meta.yaml) including R dependencies (which usually complicate environments) and post-processing analysis tools which are out-of-scope.  As a result of this, *VEBA* will not support Trinity.  However, you can easily use Trinity-based transcripts with other *VEBA* modules but may need to generate the `mapped.sorted.bam` files yourself.

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### If I have counts (i.e., reads mapped) to a Species-Level Cluster (SLC), does that mean I have a corresponding MAG in the sample?

Yes and no.
 
Let’s say there actually was *organism\_A* in *sample\_1* but *organism\_B* was at much higher abundance. In the sequencing process, you are shredding up DNA and sequencing those fragmented bits.  That means organisms that are in higher abundance will take up more of the slots available that are going to be sequenced.  The total number of slots would be analogous to sequencing depth. So in the end you’re not getting a true measure of absolute abundances but a sample of relative abundances.
 
It’s possible that not all of *organism\_A* got sequenced because more abundant organisms like *organism\_B* took up more of the available slots available by the sequencer.   In this case, it is possible that not enough reads were sampled from *organism\_A* to produce long enough of contigs or contigs that had the marker genes so *sample\_1* wasn’t able to yield any high quality MAG for that particular organism (i.e., *organism\_A*).  
 
Then consider *sample\_2* where the abundance and coverage was high enough for *organism\_A* that it could produce long contigs and these long contigs had enough marker genes to be considered high quality; this could yield the fully assembled MAG.
 
When we map the reads back using global mapping, we are mapping to ALL the MAGs not just the MAGs from the corresponding sample. That means the reads from the low abundance *organism\_A* in *sample\_1* would have reads that mapped to *organism\_A* from the *sample\_2* MAG; even though we were not able to recover a complete *organism\_A* MAG from *sample\_1* it could still be in there, albeit, fragmented.

If this is NOT what you want, then use local mapping mode instead of global mapping.

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

####  Why can't I run multiple instances of `transdecoder_wrapper.py` at the same time?

In short, the tool forces the output files to be in the current working directory (even though an output directory is specified).  I've submitted a feature request issue on GitHub (https://github.com/TransDecoder/TransDecoder/issues/169).  I considered forking and making an unofficial update but I don't know Perl and further updating *VEBA* takes priority.

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### Why am I getting a (core dumped) error for `annotate.py` when running `hmmsearch`?

```
cat annotation_output/log/2__hmmsearch-pfam.e

Fatal exception (source file p7_pipeline.c, line 697):
Target sequence length > 100K, over comparison pipeline limit.
(Did you mean to use nhmmer/nhmmscan?)
/bin/sh: line 1: 3527164 Aborted                 (core dumped) ( /expanse/projects/jcl110/anaconda3/envs/VEBA-annotate_env/bin/hmmsearch --tblout annotation_output/intermediate/2__hmmsearch-pfam/output.tsv --cut_ga --cpu 8 --seed 1 ${VEBA_DATABASE}/Annotate/Pfam/Pfam-A.hmm.gz proteins.faa > /dev/null )
```

This is likely because you have [sequences longer than 100k](https://www.biostars.org/p/487110/).  In versions after `v1.1.0` this will be addressed in the backend but in the meantime you can do the following to not trigger this error: `seqkit seq -M 100000 proteins.faa > proteins.lt100k.faa` (assuming your fasta file is called `proteins.faa`).

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### I get an error when trying to use custom options (e.g., `--assembler_options`) saying it expected one argument even though one was given.

Suppose your `assembly.py` job was prematurely canceled for some reason and you tried to use `--assembler_options '--continue'` to continue where the assembler left off.  

You would get the following error:

```
assembly.py: error: argument --assembler_options: expected one argument
```

To get around this, use an equal sign when providing the argument values.  (i.e., `--assembler_options='--continue'`)

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### I get the following error when running `featureCounts`:  `ERROR: Paired-end reads were detected in single-end read library`

This is happens with the update of `subread v2.0.1 -> v2.0.3` [issue/22](https://github.com/jolespin/veba/issues/22).  For v1.1.1, `binning-viral.py` uses `v2.0.3` before I realized they changed the functionality. The workaround was to use `--featurecounts_options='-p --countReadPairs'`.  In v1.1.2, `subread` has been updated to `v2.0.3` in `VEBA-assembly_env, VEBA-binning-*_env, and VEBA-mapping_env` which uses `-p --countReadPairs` flags as default and bypasses it if `--long_reads` flag is used. Read [this BioStars post](https://www.biostars.org/p/9561574/#9561663) for more information.

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### How can I reinstall just a single module or environment? 

Perhaps you customized your environment and broke it or it just never installed correctly and you're just noticing it now.  Regardless, it's pretty easy to patch your installation. 

Let's say you broke your assembly environment, all you have to do is the following: 

```bash
# Specify environment
ENV_NAME="VEBA-annotate_env"
# Remove the current environment
mamba env remove -n ${ENV_NAME}
# Create a new environment with the required dependencies
mamba env create -n ${ENV_NAME} -f veba/install/environments/${ENV_NAME}.yml
# Update the scripts in the environments $PATH
bash veba/install/update_environment_scripts.sh veba/
# Update the environment variables in the environments $PATH
bash veba/install/update_environment_variables.sh /path/to/veba_database/
```

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### Why did the `DAS_Tool` step of `binning-prokaryotic.py` fail?

Most likely your bins didn't have any marker genes so they failed QC.  To confirm this, check the log files: 


Notice here there were no queries aligned:
```
Total time = 1.466s
Reported 0 pairwise alignments, 0 HSPs.
0 queries aligned.
```

You may also see this error that files are empty:

```
Error: Error detecting input file format. First line seems to be blank.
verifying blast did not work
mv: cannot stat 'veba_output/binning/prokaryotic/22_UDP0074_S22_L003/intermediate/2__prodigal/gene_models.faa.scg': No such file or directory
single copy gene prediction using diamond failed. Aborting
```

If you think this is an error, take a look at your assembly quality: 

`cat [fasta] | seqkit seq -m [minimum_threshold] -a`

Are there any large contigs? What's the N50? 

<p align="right"><a href="#faq-top">^__^</a></p>

______________________


#### How was the MicroEuk_v3 database compiled from the source databases?

Check out the [*VEBA* step-by-step guide](https://github.com/jolespin/veba/tree/main/data/MicroEuk_v3/README.md) on how the MicroEuk_v3 protein database was generated.

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

#### How can I update the database from VEBA v2.1.0 (VEBA Database: VDB_v6) to VEBA v2.2.0 (VEBA Database: VDB_v7)?

After uninstalling and installing the new environments: 

```
bash veba/uninstall.sh
bash veba/install.sh
```

Next you need to update the database. In `GTDB-Tk v2.3.x` the database used was `r214.1` but with `GTDB-Tk v2.4.x` the database used is `r220`.  This update is essential because it swaps out `FastANI` for `Skani` which is already used throughout `VEBA`. We have also changed some directory names and added a few more essential files. To do this, you can follow these steps: 

```bash
# 1. Set VEBA database directory to the `DATABASE_DIRECTORY` variable in an interactive session or a script
DATABASE_DIRECTORY=/path/to/veba_database

# 2. Clear out the existing `GTDB` directory. We could use a wildcard but a single misplaced space can have drastic consequences.  Better to be safe with splitting into 2 lines.
rm -rfv ${DATABASE_DIRECTORY}/Classify/GTDB/
mkdir -p ${DATABASE_DIRECTORY}/Classify/GTDB/

# 3. Download GTDB from the data.ace.uq.edu.au mirror b/c it's way faster than data.gtdb.ecogenomic.org.
GTDB_VERSION="220"
wget -v -P ${DATABASE_DIRECTORY} https://data.ace.uq.edu.au/public/gtdb/data/releases/release${GTDB_VERSION}/${GTDB_VERSION}.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r${GTDB_VERSION}_data.tar.gz
tar xvzf ${DATABASE_DIRECTORY}/gtdbtk_r${GTDB_VERSION}_data.tar.gz -C ${DATABASE_DIRECTORY}
rm -rf ${DATABASE_DIRECTORY}/Classify/GTDB
mv ${DATABASE_DIRECTORY}/release${GTDB_VERSION} ${DATABASE_DIRECTORY}/Classify/GTDB
echo "r${GTDB_VERSION}" > ${DATABASE_DIRECTORY}/Classify/GTDB/database_version
wget -P ${DATABASE_DIRECTORY}/Classify/GTDB/ https://data.ace.uq.edu.au/public/gtdb/data/releases/release${GTDB_VERSION}/${GTDB_VERSION}.0/RELEASE_NOTES.txt 
rm -rf ${DATABASE_DIRECTORY}/gtdbtk_r${GTDB_VERSION}_data.tar.gz

# 4. Download GTDB r220 mash sketch database
GTDB_ZENODO_RECORD_ID="11494307"
wget -v -O ${DATABASE_DIRECTORY}/gtdb_r${GTDB_VERSION}.msh https://zenodo.org/records/${GTDB_ZENODO_RECORD_ID}/files/gtdb_r${GTDB_VERSION}.msh?download=1
mkdir -p ${DATABASE_DIRECTORY}/Classify/GTDB/mash/
mv ${DATABASE_DIRECTORY}/gtdb_r${GTDB_VERSION}.msh ${DATABASE_DIRECTORY}/Classify/GTDB/mash/gtdb.msh

# 5. Download the Pfam clans
wget -v -P ${DATABASE_DIRECTORY}/Annotate/Pfam http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz

# 6. Change KOFAM to KOfam
mv ${DATABASE_DIRECTORY}/Annotate/KOFAM ${DATABASE_DIRECTORY}/Annotate/KOfam
```

After you have the databases downloaded, then you need to set up your environment variables for the `VEBA-classify-prokaryotic_env`:

```bash
bash veba/install/update_environment_variables.sh /path/to/veba_database
```

For more details on this last step, see [how to update environment variables with pre-configured database](https://github.com/jolespin/veba/tree/main/install#updating-environment-variables-with-pre-configured-database).

<p align="right"><a href="#faq-top">^__^</a></p>

______________________

