#### Frequently Asked Questions

#### 1. *VEBA* has so many modules and capabilities, how can I get a feel for how to use it for my dataset(s)?

Check out the [walkthroughs](https://github.com/jolespin/veba/tree/main/walkthroughs) where there are step-by-step workflows for different types of data.

#### 2. It says the total database size is 369G but I already have some of the databases downloaded. Can I use these preexisting databases with *VEBA* instead of having redundant databases?

Yes! Just symlink them so it fits the database structure detailed out [here](https://github.com/jolespin/veba/tree/main/install#database-structure). The bulk of the database is the `Diamond` database of NCBI's NR proteins which your institute might already have on their servers.  Just make sure they compiled it with taxonomy information like *VEBA* does [here](https://github.com/jolespin/veba/blob/1755c762f3ea5626fb4cbd327b2d24e05dfc0a2f/install/download_databases.sh#L102).  Other large-ish databases you might be able to symlink are [GTDB-Tk](https://github.com/jolespin/veba/blob/1755c762f3ea5626fb4cbd327b2d24e05dfc0a2f/install/download_databases.sh#L32), [CheckV](https://github.com/jolespin/veba/blob/1755c762f3ea5626fb4cbd327b2d24e05dfc0a2f/install/download_databases.sh#L43), [CheckM](https://github.com/jolespin/veba/blob/1755c762f3ea5626fb4cbd327b2d24e05dfc0a2f/install/download_databases.sh#L53), [KOFAM](https://github.com/jolespin/veba/blob/1755c762f3ea5626fb4cbd327b2d24e05dfc0a2f/install/download_databases.sh#L80), or [Pfam](https://github.com/jolespin/veba/blob/1755c762f3ea5626fb4cbd327b2d24e05dfc0a2f/install/download_databases.sh#L89).  Make sure your databases fall in line with the specifications in the [version notes](https://github.com/jolespin/veba/blob/main/install/README.md#version-notes). *However, if you do this option it will be difficult to diagnose errors so this should only be for advanced users.* 

#### 3. I already have genomes binned and/or genes modeled from another program or downloaded from a repository (e.g., NCBI), can I use them with *VEBA*?

Yes! *VEBA* isn't restrictive with the source of the data for most of the modules.  If you have genomes or gene models derived from another source, you can still use the following modules: `coverage.py`, `cluster.py`, `annotate.py`, `phylogeny.py`, `index.py`, `mapping.py`, and any of the [utility scripts](https://github.com/jolespin/veba/tree/main/src/scripts) that apply. 

#### 4. How can I speed up the installation?

You can replace all of the `conda` references to `mamba` but this hasn't been tested yet.   With `conda`, it takes ~1.5 hours to install all the environments and `mamba` only drops the time by ~10 minutes so it's not recommended. The databases take ~4.5 to download/configure.  **Please refer to the [documentation](https://github.com/jolespin/veba/blob/main/install/README.md)** to make sure you allocate enough resources to run `Diamond` and `MMSEQS2` in the backend of the database config.  If you have connection issues to your remote server, you can always use a screen so it doesn't lose your progress when you are installing VEBA (I tend to do this for large copy jobs).  Here is an example command to ssh into your remote server and launching a screen: `ssh -t [username]@[domain] 'screen -DR'`

#### 5. Are there any database versions that are mandatory?

Yes, there a few and they are detailed out in the [version notes](https://github.com/jolespin/veba/blob/main/install/README.md#version-notes).

#### 6. Why are there different conda environments for different modules and how do I know which one to use?

This is because there are SO MANY packages throughout all the workflows that it's literally impossible to install all of them in one environment.  I tried to make the environments as straight forward as possible but I understand this could be confusing so I'm actively working on this.  The environment names are pretty straight forward (e.g., use `VEBA-annotate_env` for the `annotate.py` module) but if you have questions, they are listed out [here](https://github.com/jolespin/veba/blob/main/src/README.md).

#### 7. What is [GenoPype](https://github.com/jolespin/genopype) and why does *VEBA* use it instead of Snakemake or NextFlow?

`GenoPype` is a solution I developed to meet the needs for my personal pipelines.  It creates checkpoints, log files, intermediate directories, validates i/o, and everything under the sun.  Future versions may use `Snakemake` but right now `GenoPype` was designed specifically for *VEBA* and since its inception I've used it as the framework for many production-level pipelines.

#### 8. Can I install this via Bioconda?

Currently, not directly but the install scripts are all built around conda so you are essentially doing the same thing.  However, I will work on getting these up on bioconda soon.


#### 9. Why does `preprocess.py` use [fastq_preprocessor](https://github.com/jolespin/fastq_preprocessor) instead of *KneadData*?

`KneadData` is great and I've used it for years but it's a bit dated at the moment.  There are better tools available for the backend and I basically reimplemented the `KneadData` workflow using the following: `fastp` instead of `Trimmomatic`; automatic repairing with bbsuite's `repair.sh` (necessary for SPAdes assemblers); still uses `bowtie2`; bbsuite's `bbduk.sh` to quantify reads that match k-mers (e.g., ribosomal reads); and runs `seqkit stats` on all the steps for a full accounting of reads at the end.



#### 10. How do I report an issue or contribute?

*VEBA* is currently under active development. If you are interested in requesting features or wish to report a bug, please post a GitHub issue prefixed with the tag `[Feature Request]` and `[Bug]`, respectively.  If you want to contribute or have any other inquiries, contact me at `jespinoz[A|T]jcvi[DOT]org`

#### 11. During installation, I got *SafetyErrors* and *ClobberErrors*.  Does my *VEBA* installation work?

These are known errors that have to do with `CheckM` and `Perl` dependencies, respectively. In short, these are non-fatal errors and will not affect your installation.  For more details, check this section of the [installation manual](https://github.com/jolespin/veba/tree/main/install#common-installation-errors-that-do-not-affect-veba-functionality). 

#### 12. Why did I get a `KeyError: 'TMPDIR'`?

This is because CheckM can't handle long directory paths.  By default, the temporary directory is set to the TMPDIR environment variable.  If you don't have a TMPDIR environment variable for some reason, add a TMPDIR environment variable to your path either in the script or your ~/.bash_profile.  For example, `export TMPDIR=/path/to/temporary/directory/with/read/write/access.  

Here's information about the canonical `TMPDIR` environment variable:

>`TMPDIR` is the canonical environment variable in Unix and POSIX that should be used to specify a temporary directory for scratch space. Most Unix programs will honor this setting and use its value to denote the scratch area for temporary files instead of the common default of /tmp or /var/tmp.
>
>Source - https://en.wikipedia.org/wiki/TMPDIR

#### 13. Why did I get an `AssertionError` The following path does not exist /path/to/scaffolds\_to\_bins.tsv?

This means that you don't have any MAGs that meet the quality threshold. This is typically an empty file that throws the error.  You could always lower the completeness or completion thresholds but this may yield lower quality results.

#### 14. Why am I getting errors for `Kingfisher` saying a command was not found?

This is likely because `sra-tools` was not installed properly or at all.  Installing `sra-tools` should fix the problem: `conda install -c bioconda sra-tools --force-reinstall` from the `VEBA-preprocess_env` environment.  More details [in this walkthrough](https://github.com/jolespin/veba/blob/main/walkthroughs/download_and_preprocess_reads.md#3-create-a-directory-for-the-raw-fastq-reads-and-download-the-sequences-from-ncbi).

#### 15. Why am I getting an error at the last step of `binning-prokaryotic.py`?

There's a few reasons why you could be getting an error: 

a. You might not actually have any high quality bins.  To check this, manually inspect the `CheckM` results from the intermediate results.  For example: 

```bash
ID="*" # Change this to the ID you are curious about
for FP in veba_output/binning/prokaryotic/${ID}/intermediate/*__checkm/filtered/checkm_output.filtered.tsv; do 
	echo ${FP}
	cat ${FP}
	done
```

Some of these files might (and should) be empty and that's expected but what you're looking for are the actual results.  Recall that to properly handle CPR bacteria we allow `k__Bacteria` phyla through at any completion/contamination levels which will be filtered out at the last step.  `GTDB-Tk` is run on these MAGs and if they are CPR bacteria then they are reevaluated with the [proper marker set](https://github.com/Ecogenomics/CheckM/blob/master/custom_marker_sets/cpr_43_markers.hmm).  If you only have `k__Bacteria` MAGs then there is a chance that you didn't recover any high quality MAGs.

If you don't have any `checkm_output.filtered.tsv` then likely no MAGs passed `DAS Tool` then you probably have very low quality genomes if any. Next step is to check the `CheckM` output before filtered.

```bash
ID="*" # Change this to the ID you are curious about
for FP in veba_output/binning/prokaryotic/${ID}/intermediate/*__checkm/output.tsv; do 
	echo ${FP}
	cat ${FP}
	done
```
Are there any MAGs here? If so, how are the completeness values? What about the contamination values? Are they meeting the thresholds? If so, then submit a GitHub issue because they should pass.  If not, then you probably just have poor quality data.  If all of your samples are like this then consider doing a bona fide coassembly (not pseudo-coassembly).  [Here is a walkthrough to do that with VEBA.](https://github.com/jolespin/veba/blob/main/walkthroughs/setting_up_coassemblies.md). If that still doesn't yield results then assembly-centric metagenomics is likely not the way forward with your dataset and you should consider using a read-based profiling tool like [Kraken2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) which is not implemented in *VEBA* but very easy to use with minimal pre/post-processing.  If you use an external profiling tool, [you can still use the preprocessed reads from *VEBA*](https://github.com/jolespin/veba/blob/main/walkthroughs/download_and_preprocess_reads.md#4-perform-qualityadapter-trimming-remove-human-contamination-and-count-the-ribosomal-reads-but-dont-remove-them).
	

b. If you checked the above and it looks like you have decent MAGs (that is, non `k__Bacteria` hits with acceptable completion and contamination ratios (e.g., completeness ≥ 70, contamination < 10) then check the `GTDB-Tk` log file as you may have an error with `pplacer`.  Check the following log file:

```bash
ID="*" # Change this to the ID you are curious about

cat veba_output/binning/prokaryotic/${ID}/intermediate/*__cpr_adjustment/gtdbtk_output/gtdbtk.log	
```

If you have this error, then it's your `pplacer` installation: 

```
gtdbtk.exceptions.PplacerException: An error was encountered while running pplacer.
```

This is a known issue [GTDBTk/issues/170](https://github.com/Ecogenomics/GTDBTk/issues/170) and is related to memory issues (e.g., try 128GB).  Don't increase the `--pplacer_threads` because this requires massive amounts of memory. 


#### 16. Why am I getting compatibility issues when creating my environments?

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

#### 17. How can I use options for backend programs that are not arguments in *VEBA*?

While *VEBA* accounts for the most important parameters in the backend, it doesn't hard-code direct access to every parameter (that would be crazy!).  However, many of these options are still accessible for key programs in the backend.  For instance, if you wanted to adjust the *DAS Tool* score threshold, which is not an argument that is hard-coded in *VEBA* like `--dastool_searchengine`, you can simple add the following argument in your *VEBA* command: `--dastool_options '--score_threshold 0.6 --megabin_penalty 0.7'`.  This can be any number of additional arguments, just note that certain ones can break *VEBA* (e.g., changing basenames) so be mindful. 

**^ Please note the usage of quotes here ^** 

When using this functionality, just make sure that the argument doesn't overlap with the specified arguments for *VEBA*.  For instance, in the case of *DAS Tool* we already hard-coded access to the `--search_engine` argument via the `--dastool_searchengine` so don't use `--dastool_options '--search_engine <value>'`. Again, be mindful when using this advanced usage.

#### 18. I got an error, how can I diagnose the issue?

*VEBA* is set up to produce log files for each step.  Most steps cannot proceed with an error but to allow for convergence on iterative binning errors are allowed in some steps of the `binning-prokaryotic.py` module.  If you are experiencing an error, look at the log files in the project directory.  For instance, if you are recieving an error for `binning-prokaryotic.py` look under the following directory: `veba_output/binning/prokaryotic/${ID}/logs/` where stderr and stdout are denoted by `.e` and `.o` extensions.  Also, check out the files in the corresponding intermediate directory.  

For instance, if you received an error at `63__cpr_adjustment` then look at these files to diagnose your issues. 

1. Are there any MAGs that made it pass the filters?

`cat veba_output/binning/prokaryotic/${ID}/intermediate/*__checkm/filtered/checkm_output.filtered.tsv`

If so, then you should check the last step.  If you have 10 iterations then it will be step 63.  If you have fewer iterations, then it will be a different step that is lower.

If not, then manually inspect the CheckM results before filtering. 

`cat veba_output/binning/prokaryotic/${ID}/intermediate/*__checkm/output.tsv`

Do you have MAGs there? Do any of them look legit or are they poor quality? **If your MAGs are `≥ the --checkm_completeness` and `< the --checkm_contamination` thresholds but are not making it through the step**, then please submit a GitHub issue with your log files, scaffolds, and BAM file so I can reproduce and diagnose.

See FAQ #15 for more details on this.  

2. Check the stderr/stdout for the last step: 

`cat veba_output/binning/prokaryotic/DE_001/log/63__cpr_adjustment.*`

Did `GTDB-Tk` run or is there an memory issue w/ `pplacer`? What about an unspecified environment variable?  Is that path valid or is incorrect?


If you can't figure it out, then submit a GitHub issue ticket and provide a zipped directory of the log files. 

#### 19. How can I restart a module from a specific step?

You can do this by using the `--restart_from_checkpoint <int>` argument which is available on all of the modules.  This goes through and removes all of the checkpoints and intermediate files from that step onwards. 

#### 20. How can I allocate enough resources for `binning-prokarytic.py` module?

The penultimate step `[step-number]__cpr_adjustment]` of `binning-prokaryotic.py` is the most memory intensive stage as this uses `GTDB-Tk` to identify CPR bacteria. There are two options: 1) The first is the easiest but not the best use of resources for larger datasets.  This option is to use ~243GB memory in the initial command.  This will ensure that `GTBD-Tk` has enough resources to run to completion but if you are being charged for compute resources this easier option may be more expensive.; 2) The second option is more involved is to use a lower amount of memory (e.g., 50GB) in the initial command.  Once the penultimate step is launched, it will fail due to not having enough resources.  If this was step 63 you could rerun the program using `--restart_from_checkpoint 63` and allocate `243G` to your `SLURM/SunGridEngine` job submission.  Personally, this is what I do to minimize resource consumption but this adds an extra manual step.  Again, to get around this you can just allocate the maximum resources at the beginning and not worry about it but the former will be more resource intensive.

Versions prior to `v1.0.2a` (e.g., `v1.0.0` and `v1.0.1` using `GTDB-Tk v1.x` but `VEBA v1.0.2a` has updated to `GTDB-Tk v2.x` which uses MUCH less memory.  If you are using `VEBA v1.0.2a` then you can use ~60GB of memory while earlier versions reliant on `GTDB-Tk v1.x` will require between 128-243GB of memory.

#### 21. What can I do if `MaxBin2` is taking magnitudes longer to run than `Metabat2` and `CONCOCT` in `binning-prokaryotic.py` module?

If you have a lot of samples and a lot of contigs then `MaxBin2` is likely taking forever to run.  If this is the case, you can use the `--skip_maxbin2` flag because it takes MUCH longer to run. For the Plastisphere it was going to take 40 hours per `MaxBin2` run (there are 2 `MaxBin2` runs) per iteration. `Metabat2` and `CONCOCT` can do the heavy lifting much faster and often with better results so it's recommended to skip `MaxBin2` for larger datasets.

#### 22. The host for the microbiome I'm studying isn't human, but instead *[organism X]*.  How can I remove host contamination?

You can either make your own database or, if you are studying a model organism, you can download the precompiled index files for your organism on the [Bowtie2 website](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml). 

Here are a few shortcuts: 

* [*H. sapiens* CHM13 v2 (T2T)](https://genome-idx.s3.amazonaws.com/bt/chm13v2.0.zip)
* [*M. musculus* GRCm39](https://genome-idx.s3.amazonaws.com/bt/GRCm39.zip)
* [*A. thaliana* TAIR10](https://genome-idx.s3.amazonaws.com/bt/TAIR10.zip)

#### 23. What's the difference between a coassembly and a pseudo-coassembly?

Coassembly is when a user concatenates all forward reads into one file (e.g., `cat *_1.fastq.gz > concat_1.fastq.gz`) and all reverse reads into another file (e.g., `cat *_2.fastq.gz > concat_2.fastq.gz`) which is then input into an assembly algorithm (e.g., `metaSPAdes`) to perform "coassembly".  This is often performed when the samples are similar enough to contain similar strains of bacteria and the samples are not deep enough to yield high quality sample-specific assemblies. 

For pseudo-coassembly binning, the user first assembles all of the samples individually (i.e., sample-specific assembly) and then bins out MAGs; preferably using iterative prokaryotic binning followed by eukaryotic and viral binning if applicable.  In most pipelines, the unbinned contigs are discarded but in certain cases (e.g., when the samples are similar enough in origin such as different samples from the same location in the same study) these unbinned contigs can repurposed in a "pseudo-coassembly", a concept introduced in the *VEBA* publication, where unbinned contigs are concatenated together to produce a pseudo-coassembly (e.g., `cat */unbinned.fasta > pseudo-coassembly.fasta`).  **Note that an additional round of assembly is not performed here.** The logic for this procedure is that genomes left over after binning in each individual sample are incomplete fragments which is why they were not recovered during the sample-specific binning and pseudo-coassembly binning has the potential to combine said fragments into a complete genome with reduced likelihood of contaminated genomes than binning using the entire coassembled dataset.  

For more information on *bona fide* coassemblies and what they are, please refer to [AstrobioMike's Happy Belly Bioinformatics blogpost](https://astrobiomike.github.io/metagenomics/metagen_anvio#what-is-a-co-assembly).

#### 24. What's the difference between a bin and a MAG?

In the *VEBA* suite, we define bins as candidate genomes output by binning algorithms that have not been quality assessed and MAGs as genomes that have been quality filtered by *CheckM*, *BUSCO*, and *CheckV* for prokaryotes, eukaryotes, and viruses, respectively.

#### 25. How can I update the human reference included in my *VEBA* database?** (i.e., GRCh38  →  CHM13v2.0)

As of 2022.10.18 *VEBA* has switched from using the "GRCh38 no alt analysis set" to the "CHM13v2.0 telomore-to-telomere" build for human.  If you've installed *VEBA* before this date or are using `v1.0.0` release from [Espinoza et al. 2022](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04973-8) then you can update with the following code:

```
conda activate VEBA-database_env
wget -v -P ${VEBA_DATABASE} https://genome-idx.s3.amazonaws.com/bt/chm13v2.0.zip
unzip -d ${VEBA_DATABASE}/Contamination/ ${VEBA_DATABASE}/chm13v2.0.zip
rm -rf ${VEBA_DATABASE}/chm13v2.0.zip

# Use this if you want to remove the previous GRCh38 index
rm -rf ${VEBA_DATABASE}/Contamination/grch38/
```

#### 26. Error when installing environments through `conda` (or `mamba`) saying `Encountered problems while solving` and/or `Problem: nothing provides`?

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

#### 27. How can I make minor updates instead of reinstalling everything?

Please refer to our [patches](https://github.com/jolespin/veba/blob/main/install/PATCHES.md) section.

#### 28.  While running `assembly.py` my job errored (e.g., timed out or ran out of memory), how can I resume the assembly without starting over?

If you're using a `SPAdes`-based program (e.g., `metaSPAdes`) you can use the following option: `--assembler_options='--restart-from last'`.  

For example, the following `assembly.py` command: `source activate VEBA-assembly_env && assembly.py -1 ${R1} -2 ${R2} -n ${ID} -o ${OUT_DIR} -p ${N_JOBS} --assembler_options='--restart-from last'`

#### 29. My job errored in the middle of a step because of a minor issue, how can I continue from the middle of a step and created a checkpoint?

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

4. Now rerun the original command: `assembly.py -1 ${R1} -2 ${R2} -n ${ID} -o ${OUT_DIR} -p ${N_JOBS}

*VEBA* should register that step 1 is complete and will continue with step 2. 

#### 30. How can I install a developmental module environment on my VEBA?

Developmental/experimental environments can be installed separately.  They are not installed automatically because they use far more compute resources and time than the other environments.  These include the following: 

* `VEBA-amplicon_env`
* `VEBA-biosynthetic_env`

The below code shows how to install the `VEBA-biosynthetic_env` environment as an example:

1. Specify the path to the VEBA repository directory: 

```
VEBA_REPOSITORY_DIRECTORY=path/to/veba_repository_directory (e.g., a release or from git clone https://github.com/jolespin/veba)
```

2. Create the environment: 

```
conda env create -n VEBA-biosynthetic_env -f veba/install/environments/devel/biosynthetic_env.yml
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

For more information, please refer the the [patch update documentation](https://github.com/jolespin/veba/blob/main/install/PATCHES.md#patches).

#### 31. Why are my `unbinned.fasta` files missing or empty?

If you are running `binning-prokaryotic.py` -> `binning-eukaryotic.py` -> `binning-viral.py` or similar you may encounter either a missing `unbinned.fasta` file or an empty `unbinned.fasta` file.  These go hand in hand and are most likely because no contigs met the minimum length requirements leading to an empty file.  If there is an empty `unbinned.fasta` file then trying to use it downstream (e.g., `binning-prokaryotic.py` -> `binning-eukaryotic.py`) will throw an error during file validation and will not produce any `unbinned.fasta` files downstream of this.  You can check this by using the following command: 

```bash
M=1500 # Whatever your minimum contig length is set for when you initially ran
ID= # The ID in question

# Check to see if any sequences are long enough
cat veba_output/assembly/${ID}/output/scaffolds.fasta | seqkit seq -m ${M} | grep -c "^>"

# If the output is 0 then the explanation above holds.
```
  