# Getting started with *VEBA*

#### Basics:
Since VEBA functionality benefits from structure, it's good to have a list of identifiers that you can use for for-loops. In the examples, it will be the following: `identifiers.list`

However, for datasets with metagenomics and metatranscriptomics it's often useful to have a master list `identifiers.list` and separate `identifiers.dna.list` and `identifiers.rna.list` for metagenomic and metatranscriptomic samples, respectively. 

Our VEBA project directory is going to be `veba_output` and each step will be in a subdirectory.  

* e.g., `veba_output/preprocess` will have all of the preprocessed reads.  

In the workflows that work on specific samples, there will be sample subdirectories. 

* e.g., `veba_output/preprocess/SRR17458603/output`, `veba_output/preprocess/SRR17458606/output `, ...
* e.g., `veba_output/assembly/SRR17458603/output `, `veba_output/assembly/SRR17458606/output `, ...

Many of these jobs should be run using a job scheduler like [SunGridEngine](https://docs.oracle.com/cd/E19279-01/820-3257-12/n1ge.html) or [SLURM](https://slurm.schedmd.com/documentation.html).  This [resource](https://www.miamioh.edu/research/research-computing-support/services/hpc-cluster/sbatch-translation/) is useful for converting commands between SunGridEnginer and SLURM. I've used both and these are adaptations of the submission commands you can use as a template:

```
# Let's create some informative name. Remember we are going to create a lot of jobs and log files for the different workflows if you have multiple samples
N=preprocessing__${ID}
	
CMD="some command we want to run"
	
# SunGridEngine:
qsub -o logs/${N}.o -e logs/${N}.e -cwd -N ${N} -j y -pe threaded ${N_JOBS} "${CMD}"
	
# SLURM:
sbatch -J ${N} -N 1 -c ${N_JOBS} --ntasks-per-node=1 -o logs/${N}.o -e logs/${N}.e --export=ALL -t 12:00:00 --mem=20G --wrap="${CMD}"
```

#### Available walkthroughs:

*  [Downloading and preprocessing fastq files](download_and_preprocess_reads.md) - Explains how to download reads from NCBI and run *VEBA's* `preprocess.py` module to decontaminate either metagenomic and/or metatranscriptomic reads.
* [Complete end-to-end metagenomics analysis](end-to-end_metagenomics.md) - Goes through assembling metagenomic reads, binning, clustering, classification, and annotation.  We also show how to use the unbinned contigs in a pseudo-coassembly with guidelines on when it's a good idea to go this route.
*  [Recovering viruses from metatranscriptomics](recovering_viruses_from_metatranscriptomics.md) - Goes through assembling metatranscriptomic reads, viral binning, clustering, and classification.
*  [Read mapping and counts tables](read_mapping_and_counts_tables.md) - Read mapping and generating counts tables at the contig, MAG, SLC, ORF, and SSO levels. 
* [Phylogenetic inference](phylogenetic_inference.md) - Phylogenetic inference of eukaryotic diatoms.
* [Setting up *bona fide* coassemblies for metagenomics or metatranscriptomics](setting_up_coassemblies.md) - In the case where all samples are of low depth, it may be useful to use coassembly instead of sample-specific approaches.  This walkthrough goes through concatenating reads, creating a reads table, coassembly of concatenated reads, aligning sample-specific reads to the coassembly for multiple sorted BAM files, and mapping reads for scaffold/transcript-level counts.  Please note that a coassembly differs from the pseudo-coassembly concept introduced in the VEBA publication.  For more information regarding the differences between *bona fide* coassembly and pseud-coassembly, please refer to [*23. What's the difference between a coassembly and a pseudo-coassembly?*](https://github.com/jolespin/veba/blob/main/FAQ.md). 
* [Bioprospecting for biosynthetic gene clusters](bioprospecting_for_biosynthetic_gene_clusters.md) - Detecting biosynthetic gene clusters (BGC) with and scoring novelty of BGCs.

___________________________________________

**Coming Soon:**

* Workflow for low-depth samples with no bins
* Workflow for ASV detection from short-read amplicons
* Workflows for integrating 3rd party software with *VEBA*:
	* Using [EukHeist](https://github.com/AlexanderLabWHOI/EukHeist) for eukaryotic binning followed by *VEBA* for mapping and annotation.
	* Using [EukMetaSanity](https://github.com/cjneely10/EukMetaSanity) for modeling genes for eukaryotic genomes recovered with *VEBA*.

___________________________________________

##### Notes:

* The final output files are in the `output` subdirectory and to avoid redundant files many of these symlinked from the `intermediate` directory.  This can cause issues if you are using a ["scratch" directory](https://en.wikipedia.org/wiki/Scratch_space) where the files are deleted after a certain amount of time. If you have a [crontab](https://www.man7.org/linux/man-pages/man5/crontab.5.html) set up, make sure it also touches symlinks and not just files.
* You'll need to adjust the memory and time for different jobs.  Assembly will take much longer than preprocessing.  Annotation will require more memory than mapping. 