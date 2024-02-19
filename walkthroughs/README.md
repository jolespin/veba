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

##### Accessing SRA: 

*  **[Downloading and preprocessing fastq files](docs/download_and_preprocess_reads.md)** - Explains how to download reads from NCBI and run *VEBA's* `preprocess.py` module to decontaminate either metagenomic and/or metatranscriptomic reads.

##### End-to-end workflows:

* **[Complete end-to-end metagenomics analysis](docs/end-to-end_metagenomics.md)** - Goes through assembling metagenomic reads, binning, clustering, classification, and annotation.  We also show how to use the unbinned contigs in a pseudo-coassembly with guidelines on when it's a good idea to go this route.
* **[Recovering viruses from metatranscriptomics](docs/recovering_viruses_from_metatranscriptomics.md)** - Goes through assembling metatranscriptomic reads, viral binning, clustering, and classification.
* **[Setting up *bona fide* coassemblies for metagenomics or metatranscriptomics](docs/setting_up_coassemblies.md)** - In the case where all samples are of low depth, it may be useful to use coassembly instead of sample-specific approaches.  This walkthrough goes through concatenating reads, creating a reads table, coassembly of concatenated reads, aligning sample-specific reads to the coassembly for multiple sorted BAM files, and mapping reads for scaffold/transcript-level counts.  Please note that a coassembly differs from the pseudo-coassembly concept introduced in the VEBA publication.  For more information regarding the differences between *bona fide* coassembly and pseud-coassembly, please refer to [*23. What's the difference between a coassembly and a pseudo-coassembly?*](https://github.com/jolespin/veba/blob/main/FAQ.md#23-whats-the-difference-between-a-coassembly-and-a-pseudo-coassembly). 

##### Phylogenetics:

* **[Phylogenetic inference](docs/phylogenetic_inference.md)** - Phylogenetic inference of eukaryotic diatoms.

##### Bioprospecting:

* **[Bioprospecting for biosynthetic gene clusters](docs/bioprospecting_for_biosynthetic_gene_clusters.md)** - Detecting biosynthetic gene clusters (BGC) with and scoring novelty of BGCs.
* **[CRISPR-Cas system screening with *de novo* genomes](docs/crispr-cas_system_screening_de-novo_genomes.md)** How use `CRISPRCasTyper` as a post hoc analysis for screening genomes.

##### Read mapping, rapid profiling, feature engineering, and converting counts tables:

* **[Taxonomic profiling *de novo* genomes](docs/taxonomic_profiling_de-novo_genomes.md)** - Explains how to build and profile reads to custom `Sylph` databases from *de novo* genomes.
* **[Pathway profiling *de novo* genomes](docs/pathway_profiling_de-novo_genomes.md)** - Explains how to build and align reads to custom `HUMAnN` databases from *de novo* genomes and annotations.
* **[Read mapping and counts tables](docs/read_mapping_and_counts_tables.md)** - Traditional read mapping and generating counts tables at the contig, MAG, SLC, ORF, and SSO levels. 
* **[Merging counts tables with taxonomy](docs/merging_counts_with_taxonomy.md)** - Explains how to merge counts tables with taxonomy.
* **[Phylogenomic functional categories using *de novo* genomes](docs/phylogenomic_functional_categories.md)** - PhyloGenomic Functional Categories (PGFC) using annotations, clusters, and counts tables as implemented in [*Espinoza et al. 2022*](https://academic.oup.com/pnasnexus/article/1/5/pgac239/6762943).
* **[Converting counts tables](docs/converting_counts_tables.md)** - Convert your counts table (with or without metadata) to [anndata](https://anndata.readthedocs.io/en/latest/index.html) or [biom](https://biom-format.org/) format.  Also supports [Pandas pickle](docs/https://pandas.pydata.org/docs/reference/api/pandas.read_pickle.html) format.

##### Containerization and AWS:

* **[Adapting commands for Docker](docs/adapting_commands_for_docker.md)** - Explains how to download and use Docker for running VEBA.
* **[Adapting commands for AWS](docs/adapting_commands_for_aws.md)** - Explains how to download and use Docker for running VEBA specifically on AWS.

___________________________________________

**Coming Soon:**

* Visualizing genome-clusters with `NetworkX`
* Workflow for low-depth samples with no bins
* Assigning eukaryotic taxonomy to unbinned contigs (`metaeuk taxtocontig`)
* Bioprospecting using [`PlasticDB` database](https://plasticdb.org/)
* Targeted pathway profiling of large and complex reference databases
___________________________________________

##### Notes:

* The final output files are in the `output` subdirectory and to avoid redundant files many of these symlinked from the `intermediate` directory.  This can cause issues if you are using a ["scratch" directory](https://en.wikipedia.org/wiki/Scratch_space) where the files are deleted after a certain amount of time. If you have a [crontab](https://www.man7.org/linux/man-pages/man5/crontab.5.html) set up, make sure it also touches symlinks and not just files.
* You'll need to adjust the memory and time for different jobs.  Assembly will take much longer than preprocessing.  Annotation will require more memory than mapping. 