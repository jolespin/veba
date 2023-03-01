### Recovering viruses from metatranscriptomics
This walkthrough goes through assembling metatranscriptomic reads, viral binning, clustering, and classification.

What you'll end up with at the end of this are genomes, proteins, gene models, clusters, and *a lot* more.

Please refer to the [end-to-end metagenomics workflow](end-to-end_metagenomics.md) for details on annotating viral genes.

_____________________________________________________

#### Steps:

1. Preprocess reads and get directory set up
2. Assemble reads, map reads to assembly, calculate assembly statistics, and index the assembly
3. Recover viruses from metatranscriptomic assemblies
4. Cluster genomes and proteins
5. Classify viral genomes

#### 1. Preprocess reads and get directory set up

**Conda Environment:** `conda activate VEBA-preprocess_env`

Refer to the [downloading and preprocessing reads workflow](download_and_preprocess_reads.md).  At this point, it's assumed you have the following: 

* A file with all of your identifiers on a separate line (e.g., `identifiers.list` but you can call it whatever you want)
* A directory to keep all your logs called `logs/`
* A directory of preprocessed reads: `veba_output/preprocess/${ID}/output/cleaned_1.fastq.gz` and `veba_output/preprocess/${ID}/output/cleaned_2.fastq.gz` where `${ID}` represents the identifiers in `identifiers.list`.

#### 2. Assemble reads, map reads to assembly, calculate assembly statistics, and index the assembly

Here we are going to assemble all of the reads using `rnaSPAdes`.  

**Conda Environment:** `conda activate VEBA-assembly_env`

```
# Set the number of threads to use for each sample. Let's use 4
N_JOBS=4

# This is not the default output directory and is customized for metatranscriptomics
OUT_DIR=veba_output/transcript_assembly

# Just making sure your log directory is available
mkdir -p logs 

# Iterate through the identifiers and run the assembly command
for ID in $(cat identifiers.list); do
	
	# Set a useful name for log files
	N="assembly__${ID}";
	rm -f logs/${N}.*
	
	# Get forward and reverse reads
	R1=veba_output/preprocess/${ID}/output/cleaned_1.fastq.gz
	R2=veba_output/preprocess/${ID}/output/cleaned_2.fastq.gz
	
	# Set up command
	CMD="source activate VEBA-assembly_env && assembly.py -1 ${R1} -2 ${R2} -n ${ID} -o ${OUT_DIR} -p ${N_JOBS} -P rnaspades.py"
	
	# Either run this command or use SunGridEnginge/SLURM
	
	done

```

**The following output files will produced for each sample:** 

* featurecounts.tsv.gz - featureCounts output for contig-level counts
* mapped.sorted.bam - Sorted BAM
* mapped.sorted.bam.bai - Sorted BAM index
* genes\_to\_transcripts.tsv - Gene identifier to transcript identifier
* transcripts.fasta - Transcript assemblies
* transcripts.fasta.\*.bt2 - Bowtie2 index of transcripts
* transcripts.fasta.saf - SAF formatted file for transcript-level counts with featureCounts
* seqkit_stats.tsv.gz - Assembly statistics

The main one we need is `transcripts.fasta`  which we will use for binning.

*rnaSPAdes* transcript identifiers have the following format: 

`NODE_1_length_117535_cov_30.642667_g0_i0`

Where `g0` refers to the predicted gene and `i0` refers to the isoform transcript for `g0`. Node, length, and coverage refer to the backend node in the assembly graph, transcript length, and transcript coverage, respectively.  This can be useful for tools like [*TransDecoder*](https://github.com/TransDecoder/TransDecoder) that take a `--gene_trans_map` argument (TransDecoder.LongOrfs) and a `--single_best_only` (TransDecoder.Predict) as we implemented when investigating [coral holobiomes](https://www.science.org/doi/10.1126/sciadv.abg3088).  Although,  *VEBA* has [a wrapper script for *TransDecoder* and automated homology searches](https://github.com/jolespin/veba/blob/main/src/scripts/transdecoder_wrapper.py), this tutorial will use *Prodigal* instead.

#### 3. Recover viruses from metatranscriptomic assemblies
We use a similar approach to the metagenomics with *geNomad* and *CheckV* but using the assembled transcripts instead.  Again, the criteria for high-quality viral genomes are described by the [*CheckV* author](https://scholar.google.com/citations?user=gmKnjNQAAAAJ&hl=en) [here in this Bitbucket Issue (#38)](https://bitbucket.org/berkeleylab/checkv/issues/38/recommended-cutoffs-for-analyzing-checkv).

**Conda Environment:** `conda activate VEBA-binning-viral_env`

```
N_JOBS=4

for ID in $(cat identifiers.list);
	do N="binning-viral__${ID}";
	rm -f logs/${N}.*
	FASTA=veba_output/transcript_assembly/${ID}/output/transcripts.fasta
	BAM=veba_output/transcript_assembly/${ID}/output/mapped.sorted.bam
	CMD="source activate VEBA-binning-viral_env && binning-viral.py -f ${FASTA} -b ${BAM} -n ${ID} -p ${N_JOBS} -m 1500 -o veba_output/binning/viral -a genomad"
	
	# Either run this command or use SunGridEnginge/SLURM

	done
```


**The following output files will produced for each sample:** 

* binned.list - List of binned transcripts
* bins.list - List of MAG identifiers
* quality_summary.filtered.tsv - Filtered CheckV output
* featurecounts.orfs.tsv.gz - ORF-level counts table (If --bam file is provided)
* genome_statistics.tsv - Genome assembly statistics
* genomes/ - MAG subdirectory
* genomes/\*.fa - MAG assembly fasta
* genomes/\*.faa - MAG protein fasta
* genomes/\*.ffn - MAG CDS fasta
* genomes/\*.gff - MAG gene models
* genomes/identifier_mapping.tsv - Identifier mapping between [id_orf, id_transcript, id_mag]
* scaffolds\_to\_bins.tsv - Identifier mapping between [id_transcript, id_mag]
* unbinned.fasta - Fasta of unbinned transcripts that have passed length thresholding
* unbinned.list - List of unbinned transcripts


#### 4. Cluster genomes and proteins
To analyze these data, we are going to generate some counts tables and we want a single set of features to compare across all samples.  To achieve this, we are going to cluster the genomes into species-level clusters (SLC) and the proteins into SLC-specific protein clusters (SSPC).  Further, this clustering is dual purpose as it alleviates some of the bias from [the curse(s) of dimensionality](https://www.nature.com/articles/s41592-018-0019-x) with dimensionality reduction via feature compression - [a type of feature engineering](https://towardsdatascience.com/what-is-feature-engineering-importance-tools-and-techniques-for-machine-learning-2080b0269f10).

**Conda Environment:** `conda activate VEBA-cluster_env`


```
# We need to generate a table with the following fields:
# [organism_type]<tab>[id_sample]<tab>[id_mag]<tab>[genome]<tab>[proteins]
# The table but can include additional columns to the right (e.g., [cds]<tab>[gene_models]) but the first 5 columns are mandatory. 

compile_genomes_table.py -i veba_output/binning/ > veba_output/misc/genomes_table.tsv
	
# Let's set an intermediate amount of threads for these
N_JOBS=12

# Set up command
CMD="source activate VEBA-cluster_env && cluster.py -i veba_output/misc/genomes_table.tsv -o veba_output/cluster -p ${N_JOBS}"
	
# Either run this command or use SunGridEnginge/SLURM
```

**The following output files will produced for each domain:**

**Default genome cluster prefix:** `VSLC-1` where `V` refers to the capitalized first letter of the domain (e.g., `V`iral)

**Default protein cluster prefix:** `VSLC-1_SSPC-1` 

`veba_output/cluster/global`

`veba_output/cluster/local`

* genome\_clusters.tsv - Table with genome cluster, number of components, number of samples of origin, component identifiers, and sample identifiers of origin.
* identifier\_mapping.genomes.tsv - Identifier mapping between genome, organism-type, sample of origin, id\_genome\_cluster, and number of proteins
* identifier\_mapping.proteins.tsv - Identifier mapping between id\_protein, organism-type, id\_genome, sample of origin, id\_genome\_cluster, and id\_protein\_cluster
* identifier\_mapping.scaffolds.tsv - Identifier mapping between id_scaffold, organism-type, id\_genome, sample of origin, and id\_genome\_cluster
* mags\_to\_slcs.tsv - Identifier mapping between MAGs and genome clusters [id\_genome, id\_genome\_cluster]
* feature\_compression\_ratios.tsv - Genomic and functional feature compression ratios (FCR) for each domain
* protein\_clusters.tsv - Table with protein cluster, number of components, number of samples of origin, component identifiers, and sample identifiers of origin.
* proteins\_to\_orthogroups.tsv - Identifier mapping between id\_protein and id\_protein\_cluster
* scaffolds\_to\_slcs.tsv - Identifier mapping between id_scaffold and id\_genome\_cluster



#### 5. Classify viral genomes
Viral classification is performed using `geNomad`.  Classification can be performed using the intermediate binning results which is much quicker.  Alternatively, if you have viruses identified elsewhere you can still classify using the `--genomes` argument instead.

**Conda Environment:** `conda activate VEBA-classify_env`


```
N=classify-viral

# Specify the binning directory
BINNING_DIRECTORY=veba_output/binning/viral

# Clusters are optional but if you have already clustered then use it
CLUSTERS=veba_output/cluster/viral/output/clusters.tsv

rm -rf logs/${N}.*

# Set up the command
CMD="source activate VEBA-classify_env && classify-viral.py -i ${BINNING_DIRECTORY} -c ${CLUSTERS} -o veba_output/classify/viral"

# Either run this command or use SunGridEnginge/SLURM

```

**The following output files will produced:** 

* viral\_taxonomy.tsv - Viral genome classification based on geNomad classifications
* viral\_taxonomy.clusters.tsv - Viral cluster classification (If `--clusters` are provided)


_____________________________________________________

#### Next steps:

Now it's time to map the reads so we can analyze the data using compositional data analysis (CoDA).  Please refer to the [read mapping and counts tables workflow](read_mapping_and_counts_tables.md). 
