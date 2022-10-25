### Complete end-to-end metagenomics analysis
This walkthrough goes through assembling metagenomic reads, binning, clustering, classification, and annotation.  We also show how to use the unbinned contigs in a pseudo-coassembly with guidelines on when it's a good idea to go this route. 

What you'll end up with at the end of this are genomes, proteins, gene models, clusters, gene annotations, and *a lot* more.

*VEBA* is modular so feel free to skip a domain or two if you see fit.

_____________________________________________________

#### Steps:

1. Preprocess reads and get directory set up
2. Assemble reads, map reads to assembly, calculate assembly statistics, and index the assembly
3. Recover prokaryotes from metagenomic assemblies
4. Recover eukaryotes from metagenomic assemblies
5. Recover viruses from metagenomic assemblies
6. Concatenating unbinned contigs from assemblies into a pseudo-assembly [Optional]
7. Recover prokaryotes from pseudo-coassembly [Optional]
8. Recover eukaryotes from pseudo-coassembly [Optional]
9. Cluster genomes and proteins
10. Classify prokaryotic genomes
11. Classify eukaryotic genomes
12. Classify viral genomes
13. Annotate proteins

_____________________________________________________

#### 1. Preprocess reads and get directory set up

**Conda Environment:** `conda activate VEBA-preprocess_env`

Refer to the [downloading and preprocessing reads walkthrough](download_and_preprocess_reads.md).  At this point, it's assumed you have the following: 

* A file with all of your identifiers on a separate line (e.g., `identifiers.list` but you can call it whatever you want)
* A directory to keep all your logs called `logs/`
* A directory of preprocessed reads: `veba_output/preprocessed/${ID}/output/cleaned_1.fastq.gz` and `veba_output/preprocessed/${ID}/output/cleaned_2.fastq.gz` where `${ID}` represents the identifiers in `identifiers.list`.

#### 2. Assemble reads, map reads to assembly, calculate assembly statistics, and index the assembly

Here we are going to assemble all of the reads using `metaSPAdes`.  If you have metatranscriptomics, then check out the [walkthrough for recovering viruses from metatranscriptomics](recovering_viruses_from_metatranscriptomics.md).

**Conda Environment:** `conda activate VEBA-assembly_env`

```
# Set the number of threads to use for each sample. Let's use 4
N_JOBS=4

# This is the default output directory
OUT_DIR=veba_output/assembly

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
	CMD="source activate VEBA-assembly_env && assembly.py -1 ${R1} -2 ${R2} -n ${ID} -o ${OUT_DIR} -p ${N_JOBS} --spades_program metaspades.py"
	
	# Either run this command or use SunGridEnginge/SLURM
	
	done

```

The following output files will produced for each sample: 

* featurecounts.tsv.gz - featureCounts output for contig-level counts
* mapped.sorted.bam - Sorted BAM
* mapped.sorted.bam.bai - Sorted BAM index
* scaffolds.fasta - Assembly scaffolds (preferred over contigs by SPAdes documentation)
* scaffolds.fasta.\*.bt2 - Bowtie2 index of scaffolds
* scaffolds.fasta.saf - SAF formatted file for contig-level counts with featureCounts
* seqkit_stats.tsv.gz - Assembly statistics

The main ones we need are `scaffolds.fasta` and `mapped.sorted.bam` which we will use for binning.

*metaSPAdes* contig identifiers have the following format: 

`NODE_1_length_187734_cov_10.705359`

Node, length, and coverage refer to the backend node in the assembly graph, scaffold length, and scaffold coverage, respectively.


#### 3. Recover prokaryotes from metagenomic assemblies
Here we are going to perform iterative prokaryotic binning.  It's difficult to say how many iterations to use but I've found that between 3-10 usually does the job. This dataset is fairly complex and has decent depth so let's use 10 iterations.  Remember, we are analyzing each sample individually at this stage.  

**⚠️Notes:** 

1) If you have a lot of samples and a lot of contigs then use the `--skip_maxbin2` flag because it takes MUCH longer to run.  For the *Plastisphere* it was going to take 40 hours per `MaxBin2` run (there are 2 `MaxBin2` runs) per iteration.  `Metabat2` and `CONCOCT` can do the heavy lifting much faster and often with better results so it's recommended to skip `MaxBin2` for larger datasets.

2) The penultimate step `step-number]__cpr_adjustment]` is the most memory intensive stage as this uses `GTDB-Tk` to identify CPR bacteria. Approaches to properly allocate resources are explained in FAQ. 

Please refer to [FAQ](https://github.com/jolespin/veba/blob/main/FAQ.md) for more details.  In particular, refer to items: 15, 18, 20, and 21.

**Conda Environment:** `conda activate VEBA-binning-prokaryotic_env`

```
N_JOBS=4
N_ITER=10

# This is the default output
OUT_DIR=veba_output/binning/prokaryotic/

for ID in $(cat identifiers.list); do 
	N="binning-prokaryotic__${ID}";
	rm -f logs/${N}.*
	
	# Get the assembly and bam files
	FASTA=veba_output/assembly/${ID}/output/scaffolds.fasta
	BAM=veba_output/assembly/${ID}/output/mapped.sorted.bam
	
	# If your symlinks got deleted, then use these:
	#FASTA=veba_output/assembly/${ID}/intermediate/1__assembly/scaffolds.fasta
	#BAM=veba_output/assembly/${ID}/intermediate/2__alignment/mapped.sorted.bam
	CMD="source activate VEBA-binning-prokaryotic_env && binning-prokaryotic.py -f ${FASTA} -b ${BAM} -n ${ID} -p ${N_JOBS} -o ${OUT_DIR} -m 1500 -I ${N_ITER}"

	# Either run this command or use SunGridEnginge/SLURM

	done
```

The following output files will produced for each sample: 

* binned.list - List of binned contigs
* bins.list - List of MAG identifiers
* checkm_output.filtered.tsv - Filtered CheckM output
* featurecounts.orfs.tsv.gz - ORF-level counts table
* genome_statistics.tsv - Genome assembly statistics
* genomes/ - MAG subdirectory
* genomes/\*.fa - MAG assembly fasta
* genomes/\*.faa - MAG protein fasta
* genomes/\*.ffn - MAG CDS fasta
* genomes/\*.gff - MAG gene models
* genomes/identifier\_mapping.tsv - Identifier mapping between [id\_orf, id\_contig, id\_mag]
* gtdbtk\_output.filtered.tsv - Filtered GTDBTk output
* scaffolds_to_bins.tsv - Identifier mapping between [id_contig, id_mag]
* unbinned.fasta - Fasta of unbinned contigs that have passed length thresholding
* unbinned.list - List of unbinned contigs


#### 4. Recover eukaryotes from metagenomic assemblies
Let's take the unbinned contigs from the prokaryotic binning and recover eukayoritc genomes.  Unfortunately, we aren't going to do iterative binning here because there aren't any tools that can handle consensus genome binning as there is with prokaryotes (e.g., *DAS Tool*).  We have the option to use either *Metabat2* or *CONCOCT*.  In our experience, *Metabat2* works better for recovering eukaryotic genomes from metagenomes and it's also faster as well.

**Conda Environment:** `conda activate VEBA-binning-eukaryotic_env`

```
N_JOBS=4

# This is the default output
OUT_DIR=veba_output/binning/eukaryotic/

for ID in $(cat identifiers.list); do
	N="binning-eukaryotic__${ID}";
	rm -f logs/${N}.*
	FASTA=veba_output/binning/prokaryotic/${ID}/output/unbinned.fasta
	BAM=veba_output/assembly/${ID}/output/mapped.sorted.bam
	CMD="source activate VEBA-binning-eukaryotic_env && binning-eukaryotic.py -f ${FASTA} -b ${BAM} -n ${ID} -p ${N_JOBS} -m 1500 -a metabat2 -o ${OUT_DIR}"
	
	# Either run this command or use SunGridEnginge/SLURM

	done

```
The following output files will produced for each sample: 

* binned.list - List of binned contigs
* bins.list - List of MAG identifiers
* busco_results.filtered.tsv - Filtered BUSCO output
* featurecounts.orfs.tsv.gz - ORF-level counts table
* genome_statistics.tsv - Genome assembly statistics
* genomes/ - MAG subdirectory
* genomes/\*.fa - MAG assembly fasta
* genomes/\*.faa - MAG protein fasta
* genomes/\*.ffn - MAG CDS fasta
* genomes/\*.gff - MAG gene models
* genomes/identifier_mapping.tsv - Identifier mapping between [id_orf, id_contig, id_mag]
* identifier_mapping.metaeuk.tsv - Identifier mapping between original MetaEuk identifiers and modified identifiers.  Includes fully parsed MetaEuk identifiers.
* scaffolds_to_bins.tsv - Identifier mapping between [id_contig, id_mag]
* unbinned.fasta - Fasta of unbinned contigs that have passed length thresholding
* unbinned.list - List of unbinned contigs

#### 5. Recover viruses from metagenomic assemblies
Let's take the unbinned contigs from the eukaryotic binning and recover viral genomes.  Iterative binning doesn't make sense here because *VirFinder* and *CheckV* work on the contig level so we are only pulling out viruses that are on a single contig.  The criteria for high-quality viral genomes are described by the [*CheckV* author](https://scholar.google.com/citations?user=gmKnjNQAAAAJ&hl=en) [here in this Bitbucket Issue (#38)](https://bitbucket.org/berkeleylab/checkv/issues/38/recommended-cutoffs-for-analyzing-checkv).

**Conda Environment:** `conda activate VEBA-binning-viral_env`

```
N_JOBS=4

for ID in $(cat identifiers.list);
	do N="binning-viral__${ID}";
	rm -f logs/${N}.*
	FASTA=veba_output/binning/eukaryotic/${ID}/output/unbinned.fasta
	BAM=veba_output/assembly/${ID}/output/mapped.sorted.bam
	CMD="source activate VEBA-binning-viral_env && binning-viral.py -f ${FASTA} -b ${BAM} -n ${ID} -p ${N_JOBS} -m 1500 -o veba_output/binning/viral"
	# Either run this command or use SunGridEnginge/SLURM

	done
```


The following output files will produced for each sample: 

* binned.list - List of binned contigs
* bins.list - List of MAG identifiers
* quality_summary.filtered.tsv - Filtered CheckV output
* featurecounts.orfs.tsv.gz - ORF-level counts table (If --bam file is provided)
* genome_statistics.tsv - Genome assembly statistics
* genomes/ - MAG subdirectory
* genomes/\*.fa - MAG assembly fasta
* genomes/\*.faa - MAG protein fasta
* genomes/\*.ffn - MAG CDS fasta
* genomes/\*.gff - MAG gene models
* genomes/identifier_mapping.tsv - Identifier mapping between [id_orf, id_contig, id_mag]
* scaffolds_to_bins.tsv - Identifier mapping between [id_contig, id_mag]
* unbinned.fasta - Fasta of unbinned contigs that have passed length thresholding
* unbinned.list - List of unbinned contigs

#### ⚠️ Steps 6-8 are "pseudo-coassembly binning" and require a bit of "common-sense" logic.  Read below to make sure your data qualifies:
This part is optional and should be used with some domain knowledge.  If the dataset you're working with has marine samples from the Pacific and terrestrial lakes then you probably don't want to do this.  However, if you have a bunch of human microbiome samples all within the same cohort from the same body site then you might want to consider running these steps to pull out some more MAGs.  If I wasn't clear, only use this if the biosamples you're analyzing are similar to each other and are assumed to be from a similar population. 

Another bit to consider is specs on your unbinned samples.  If the longest contig is 15k, most are around your lower limit ~1500, and your file size is 10MB it's probably not worth running this.  Similarly, you have a few longer contigs and a higher N50 it could be worth running prokaryotic and not eukaryotic.  Regardless, the binning modules are pretty quick to run when there's low information assemblies so it's up to you whether you feel like it's worth running this.

In the case of the Plastisphere dataset, we found a lot more prokaryotic MAGs using the pseudo-coassembly method.  We also recovered some fragmented eukaryotic MAGs but they didn't pass our QC (i.e., completeness ≥ 50, contamination < 10) so we removed them (by we I mean the program did based on our thresholds).

If you're unsure about how similar you're biosamples are to each other then just skip this to be safe.

#### ⚠️ 6. Concatenating unbinned contigs from assemblies into a pseudo-assembly [Optional]

That said, if you decide to move forward with the multi-sample approach then the next move we want to make is to use the unbinned contigs as they might have fragmented genomes that could be combined with other fragmented genomes from related samples to recover a (nearly) complete genome. We can take the unbinned contigs from the viral binning but we prefer to be more explicit at this step. Once we get the unbinned contigs (i.e., pseudo-coassembly), let's map the reads to them using `coverage.py`.

**Conda Environment:** `conda activate VEBA-assembly_env`

```
# Let's create a sub-directory to store some these intermediate files
mkdir -p veba_output/misc

# Get all the identifiers from all of the MAGs

cat veba_output/binning/*/*/output/genomes/*.fa | grep "^>" | cut -c2- > veba_output/misc/all_sample_specific_mags.binned_contigs.list

# Most VEBA environments should have SeqKit installed.  
# I recommend having this light-weight program in base environment 
# if you do a lot of fasta manipulation. 

conda activate VEBA-preprocess_env

# Now concatenate all the assemblies, pipe the merged assembly into SeqKit, filter out contigs less than the minimum contig length (e.g., 1500 in this case), and grep out the contigs we want using the list file we created in the command above (all_sample_specific_mags.binned_contigs.list) using 4 threads.  Note the `-v` flag in `seqkit grep` as this "invert[s] the sense of matching, to select non-matching records"

N_JOBS=4

cat veba_output/assembly/*/output/scaffolds.fasta | seqkit seq -m 1500 -j ${N_JOBS} | seqkit grep -v -j ${N_JOBS} -f veba_output/misc/all_sample_specific_mags.binned_contigs.list > veba_output/misc/all_sample_specific_mags.unbinned_contigs.gt1500.fasta

# Now create a read table of the following format: 
# [id_sample]<tab>[path/to/r1.fastq.gz]<tab>[path/to/r2.fastq.gz]
# compile_reads_table.py is a script installed in every VEBA environment. 
# If for some reason it's unavailable then just download it here: 
# https://github.com/jolespin/veba/tree/main/src/scripts
# If you want, you can use `-r` or `--relative` flag for relative paths. 
# Note, functionality prior to 2022.10.24 uses `-a` and `--absolute` but now absolute paths are default.

compile_reads_table.py -i veba_output/preprocess/ > veba_output/misc/reads_table.tsv

# Now let's map all the reads to the pseudo-coassembly (i.e., all_sample_specific_mags.unbinned_contigs.gt1500.fasta)

N=pseudo-coassembly

N_JOBS=16 # Let's use more threads here because we are going to be handling multiple samples at once

CMD="source activate VEBA-assembly_env && coverage.py -f veba_output/misc/all_sample_specific_mags.unbinned_contigs.gt1500.fasta -r reads_table.tsv -p ${N_JOBS} -o veba_output/assembly/pseudo-coassembly -m 1500"

# Either run this command or use SunGridEnginge/SLURM
```

The following output files will produced: 

* featurecounts.tsv.gz - featureCounts counts table of all samples
* [sample_id]/mapped.sorted.bam - Sorted BAM file under a subdirectory for each sample
* reference.fasta - Reference fasta (typically this would be the pseudo-coassembly of unbinned contigs)
* reference.fasta.saf - SAF formatted file for contig-level counts with featureCounts
* seqkit_stats.tsv - Assembly statistics

#### ⚠️ 7. Recover prokaryotes from pseudo-coassembly [Optional]
Let's try to recover some prokaryotes using the concatenated unbinned contigs.  

*Note: If you have a lot of samples and a lot of contigs then use the --skip_maxbin2 flag because it takes MUCH longer to run.  For the Plastisphere it was going to take 40 hours per MaxBin2 run (there are 2 MaxBin2 runs) per iteration.  Metabat2 and CONCOCT can do the heavy lifting much faster and often with better results so it's recommended to skip MaxBin2 for larger datasets.*

**Conda Environment:** `conda activate VEBA-binning-prokaryotic_env`

```
# Setting more threads since we are only running this once
N_JOBS=32

# Set the number of iterations.  It's recommended to use the same amount you used for the sample-specific binning
N_ITER=5

# Set up filepaths and names
NAME="pseudo-coassembly"
N="binning-prokaryotic__${NAME}"
rm -f logs/${N}.*
FASTA=veba_output/assembly/pseudo-coassembly/output/reference.fasta
BAMS=veba_output/assembly/pseudo-coassembly/output/*/mapped.sorted.bam

# Set up command
CMD="source activate VEBA-binning-prokaryotic_env && binning-prokaryotic.py -f ${FASTA} -b ${BAMS} -n ${NAME} -p ${N_JOBS} -m 1500 -I ${N_ITER} --skip_maxbin2"

# Either run this command or use SunGridEnginge/SLURM

```

Check Step 3 for the output file descriptions. 


#### ⚠️ 8. Recover eukaryotes from pseudo-coassembly [Optional]
Let's try to recover some eukaryotes using the updated concatenated unbinned contigs. 

**Conda Environment:** `conda activate VEBA-binning-eukaryotic_env`


```
# Setting more threads since we are only running this once
N_JOBS=32


# Set up filepaths and names
NAME="pseudo-coassembly"
N="binning-eukaryotic__${NAME}"
rm -f logs/${N}.*
FASTA=veba_output/binning/prokaryotic/${NAME}/output/unbinned.fasta
BAMS=veba_output/assembly/${NAME}/output/*/mapped.sorted.bam

# Set up command
CMD="source activate VEBA-binning-eukaryotic_env && binning-eukaryotic.py -f ${FASTA} -b ${BAMS} -n ${NAME} -p ${N_JOBS} -m 1500 -a metabat2  -o veba_output/binning/eukaryotic"

# Either run this command or use SunGridEnginge/SLURM

```

Check Step 4 for the output file descriptions. 


#### 9. Cluster genomes and proteins
To analyze these data, we are going to generate some counts tables and we want a single set of features to compare across all samples.  To achieve this, we are going to cluster the genomes into species-level clusters (SLC) and the proteins into SLC-specific orthogroups (SSO).  Further, this clustering is dual purpose as it alleviates some of the bias from [the curse(s) of dimensionality](https://www.nature.com/articles/s41592-018-0019-x) with dimensionality reduction via feature compression - [a type of feature engineering](https://towardsdatascience.com/what-is-feature-engineering-importance-tools-and-techniques-for-machine-learning-2080b0269f10).

**Conda Environment:** `conda activate VEBA-cluster_env`


```
# We need to get genome and protein lists for each domain
# (below assumes you recovered eukaryotic.  If you are missing
# one of the domains below then just remove for for-loop. Note, these paths # will be relative just like the `reads_table.tsv` example in step 6 so as # long as you run from project directory (i.e., the directory in which 
# `veba_output` is a sub-directory) then you're good.

for DOMAIN in prokaryotic eukaryotic viral; do
	ls veba_output/binning/${DOMAIN}/*/output/genomes/*.fa > veba_output/misc/${DOMAIN}__genomes.list
	ls veba_output/binning/${DOMAIN}/*/output/genomes/*.faa > veba_output/misc/${DOMAIN}__proteins.list
	cat veba_output/binning/${DOMAIN}/*/output/scaffolds_to_bins.tsv > veba_output/misc/${DOMAIN}__scaffolds_to_bins.tsv
	done
	
# Let's set an intermediate amount of threads for these
N_JOBS=12

# Now lets run clustering for all the domains separately.
for DOMAIN in prokaryotic eukaryotic viral; do
	N=cluster-${DOMAIN}
	rm -f logs/${N}.*
	
	# This gets the first character from each domain that we are going to use for the clustering
	# (e.g., prokaryotic clustering will be PSLC, eukaryotic as ESLC)
	PREFIX_CHAR=$(echo $DOMAIN | cut -c1)
	
	# Set up the command
	# Note that ^^ just uppercases the prefix
	CMD="source activate VEBA-cluster_env && cluster.py -i veba_output/misc/${DOMAIN}__scaffolds_to_bins.tsv -m veba_output/misc/${DOMAIN}__genomes.list -a veba_output/misc/${DOMAIN}__proteins.list -o veba_output/cluster/${DOMAIN} -p ${N_JOBS} --cluster_prefix ${PREFIX_CHAR^^}SLC"
	
	# Either run this command or use SunGridEnginge/SLURM

	done
```

The following output files will produced for each domain: 


* clusters.tsv - Identifier mapping between MAGs and clusters [id_mag, id_cluster]
* identifier_mapping.orthogroups.tsv - Identifier mapping between ORFs, contigs, orthogroups, and clusters
* proteins_to_orthogroups.tsv - Identifier mapping between ORFs and orthogroups [id_orf, id_orthogroup]
* scaffolds_to_clusters.tsv - Identifier mapping between contigs and clusters [id_contig, id_cluster]

#### 10. Classify prokaryotic genomes
*VEBA* already ran *GTDB-Tk* for you so this module will be quick as it just needs to parse the data and use the clustering information we just computed.  However, you can run *GTDB-Tk* through this as well.

**Conda Environment:** `conda activate VEBA-classify_env`

```
# This isn't threaded so just use 1 or don't even specify
N_JOBS=1

N=classify-prokaryotic
rm -f logs/${N}.*

# Specify the binning directory
BINNING_DIRECTORY=veba_output/binning/prokaryotic

# Clusters are optional but if you have already clustered then use it
CLUSTERS=veba_output/cluster/prokaryotic/output/clusters.tsv

# Set up the command
CMD="source activate VEBA-classify_env && classify-prokaryotic.py -i ${BINNING_DIRECTORY} -c ${CLUSTERS} -p ${N_JOBS} -o veba_output/classify/prokaryotic"

# Either run this command or use SunGridEnginge/SLURM

```

The following output files will produced: 

* prokaryotic_taxonomy.tsv - Prokaryotic genome classification based on GTDBTk
* prokaryotic_taxonomy.clusters.tsv - Prokaryotic cluster classification (If --clusters are provided)

#### 11. Classify eukaryotic genomes
*VEBA* is going to use the *MetaEuk/MMSEQS2* protein alignments based on [*VEBA's* microeukaryotic protein database](https://doi.org/10.6084/m9.figshare.19668855.v1).  The default is to use [BUSCO's eukaryota_odb10](https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2020-09-10.tar.gz) marker set but you can use the annotations from all proteins if you want by providing the `--include_all_genes` flag. The former will take a little bit longer since it needs to run *hmmsearch* but it's more robust and doesn't take that much longer.

**Conda Environment:** `conda activate VEBA-classify_env`

```
# This is threaded if you use the default (i.e., core marker detection)
N_JOBS=4

N=classify-eukaryotic
rm -rf logs/${N}.*

# Specify the binning directory
BINNING_DIRECTORY=veba_output/binning/eukaryotic

# Clusters are optional but if you have already clustered then use it
CLUSTERS=veba_output/cluster/eukaryotic/output/clusters.tsv

# Set up the command
CMD="source activate VEBA-classify_env && classify-eukaryotic.py -i ${BINNING_DIRECTORY} -c ${CLUSTERS} -o veba_output/classify/eukaryotic -p ${N_JOBS}"

# Either run this command or use SunGridEnginge/SLURM

```

The following output files will produced: 

* eukaryotic_taxonomy.tsv - Eukaryotic genome classification based on microeukaryotic protein database and BUSCO's eukaryota_odb10 marker set
* eukaryotic_taxonomy.clusters.tsv - Eukaryotic cluster classification (If --clusters are provided)
* gene-source_lineage.tsv - Gene source lineage and scores for classifying MAGs [id_gene, id_scaffold, id_mag, id_target, id_source, lineage, bitscore]

#### 12. Classify viral genomes
Viral classification  is a bit unchartered territory so the classifications are quite low level.  We use the VOG classification from *CheckV* because it is more common.  We also use the "isolation source" information when clusters are provided.

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

The following output files will produced: 

* viral_taxonomy.tsv - Viral genome classification based on CheckV classifications
* viral_taxonomy.clusters.tsv - Viral cluster classification (If --clusters are provided)
* viral_isolation_source.clusters.tsv - Viral isolation source consensus (If --clusters are provided)

#### 13. Annotate proteins
Now that allf of the MAGs are recovered and classified, let's annotate the proteins using best-hit against NR, Pfam, and KOFAM.

**Conda Environment:** `conda activate VEBA-annotate_env`

```
# Let's merge all of the proteins.  

cat veba_output/binning/*/*/output/genomes/*.faa > veba_output/misc/all_genomes.proteins.faa

# Let's merge all of the identifier mappings. 
# Note, this is optional but it's helpful for diving further into annotations
# especially for viruses
cat veba_output/binning/*/*/output/genomes/identifier_mapping.tsv > veba_output/misc/all_genomes.identifier_mapping.tsv

# If you are restricted by resources or time you may want to do this in batches.  
# You can split up the proteins into separate batches with SeqKit.
# This would make 100 files: stdin.part_001.fasta - stdin.part_100.fasta

N_PARTITIONS=100
mkdir -p veba_output/misc/partitions/
veba_output/misc/all_genomes.proteins.faa | seqkit split2 -p ${N_PARTITIONS} -O veba_output/misc/partitions/

# However, we aren't doing that for sake of simplicity.  Let's just annotate everything at once:

# Let's use a higher number of threads here
N_JOBS=32

# Set up log names and output paths
N="annotate"
rm -f logs/${N}.*
OUT_DIR="veba_output/annotation/"
mkdir -p ${OUT_DIR}

# Input filepaths
PROTEINS=veba_output/misc/all_genomes.proteins.faa
IDENTIFIER_MAPPING=veba_output/misc/all_genomes.identifier_mapping.tsv
CMD="source activate VEBA-annotate_env && annotate.py -a ${PROTEINS} -i ${IDENTIFIER_MAPPING} -o ${OUT_DIR} -p ${N_JOBS}"

# Either run this command or use SunGridEnginge/SLURM

```

The following output files will produced: 

* annotations.orfs.tsv.gz - Concatenated annotations from Diamond (NR), HMMSearch (Pfam), and KOFAMSCAN (KOFAM)
* lineage_predictions.contigs.tsv.gz - [Experimental] Lineage predictions for contigs based on NR annotations.  This should be used only for experimentation as many lineages aren't defined properly. Contig-level classifications.
* lineage_predictions.mags.tsv.gz - [Experimental] Lineage predictions for contigs based on NR annotations.  This should be used only for experimentation as many lineages aren't defined properly. MAG-level classifications.

_____________________________________________________

#### Next steps:

Now it's time to map the reads so we can analyze the data using compositional data analysis (CoDA).  Please refer to the [read mapping and counts tables walkthrough](read_mapping_and_counts_tables.md). 

