### Complete end-to-end metagenomics analysis
This walkthrough goes through assembling metagenomic reads, binning, clustering, classification, and annotation.  We also show how to use the unbinned contigs in a pseudo-coassembly with guidelines on when it's a good idea to go this route. 

What you'll end up with at the end of this are genomes, proteins, gene models, clusters, gene annotations, and *a lot* more.

*VEBA* is modular so feel free to skip a domain or two if you see fit.

_____________________________________________________

#### Steps:

1. Preprocess reads and get directory set up
2. Assemble reads, map reads to assembly, calculate assembly statistics, and index the assembly
3. Recover viruses from metagenomic assemblies
4. Recover prokaryotes from metagenomic assemblies
5. Recover eukaryotes from metagenomic assemblies
6. Concatenating unbinned contigs from assemblies into a pseudo-assembly [Optional]
7. Recover prokaryotes from pseudo-coassembly [Optional]
8. Recover eukaryotes from pseudo-coassembly [Optional]
9. Cluster genomes and proteins
10. Classify viral genomes
11. Classify prokaryotic genomes
12. Classify eukaryotic genomes
13. Annotate proteins

_____________________________________________________

#### 1. Preprocess reads and get directory set up

**Conda Environment:** `conda activate VEBA-preprocess_env`

Refer to the [downloading and preprocessing reads walkthrough](download_and_preprocess_reads.md).  At this point, it's assumed you have the following: 

* A file with all of your identifiers on a separate line (e.g., `identifiers.list` but you can call it whatever you want)
* A directory to keep all your logs called `logs/`
* A directory of preprocessed reads: `veba_output/preprocess/${ID}/output/cleaned_1.fastq.gz` and `veba_output/preprocess/${ID}/output/cleaned_2.fastq.gz` where `${ID}` represents the identifiers in `identifiers.list`.

#### 2. Assemble reads, map reads to assembly, calculate assembly statistics, and index the assembly

Here we are going to assemble all of the reads using `metaSPAdes`.  If you have metatranscriptomics, then check out the [walkthrough for recovering viruses from metatranscriptomics](recovering_viruses_from_metatranscriptomics.md).  

**Recommended memory request:** For this *Plastisphere* dataset, I requested `64GB` of memory from my HPC.  Though, this will change depending on how deep your samples are sequenced.

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
	CMD="source activate VEBA-assembly_env && assembly.py -1 ${R1} -2 ${R2} -n ${ID} -o ${OUT_DIR} -p ${N_JOBS} -P metaspades.py"
	
	# Either run this command or use SunGridEnginge/SLURM
	
	done

```

**The following output files will produced for each sample:** 

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

#### 3. Recover viruses from metagenomic assemblies
In previous tutorials, we performed prokaryotic and eukaryotic binning first but removing viral contigs first should make the prokaryotic and eukaryotic binning more refined since we are removing noise.

Let's start the binning with viruses since this is performed on a per-contig bases (instead of a collection of contigs).  Iterative binning doesn't make sense here because *geNomad*, *VirFinder* and *CheckV* work on the contig level so we are only pulling out viruses that are on a single contig.  The criteria for high-quality viral genomes are described by the [*CheckV* author](https://scholar.google.com/citations?user=gmKnjNQAAAAJ&hl=en) [here in this Bitbucket Issue (#38)](https://bitbucket.org/berkeleylab/checkv/issues/38/recommended-cutoffs-for-analyzing-checkv).  Though, these values are quite stringent and can be lowered to recover more candidates.  

**Recommended memory request:** `16 GB`

**Conda Environment:** `conda activate VEBA-binning-viral_env`

```
N_JOBS=4

for ID in $(cat identifiers.list);
	do N="binning-viral__${ID}";
	rm -f logs/${N}.*
	FASTA=veba_output/assembly/${ID}/output/scaffolds.fasta
	BAM=veba_output/assembly/${ID}/output/mapped.sorted.bam
	CMD="source activate VEBA-binning-viral_env && binning-viral.py -f ${FASTA} -b ${BAM} -n ${ID} -p ${N_JOBS} -m 1500 -o veba_output/binning/viral"
	# Either run this command or use SunGridEnginge/SLURM

	done
```


**The following output files will produced for each sample:** 

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



#### 4. Recover prokaryotes from metagenomic assemblies
Here we are going to perform iterative prokaryotic binning.  It's difficult to say how many iterations to use but I've found that between 3-10 usually does the job. This dataset is fairly complex and has decent depth so let's use 10 iterations.  Remember, we are analyzing each sample individually at this stage.  

If you have a lot of samples and a lot of contigs then use the `--skip_maxbin2` flag because it takes MUCH longer to run.  For the *Plastisphere* it was going to take 40 hours per `MaxBin2` run (there are 2 `MaxBin2` runs) per iteration.  `Metabat2` and `CONCOCT` can do the heavy lifting much faster and often with better results so it's recommended to skip `MaxBin2` for larger datasets.  

**Recommended memory request:** `10GB` 

*Versions prior to `v1.1.0` were reliant on `GTDB-Tk` which needed at least `60GB`.  `GTDB-Tk` is no longer required with the update of `CheckM` to `CheckM2`.*

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
	FASTA=veba_output/binning/viral/${ID}/output/unbinned.fasta
	BAM=veba_output/assembly/${ID}/output/mapped.sorted.bam
	
	CMD="source activate VEBA-binning-prokaryotic_env && binning-prokaryotic.py -f ${FASTA} -b ${BAM} -n ${ID} -p ${N_JOBS} -o ${OUT_DIR} -m 1500 -I ${N_ITER}"

	# Either run this command or use SunGridEnginge/SLURM

	done
```

**The following output files will produced for each sample:** 

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


#### 5. Recover eukaryotes from metagenomic assemblies
Let's take the unbinned contigs from the prokaryotic binning and recover eukayoritc genomes.  Unfortunately, we aren't going to do iterative binning here because there aren't any tools that can handle consensus genome binning as there is with prokaryotes (e.g., *DAS Tool*).  We have the option to use either *Metabat2* or *CONCOCT*.  In our experience, *Metabat2* works better for recovering eukaryotic genomes from metagenomes and it's also faster as well.

**Recommended memory request:** `128GB` 

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
**The following output files will produced for each sample:** 

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

___________________________________________________________________________________

#### ⚠️ Steps 6-8 are "pseudo-coassembly binning" and require a bit of "common-sense" logic.  Read below to make sure your data qualifies:
This part is optional and should be used with some domain knowledge.  If the dataset you're working with has marine samples from the Pacific and terrestrial lakes then you probably don't want to do this.  However, if you have a bunch of human microbiome samples all within the same cohort from the same body site then you might want to consider running these steps to pull out some more MAGs.  If I wasn't clear, only use this if the biosamples you're analyzing are similar to each other and are assumed to be from a similar population. 

Another bit to consider is specs on your unbinned samples.  If the longest contig is 15k, most are around your lower limit ~1500, and your file size is 10MB it's probably not worth running this.  Similarly, you have a few longer contigs and a higher N50 it could be worth running prokaryotic and not eukaryotic.  Regardless, the binning modules are pretty quick to run when there's low information assemblies so it's up to you whether you feel like it's worth running this.

In the case of the Plastisphere dataset, we found a lot more prokaryotic MAGs using the pseudo-coassembly method.  We also recovered some fragmented eukaryotic MAGs but they didn't pass our QC (i.e., completeness ≥ 50, contamination < 10) so we removed them (by we I mean the program did based on our thresholds).

If you're unsure about how similar you're biosamples are to each other then just skip this to be safe.

#### ⚠️ 6. Concatenating unbinned contigs from assemblies into a pseudo-assembly [Optional]

That said, if you decide to move forward with the multi-sample approach then the next move we want to make is to use the unbinned contigs as they might have fragmented genomes that could be combined with other fragmented genomes from related samples to recover a (nearly) complete genome. We can take the unbinned contigs from the viral binning but we prefer to be more explicit at this step. Once we get the unbinned contigs (i.e., pseudo-coassembly), let's map the reads to them using `coverage.py`.

**Recommended memory request:** `24 GB`

**Conda Environment:** `conda activate VEBA-assembly_env`

```


# Let's create a sub-directory to store some these intermediate files
mkdir -p veba_output/misc

# Activate any VEBA environment
# Most VEBA environments should have SeqKit installed.  
# I recommend having this light-weight program in base environment 
# if you do a lot of fasta manipulation. 

conda activate VEBA-preprocess_env

# --------------------------------------------------------------------

# Method 1) Shortcut
cat veba_output/binning/eukaryotic/*/output/unbinned.fasta > veba_output/misc/all_sample_specific_mags.unbinned_contigs.gt1500.fasta

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

# Method 2) Longer more explicit way
# Get all the identifiers from all of the MAGs
cat veba_output/binning/*/*/output/genomes/*.fa | grep "^>" | cut -c2- | cut -f1 -d " " > veba_output/misc/all_sample_specific_mags.binned_contigs.list

# Now concatenate all the assemblies, pipe the merged assembly into SeqKit, filter out contigs less than the minimum contig length (e.g., 1500 in this case), and grep out the contigs we want using the list file we created in the command above (all_sample_specific_mags.binned_contigs.list) using 4 threads.  Note the `-v` flag in `seqkit grep` as this "invert[s] the sense of matching, to select non-matching records"

N_JOBS=4

cat veba_output/assembly/*/output/scaffolds.fasta | seqkit seq -m 1500 -j ${N_JOBS} | seqkit grep -v -j ${N_JOBS} -f veba_output/misc/all_sample_specific_mags.binned_contigs.list > veba_output/misc/all_sample_specific_mags.unbinned_contigs.gt1500.fasta

# --------------------------------------------------------------------

# Now create a read table of the following format: 
# [id_sample]<tab>[path/to/r1.fastq.gz]<tab>[path/to/r2.fastq.gz]
# compile_reads_table.py is a script installed in every VEBA environment. 
# If for some reason it's unavailable then just download it here: 
# https://github.com/jolespin/veba/tree/main/src/scripts
# If you want, you can use `-r` or `--relative` flag for relative paths. 

compile_reads_table.py -i veba_output/preprocess/ -r > veba_output/misc/reads_table.tsv

# Now let's map all the reads to the pseudo-coassembly (i.e., all_sample_specific_mags.unbinned_contigs.gt1500.fasta)

N=pseudo-coassembly

N_JOBS=16 # Let's use more threads here because we are going to be handling multiple samples at once

CMD="source activate VEBA-assembly_env && coverage.py -f veba_output/misc/all_sample_specific_mags.unbinned_contigs.gt1500.fasta -r veba_output/misc/reads_table.tsv -p ${N_JOBS} -o veba_output/assembly/pseudo-coassembly -m 1500"

# Either run this command or use SunGridEnginge/SLURM
```

**The following output files will produced:** 

* featurecounts.tsv.gz - featureCounts counts table of all samples
* [sample_id]/mapped.sorted.bam - Sorted BAM file under a subdirectory for each sample
* reference.fasta - Reference fasta (typically this would be the pseudo-coassembly of unbinned contigs)
* reference.fasta.saf - SAF formatted file for contig-level counts with featureCounts
* seqkit_stats.tsv - Assembly statistics

#### ⚠️ 7. Recover prokaryotes from pseudo-coassembly [Optional]
Let's try to recover some prokaryotes using the concatenated unbinned contigs.  

*Note: If you have a lot of samples and a lot of contigs then use the `--skip_maxbin2` flag because it takes MUCH longer to run.  For the Plastisphere it was going to take 40 hours per MaxBin2 run (there are 2 MaxBin2 runs) per iteration.  Metabat2 and CONCOCT can do the heavy lifting much faster and often with better results so it's recommended to skip MaxBin2 for larger datasets. If you have a lot of samples compared to contigs then CONCOCT may take a lot longer, if you want to skip then you can use the `--skip_concoct` flag.*

**Recommended memory request:** `10 - 24GB` 


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

Check Step 4 for the output file descriptions. 


#### ⚠️ 8. Recover eukaryotes from pseudo-coassembly [Optional]
Let's try to recover some eukaryotes using the updated concatenated unbinned contigs. 

**Recommended memory request:** `128GB` 

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

Check Step 5 for the output file descriptions. 

___________________________________________________________________________________


#### 9. Cluster genomes and proteins
To analyze these data, we are going to generate some counts tables and we want a single set of features to compare across all samples.  To achieve this, we are going to cluster the genomes into species-level clusters (SLC) and the proteins into SLC-specific protein clusters (SSPC).  Further, this clustering is dual purpose as it alleviates some of the bias from [the curse(s) of dimensionality](https://www.nature.com/articles/s41592-018-0019-x) with dimensionality reduction via feature compression - [a type of feature engineering](https://towardsdatascience.com/what-is-feature-engineering-importance-tools-and-techniques-for-machine-learning-2080b0269f10).  In versions prior to `v1.1.0`, clustering needed to be performed for each domain individually.  In `v1.1.0+`, clustering has been optimized with using `MMSEQS2` instead of `OrthoFinder` and is performed all at once (both globally and locally).  This has been made even easier with the usage of `compile_genomes_table.py` script that comes equipped with `VEBA`.

**Recommended memory request:** `24 GB` should work for most datasets but you may need to increase for much larger datasets.

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

**Default genome cluster prefix:** `ESLC-1` where `E` refers to the capitalized first letter of the domain (e.g., `E`ukaryotic)

**Default protein cluster prefix:** `ESLC-1_SSPC-1` 

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


#### 10. Classify viral genomes
Viral classification is performed using `geNomad`.  Classification can be performed using the intermediate binning results which is much quicker.  Alternatively, if you have viruses identified elsewhere you can still classify using the `--genomes` argument instead.

**Recommended memory request:** `1 GB` will work if you've performed viral binning via *VEBA*.  If not, these use `16 GB` for external genomes. 

**Conda Environment:** `conda activate VEBA-classify_env`

```
N=classify-viral

N_JOBS=1

# Specify the binning directory
BINNING_DIRECTORY=veba_output/binning/viral

# Clusters are optional but if you have already clustered then use it.  
# The clusters below contain prokaryotic, eukaryotic, and viral clusters but it just uses the viral clusters
CLUSTERS=veba_output/cluster/output/global/mags_to_slcs.tsv

rm -rf logs/${N}.*

# Set up the command
CMD="source activate VEBA-classify_env && classify-viral.py -i ${BINNING_DIRECTORY} -c ${CLUSTERS} -o veba_output/classify/viral -p ${N_JOBS}"

# Either run this command or use SunGridEnginge/SLURM

```

**The following output files will produced:** 

* taxonomy.tsv - Viral genome classification based on geNomad classifications
* taxonomy.clusters.tsv - Viral cluster classification (If `--clusters` are provided)

#### 11. Classify prokaryotic genomes
Prokaryotic classification is performed using `GTDB-Tk`.  Classification can be performed using the intermediate binning results which is easier.  Alternatively, if you have prokaryotes identified elsewhere you can still classify using the `--genomes` argument instead.

**Recommended memory request:** `72 GB`

**Conda Environment:** `conda activate VEBA-classify_env`

```
N_JOBS=16

N=classify-prokaryotic
rm -f logs/${N}.*

# Specify the binning directory
BINNING_DIRECTORY=veba_output/binning/prokaryotic

# Clusters are optional but if you have already clustered then use it.  
# The clusters below contain prokaryotic, eukaryotic, and viral clusters but it just uses the prokaryotic clusters
CLUSTERS=veba_output/cluster/output/global/mags_to_slcs.tsv

# Set up the command
CMD="source activate VEBA-classify_env && classify-prokaryotic.py -i ${BINNING_DIRECTORY} -c ${CLUSTERS} -p ${N_JOBS} -o veba_output/classify/prokaryotic"

# Either run this command or use SunGridEnginge/SLURM

```

The following output files will produced: 

* taxonomy.tsv - Prokaryotic genome classification based on GTDB-Tk
* taxonomy.clusters.tsv - Prokaryotic cluster classification (If --clusters are provided)

#### 12. Classify eukaryotic genomes
*VEBA* is going to use the *MetaEuk/MMSEQS2* protein alignments based on [*VEBA's* microeukaryotic protein database](https://doi.org/10.6084/m9.figshare.19668855.v1).  The default is to use [BUSCO's eukaryota_odb10](https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2020-09-10.tar.gz) marker set but you can use the annotations from all proteins if you want by providing the `--include_all_genes` flag. The former will take a little bit longer since it needs to run *hmmsearch* but it's more robust and doesn't take that much longer.

**Recommended memory request:** `12 GB`

**Conda Environment:** `conda activate VEBA-classify_env`

```
# This is threaded if you use the default (i.e., core marker detection)
N_JOBS=4

N=classify-eukaryotic
rm -rf logs/${N}.*

# Specify the binning directory
BINNING_DIRECTORY=veba_output/binning/eukaryotic

# Clusters are optional but if you have already clustered then use it
CLUSTERS=veba_output/cluster/output/global/mags_to_slcs.tsv

# Set up the command
CMD="source activate VEBA-classify_env && classify-eukaryotic.py -i ${BINNING_DIRECTORY} -c ${CLUSTERS} -o veba_output/classify/eukaryotic -p ${N_JOBS}"

# Either run this command or use SunGridEnginge/SLURM

```

The following output files will produced: 

* taxonomy.tsv - Eukaryotic genome classification based on microeukaryotic protein database and BUSCO's eukaryota_odb10 marker set
* taxonomy.clusters.tsv - Eukaryotic cluster classification (If --clusters are provided)
* gene-source\_lineage.tsv - Gene source lineage and scores for classifying MAGs [id_gene, id_scaffold, id_mag, id_target, id_source, lineage, bitscore]

#### 13. Merge classifications
Instead of having 3 separate classification tables, it would be much more useful to have 1 single classifcation table for viruses, prokaryotes, and eukaryotes.  To do this, use the following wrapper script: 

**Recommended memory request:** `1 GB`


**Conda Environment:** `conda activate VEBA-classify_env`

```
merge_taxonomy_classifications.py -i veba_output/classify -o veba_output/classify
```

The following output files will produced: 

* taxonomy_classifications.tsv - Taxonomy with respect to genomes 
* taxonomy_classifications.clusters.tsv - Taxonomy with respect to genome clusters

#### 14. Annotate proteins
Now that allf of the MAGs are recovered and classified, let's annotate the proteins using best-hit against UniRef, Pfam, and KOFAM.

**Conda Environment:** `conda activate VEBA-annotate_env`

```
# Let's merge all of the proteins.  

cat veba_output/binning/*/*/output/genomes/*.faa > veba_output/misc/all_genomes.proteins.faa

# Let's merge all of the identifier mappings. 
# Note, this is optional but it's helpful for diving further into annotations
# especially for viruses
cat veba_output/binning/*/*/output/genomes/identifier_mapping.tsv > veba_output/misc/all_genomes.identifier_mapping.tsv

```

**Simple usage by running all at once**

For sake of simplicity, let's just annotate everything at once (see next section to speed this up):

```

# Let's use a higher number of threads here
N_JOBS=32

# Set up log names and output paths
N="annotate"

rm -f logs/${N}.*

# Input filepaths
PROTEINS=veba_output/misc/all_genomes.proteins.faa
IDENTIFIER_MAPPING=veba_output/misc/all_genomes.identifier_mapping.tsv
CMD="source activate VEBA-annotate_env && annotate.py -a ${PROTEINS} -i ${IDENTIFIER_MAPPING} -o veba_output/annotation -p ${N_JOBS}"

# Either run this command or use SunGridEnginge/SLURM

```

**Advanced usage by splitting up fasta into multiple files**

If you are restricted by resources or time you may want to do this in batches.  
You can split up the proteins into separate batches with `SeqKit`.
This would make 100 files: stdin.part_001.fasta - stdin.part_100.fasta

```
# Split fasta file into multiple partitions
N_PARTITIONS=100
PARTITION_DIRECTORY=veba_output/misc/partitions/
mkdir -p ${PARTITION_DIRECTORY}
cat veba_output/misc/all_genomes.proteins.faa | seqkit split2 -p ${N_PARTITIONS} -O ${PARTITION_DIRECTORY}

# Now get the identifiers for each partition
IDENTIFIER_MAPPING=veba_output/misc/all_genomes.identifier_mapping.tsv

for i in $(seq -f "%03g" 1 ${N_PARTITIONS}); do # This iterates through 1-100 and zero pads
	PROTEINS=${PARTITION_DIRECTORY}/stdin.part_${i}.fasta
	IDENTIFIERS=${PARTITION_DIRECTORY}/stdin.part_${i}.list
	cat ${PROTEINS} | grep "^>" | cut -c2- | cut -f1 -d " " > ${IDENTIFIERS}
	subset_table.py -i ${IDENTIFIERS} -t ${IDENTIFIER_MAPPING} -o ${PARTITION_DIRECTORY}/stdin.part_${i}.identifier_mapping.tsv
	done

# Drop down the number of threads used since we are running 100 jobs now
N_JOBS=1

# Set up log names and output paths
OUT_DIR="veba_output/annotations_partitioned"
mkdir -p ${OUT_DIR}

# Now run annotate for each partition
for i in $(seq -f "%03g" 1 ${N_PARTITIONS}); do 
	N="annotate-${i}"
	rm -f logs/${N}.*
	FAA=${PARTITION_DIRECTORY}/stdin.part_${i}.fasta
	IDS=${PARTITION_DIRECTORY}/stdin.part_${i}.identifier_mapping.tsv
	CMD="source activate VEBA-annotate_env && annotate.py -a ${FAA} -i ${IDS} -o ${OUT_DIR}/${i} -p ${N_JOBS}"

	# Either run this command or use SunGridEnginge/SLURM
	
	done

```

The following output files will produced: 

* annotations.tsv.gz - Concatenated annotations from Diamond (NR), HMMSearch (Pfam), HMMSearch (NCBIfam-AMRFinder), HMMSearch (AntiFam), and KOFAMSCAN (KOFAM)
* lineage.weighted_majority_vote.contigs.tsv.gz - [Experimental] Lineage predictions for contigs based on NR annotations.  This should be used only for experimentation as lineages are not determined using core markers. Contig-level classifications.
* lineage.weighted_majority_vote.genomes.tsv.gz - [Experimental] Lineage predictions for contigs based on NR annotations.  This should be used only for experimentation as lineages are not determined using core markers. MAG-level classifications.

_____________________________________________________

#### Next steps:

Now it's time to map the reads so we can analyze the data using compositional data analysis (CoDA).  Please refer to the [read mapping and counts tables walkthrough](read_mapping_and_counts_tables.md). 

