### Pathway profiling of *de novo* genomes
If you build a comprehensive database, you may want to use a read-based approach to functionally profile a large set of samples.  This tutorial will show you how to build a custom HUMAnN database from your annotations and how to profile your samples where there is full accounting of reads and your genomes.

What you'll end up with at the end of this is a merged taxonomy table, a custom HUMAnN annotation table, and HUMAnN profiles.

Please refer to the [end-to-end metagenomics](end-to-end_metagenomics.md) or [recovering viruses from metatranscriptomics](recovering_viruses_from_metatranscriptomics.md) workflows for details on binning, clustering, and annotation.

_____________________________________________________

#### Steps:

1. Merge taxonomy from all domains
2. Compile custom `HUMAnN` database from *de novo* genomes
3. Functional profiling using `HUMAnN` of custom database
4. Merge the tables

**Conda Environment:** `conda activate VEBA`. Use this for intermediate scripts.
_______________________________________________________

#### 1.  Merge taxonomy from all domains

At this point, it's assumed you have the following: 

* A concatenated protein fasta (preferably with proteins ≥ 100k removed `seqkit seq -M 99999`) called `veba_output/misc/all_genomes.all_proteins.lt100k.faa`
* Taxonomy classifications for all domains
* Annotations results from the `annotate.py` module
* Clustering results from the `cluster.py` module
* A directory of preprocessed reads: `veba_output/preprocess/${ID}/output/cleaned_1.fastq.gz` and `veba_output/preprocess/${ID}/output/cleaned_2.fastq.gz` where `${ID}` represents the identifiers in `identifiers.list`.

```bash
merge_taxonomy_classifications.py -i veba_output/classify/ -o veba_output/classify/ --no_header --no_domain
```

#### 2. Compile custom `HUMAnN` database from *de novo* genomes

Here we are going to compile all the `UniRef` annotations and taxonomy identifiers to build a custom `HUMAnN` database.  The script supports piping to stdin so we are going to cut some columns from `identifier_mapping.proteins.tsv.gz`


```
zcat veba_output/cluster/output/global/identifier_mapping.proteins.tsv.gz  | cut -f1,3 | tail -n +2 | compile_custom_humann_database_from_annotations.py -a veba_output/annotation/output/annotations.proteins.tsv.gz -s veba_output/misc/all_genomes.all_proteins.lt100k.faa -o veba_output/misc/humann_uniref_annotations.tsv -t veba_output/classify/taxonomy_classifications.tsv

```

This generates a table that has the following columns (no header):

* `[id_protein]<tab>[id_uniref]<tab>[length]<tab>[lineage]`
* e.g., `S1__NODE_1648_length_2909_cov_4.311843_4        UniRef50_A0A1C7D7V0     214     d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Erythrobacter;s__Erythrobacter sp003149575;t__S1__CONCOCT__P.1__24`

#### 3. Functional profiling using `HUMAnN` of custom database

Now it's time to align reads.  Since `HUMAnN` takes in single reads, we join the reads using `bbmerge.sh` in the backend.  You could also use the aligned reads in the form of a `bam` file from the `mapping.py` module.

```
N_JOBS=4
OUT_DIR=veba_output/profiling/pathways
UNIREF_ANNOTATIONS=veba_output/misc/humann_uniref_annotations.tsv
FASTA=veba_output/misc/all_genomes.all_proteins.lt100k.faa

mkdir -p logs

for ID in $(cat identifiers.list);
do 
	N="profile-pathway__${ID}";
	rm -f logs/${N}.*
	R1=veba_output/preprocess/${ID}/output/cleaned_1.fastq.gz
	R2=veba_output/preprocess/${ID}/output/cleaned_2.fastq.gz
	CMD="source activate VEBA && veba --module profile-pathway --params \"-1 ${R1} -2 ${R2} -n ${ID} -o ${OUT_DIR} -p ${N_JOBS} -i ${UNIREF_ANNOTATIONS} -f ${FASTA}\""
	
	# Either run this command or use SunGridEnginge/SLURM

done

```


The following output files will produced for each sample: 

* reads.seqkit_stats.tsv - Sequence stats for input reads
* humann\_pathabundance.tsv - Stratified abundance of taxa-specific metabolic pathways
* humann\_pathcoverage.tsv - Stratified pathway completion ratio of taxa-specific metabolic pathways
* humann\_genefamilies.tsv - Stratified abundance of taxa-specific gene families
* humann\_diamond\_unaligned.fa.gz - Joined reads that did not align to database
* humann\_diamond\_aligned.tsv.gz - Aligned reads from translated blast search to database (blast6 format)


#### 4. Merge the tables

```
merge_generalized_mapping.py -o veba_output/profiling/pathways/merged.humann_pathcoverage.tsv veba_output/profiling/pathways/*/output/humann_pathcoverage.tsv

merge_generalized_mapping.py -o veba_output/profiling/pathways/merged.humann_pathabundance.tsv veba_output/profiling/pathways/*/output/humann_pathabundance.tsv
 
merge_generalized_mapping.py -o veba_output/profiling/pathways/merged.humann_genefamilies.tsv veba_output/profiling/pathways/*/output/humann_genefamilies.tsv
```

The following output files will produced for each sample: 

* merged.humann\_pathabundance.tsv - Stratified abundance of taxa-specific metabolic pathways of all samples
* meged.humann\_pathcoverage.tsv - Stratified pathway completion ratio of taxa-specific metabolic of all samplespathways
* meged.humann\_genefamilies.tsv - Stratified abundance of taxa-specific gene families of all samples


_____________________________________________________

#### Next steps:

Subset stratified tables by their respective levels.
