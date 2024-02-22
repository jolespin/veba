### Taxonomic profiling of *de novo* genomes
If you build a comprehensive database, you may want to use a read-based approach to taxonomicly profile a large set of samples.  This tutorial will show you how to build a custom `Sylph` database from your genomes and how to profile your samples for taxonomic abundance.

What you'll end up with at the end of this is a `Sylph` database and taxonomic abundance profiles.

Please refer to the [end-to-end metagenomics](end-to-end_metagenomics.md) or [recovering viruses from metatranscriptomics](recovering_viruses_from_metatranscriptomics.md) workflows for details on binning, clustering, and annotation.

_____________________________________________________

#### Steps:
1. Compile custom `Sylph` database from *de novo* genomes
2. Taxonomic profiling using `Sylph ` of custom database
3. Merge the tables

**Conda Environment:** `conda activate VEBA`. Use this for intermediate scripts.
_______________________________________________________

#### 1. Compile custom `Sylph` database from *de novo* genomes

At this point, it's assumed you have the following: 

* Clustering results from the `cluster.py` module
* A directory of preprocessed reads: `veba_output/preprocess/${ID}/output/cleaned_1.fastq.gz` and `veba_output/preprocess/${ID}/output/cleaned_2.fastq.gz` where `${ID}` represents the identifiers in `identifiers.list`.
* Genome assemblies.  These can either be MAGs binned with VEBA, binned elsewhere, or even reference genomes you downloaded. 


Here we are going to build 2 databases, one for viral genomes and one for non-viral genomes (i.e., prokaryotes and eukaryotes).  The reason for 2 separate databases is because there are presets used for small genomes that are different than medium to large genomes.  We need a table that has `[organism_type]<tab>[path/to/genome.fa]` with no headers.  We already have some version of this with the `veba_output/misc/genomes_table.tsv` we made for clustering.  We can pipe this into stdin for the database build script:


```
cat veba_output/misc/genomes_table.tsv | cut -f1,4 | compile_custom_sylph_sketch_database_from_genomes.py -o veba_output/profiling/databases
```

This generates 2 `Sylph` databases (assuming you have viruses and non-viruses):

* `veba_output/profiling/databases/genome_database-nonviral.syldb`
* `veba_output/profiling/databases/genome_database-viral.syldb`


#### 2. Taxonomic profiling using `Sylph` of custom database

Now it's time to profile the reads against the `Sylph` databases.  Since `Sylph` takes in a sketch of reads, we can either use a precompute reads sketch with `-s` or with paired-end reads (`-1` and `-2`) to compute the sketch in the backend.  

```
N_JOBS=4
OUT_DIR=veba_output/profiling/taxonomy
DATABASES=veba_output/profiling/databases/*.syldb
MAGS_TO_SLCS=veba_output/cluster/output/global/mags_to_slcs.tsv # Assuming you have clustering results

mkdir -p logs

for ID in $(cat identifiers.list);
do 
	N="profile-taxonomy__${ID}";
	rm -f logs/${N}.*
	R1=veba_output/preprocess/${ID}/output/cleaned_1.fastq.gz
	R2=veba_output/preprocess/${ID}/output/cleaned_2.fastq.gz
	CMD="source activate VEBA && veba --module profile-taxonomy --params \"-1 ${R1} -2 ${R2} -n ${ID} -o ${OUT_DIR} -p ${N_JOBS} -d ${DATABASES} -c ${MAGS_TO_SLCS}\""
	
	# Either run this command or use SunGridEnginge/SLURM

done

```

The following output files will produced for each sample: 

* reads.sylsp - Reads sketch if paired-end reads were provided
* sylph\_profile.tsv.gz - Output of `sylph profile`
* taxonomic_abundance.tsv.gz - Genome-level taxonomic abundance (No header)
* taxonomic_abundance.clusters.tsv.gz - SLC-level taxonomic abundance (No header, if --genome_clusters wer provided)

#### 3. Merge the tables

```
merge_generalized_mapping.py -o veba_output/profiling/taxonomy/merged.taxonomic_abundance.tsv.gz veba_output/profiling/taxonomy/*/output/taxonomic_abundance.tsv.gz

merge_generalized_mapping.py -o veba_output/profiling/taxonomy/merged.taxonomic_abundance.clusters.tsv.gz veba_output/profiling/taxonomy/*/output/taxonomic_abundance.clusters.tsv.gz
```

The following output files will produced for each sample: 

* merged.taxonomic\_abundance.tsv.gz - Merged genome-level taxonomic abundance matrix
* merged.taxonomic\_abundance.clusters.tsv.gz - Merged SLC-level taxonomic abundance matrix

_____________________________________________________

#### Next steps:

Subset stratified tables by their respective levels.
