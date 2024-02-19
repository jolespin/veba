### Merging counts tables with taxonomy classifications
This walkthrough goes through merging  your counts table with your taxonomy classifications.

What you'll end up with at the end of this is either a `.tsv.gz` with your counts tables merged to your taxonomy classifications.

It is assumed you've completed either the [end-to-end metagenomics](end-to-end_metagenomics.md) or the [recovering viruses from metatranscriptomics](recovering_viruses_from_metatranscriptomics.md) and either the [read mapping/counts tables](read_mapping_and_counts_tables.md) or [taxonomic profiling](docs/taxonomic_profiling_de-novo_genomes.md) walkthroughs.


_____________________________________________________

#### Steps:

1. Merge taxonomy from different organism types
2. Merge counts tables
3. Merge counts table with taxonomy classifications

**Conda Environment:** `conda activate VEBA`. Use this for intermediate scripts.

_____________________________________________________


#### 1. Merge taxonomy from different organism types

```bash
merge_taxonomy_classifications.py -i veba_output/classify/ -o veba_output/classify/
```

This command outputs the following tables: 

* veba_output/classify/taxonomy_classifications.tsv - Table with [id_genome]<tab>[id_organism_type]<tab>[classification]
* veba_output/classify/taxonomy_classifications.clusters.tsv - Table with [id_genome]<tab>[id_organism_type]<tab>[classification]<tab>[more fields with clustering information]


#### 2. Merge counts table

Now let's merge the counts tables.  For counts tables from the `mapping` module, if you haven't done so already you can merge them with this script: 

```bash
# Get mapping directory
MAPPING_DIRECTORY=veba_output/mapping/global

# Set output directory (this is default)
OUT_DIR=veba_output/counts
SCAFFOLDS_TO_MAGS=veba_output/cluster/output/global/scaffolds_to_mags.tsv
SCAFFOLDS_TO_SLCS=veba_output/cluster/output/global/scaffolds_to_slcs.tsv

merge_contig_mapping.py -m ${MAPPING_DIRECTORY} -c ${MAGS_TO_SLCS}  -i ${SCAFFOLDS_TO_MAGS} -o ${OUT_DIR}
```

You'll get the following output files: 

* X_mags.tsv.gz - Counts tables gzipped and tab-delimited (samples, MAGs)
* X_slcs.tsv.gz - Counts tables gzipped and tab-delimited (samples, SLCs)

If you have taxonomic profiling abundances instead: 

```bash
merge_generalized_mapping.py -o veba_output/profiling/taxonomy/merged.taxonomic_abundance.tsv.gz veba_output/profiling/taxonomy/*/output/taxonomic_abundance.tsv.gz

merge_generalized_mapping.py -o veba_output/profiling/taxonomy/merged.taxonomic_abundance.clusters.tsv.gz veba_output/profiling/taxonomy/*/output/taxonomic_abundance.clusters.tsv.gz
```

You'll get the following output files: 

* merged.taxonomic\_abundance.tsv.gz - Merged genome-level taxonomic abundance matrix
* merged.taxonomic\_abundance.clusters.tsv.gz - Merged SLC-level taxonomic abundance matrix

 
#### 3. Merge counts table with taxonomy classifications


Now that the taxonomy is merged from different organism types and counts tables are merged, we can join them. 

If you have `mapping` module output from above or similar you can use this: 

```bash

merge_counts_with_taxonomy.py -X veba_output/counts/X_mags.tsv.gz -t veba_output/classify/taxonomy_classifications.tsv -o veba_output/counts/X_mags.with_taxonomy.tsv.gz

merge_counts_with_taxonomy.py -X veba_output/counts/X_slcs.tsv.gz -t veba_output/classify/taxonomy_classifications.clusters.tsv -o veba_output/counts/X_slcs.with_taxonomy.tsv.gz
```

You'll get the following output which will contain transposed counts (rows=features) with organism type and classifications prepended:

* veba_output/counts/X_mags.with_taxonomy.tsv.gz 
* veba_output/counts/X_slcs.with_taxonomy.tsv.gz

If you have `profile-taxonomy` module output from above or similar you can use this: 

```bash
merge_counts_with_taxonomy.py -X veba_output/profiling/taxonomy/merged.taxonomic_abundance.tsv.gz -t veba_output/classify/taxonomy_classifications.tsv -o veba_output/profiling/taxonomy/merged.taxonomic_abundance.with_taxonomy.tsv.gz

merge_counts_with_taxonomy.py -X veba_output/profiling/taxonomy/merged.taxonomic_abundance.clusters.tsv.gz -t veba_output/classify/taxonomy_classifications.clusters.tsv -o veba_output/profiling/taxonomy/merged.taxonomic_abundance.clusters.with_taxonomy.tsv.gz
```

You'll get the following output which will contain transposed counts (rows=features) with organism type and classifications prepended:

* merged.taxonomic_abundance.with_taxonomy.tsv.gz
* merged.taxonomic_abundance.clusters.with_taxonomy.tsv.gz


#### Next steps:

Now it's time to analyze the data (recommend using compositional data analysis (CoDA) approach).  

If you have a `pickle` format, then you probably know what to do already.

If you have a `biom` format look into the [QIIME2 ecosystem](https://docs.qiime2.org/2023.2/tutorials/) or [BIRDMAn](https://birdman.readthedocs.io/en/stable/?badge=stable) for CoDA-aware bayesian differential abundance. 

If you have an `anndata` format, look into the [Scanpy ecosystem](https://scanpy.readthedocs.io/en/stable/tutorials.html).  Though Scanpy does not nativly support CoDA-aware methodologies, they are easy to adapt.

##### Recommended reading for CoDA:

* Thomas P Quinn, Ionas Erb, Mark F Richardson, Tamsyn M Crowley, Understanding sequencing data as compositions: an outlook and review, Bioinformatics, Volume 34, Issue 16, 15 August 2018, Pages 2870–2878, https://doi.org/10.1093/bioinformatics/bty175
* Thomas P Quinn, Ionas Erb, Greg Gloor, Cedric Notredame, Mark F Richardson, Tamsyn M Crowley, A field guide for the compositional analysis of any-omics data, GigaScience, Volume 8, Issue 9, September 2019, giz107, https://doi.org/10.1093/gigascience/giz107
* Espinoza, J.L., Shah, N., Singh, S., Nelson, K.E. and Dupont, C.L. (2020), Applications of weighted association networks applied to compositional data in biology. Environ Microbiol, 22: 3020-3038. https://doi.org/10.1111/1462-2920.15091