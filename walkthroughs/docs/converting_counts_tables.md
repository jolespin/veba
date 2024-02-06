### Converting counts table to biom or anndata format
This walkthrough goes through converting your counts table (with or without metadata) to [anndata](https://anndata.readthedocs.io/en/latest/index.html) or [biom](https://biom-format.org/) format.

What you'll end up with at the end of this is either a `.h5ad` or `.biom` format file for `anndata` and `biom`, respecitvely. 

It is assumed you've completed either the [end-to-end metagenomics](end-to-end_metagenomics.md) or the [recovering viruses from metatranscriptomics](recovering_viruses_from_metatranscriptomics.md) and the [read mapping/counts tables](read_mapping_and_counts_tables.md) walkthroughs.

`anndata` is developed for single-cell genomics but I've been using it for my metagenomics workflow.  I recommend trying it out.  

_____________________________________________________

#### Steps:

1. Provide just the counts table
2. Provide a counts table and sample metadata
3. Provide a counts table, sample metadata, and 

**Conda Environment:** `conda activate VEBA`. Use this for intermediate scripts.

_____________________________________________________


#### 1. Let's convert to a Python pickle object without any metadata

```bash
convert_counts_table.py -i veba_output/counts/X_slcs.tsv.gz -f pickle -o veba_output/counts/X_slcs.no_metadata.pkl
```

Now in Python we can load it with `Pandas`

```python
import pandas as pd
X_slcs = pd.read_pickle("veba_output/counts/X_slcs.no_metadata.pkl")
```


#### 2. Let's convert to a biom object with sample metadata

You probably have sample metadata so use that instead.  Make sure that all the sample identifiers in your counts table are accounted for in your metadata table.  

Here we don't have any so I'm going to just get the reads manifest for a toy sample metadata table. 

```bash
# Generate a sample metadata for example but use your own instead
# Fastq/
# Fastq/S1_1.fastq.gz
# Fastq/S1_2.fastq.gz
# Fastq/etc.

compile_reads_table.py -f Fastq/ -n [ID]_[DIRECTION].fastq.gz --header > veba_output/misc/reads_table.tsv
```

Convert to biom format using the sample metadata

```bash

# Convert table
convert_counts_table.py -i veba_output/counts/X_slcs.tsv.gz -f biom -o veba_output/counts/X_slcs.with_sample_metadata.biom --sample_metadata veba_output/misc/reads_table.tsv

```

Now in Python we can load it with `biom`

```python
import biom
table = biom.load_table("veba_output/counts/X_slcs.with_sample_metadata.biom")
# 16 x 4 <class 'biom.table.Table'> with 59 nonzero entries (92% dense)
```
 
#### 3. Let's convert to a anndata object with sample and feature metadata

We are going to use the sample metadata we generated in previous step. We are going to use the taxonomy classifications from the `classify.py` module for the feature metadata.

```bash
convert_counts_table.py -i veba_output/counts/X_slcs.tsv.gz -f anndata -o veba_output/counts/X_slcs.with_sample_and_feature.metadata.h5ad --sample_metadata veba_output/misc/reads_table.tsv --feature_metadata veba_output/classify/taxonomy_classifications.clusters.tsv
```

Now in Python we can load it with `anndata`

```python
import anndata as ad
adata = ad.read_h5ad("veba_output/counts/X_slcs.with_sample_and_feature.metadata.h5ad")
#AnnData object with n_obs × n_vars = 4 × 16
#    obs: 'forward-absolute-filepath', 'reverse-absolute-filepath'
#    var: 'domain', 'consensus_classification', 'homogeneity', 'number_of_unique_classifications', #'number_of_genomes', 'genomes', 'classifications', 'weights', 'score', #'number_of_unique_classification', 'consensus_taxon_id'
```


#### Next steps:

Now it's time to analyze the data (recommend using compositional data analysis (CoDA) approach).  

If you have a `pickle` format, then you probably know what to do already.

If you have a `biom` format look into the [QIIME2 ecosystem](https://docs.qiime2.org/2023.2/tutorials/) or [BIRDMAn](https://birdman.readthedocs.io/en/stable/?badge=stable) for CoDA-aware bayesian differential abundance. 

If you have an `anndata` format, look into the [Scanpy ecosystem](https://scanpy.readthedocs.io/en/stable/tutorials.html).  Though Scanpy does not nativly support CoDA-aware methodologies, they are easy to adapt.

##### Recommended reading for CoDA:

* Thomas P Quinn, Ionas Erb, Mark F Richardson, Tamsyn M Crowley, Understanding sequencing data as compositions: an outlook and review, Bioinformatics, Volume 34, Issue 16, 15 August 2018, Pages 2870–2878, https://doi.org/10.1093/bioinformatics/bty175
* Thomas P Quinn, Ionas Erb, Greg Gloor, Cedric Notredame, Mark F Richardson, Tamsyn M Crowley, A field guide for the compositional analysis of any-omics data, GigaScience, Volume 8, Issue 9, September 2019, giz107, https://doi.org/10.1093/gigascience/giz107
* Espinoza, J.L., Shah, N., Singh, S., Nelson, K.E. and Dupont, C.L. (2020), Applications of weighted association networks applied to compositional data in biology. Environ Microbiol, 22: 3020-3038. https://doi.org/10.1111/1462-2920.15091