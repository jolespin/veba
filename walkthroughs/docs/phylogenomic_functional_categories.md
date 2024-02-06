### Phylogenomic functional categories using *de novo* genomes
It's common for HUMAnN-based analysis to result in features without any resolved taxonomy.  If this is an issue for downstream analysis and you already have ORF/SSPC-level counts tables then we provide an alternative.  This tutorial will show you how to engineer your features into amalgamations known as PhyoGenomic Functional Catergories (PGFC) as implemented in [*Espinoza et al. 2022*](https://academic.oup.com/pnasnexus/article/1/5/pgac239/6762943).

What you'll end up with at the end of this is PGFC counts table and MCR ratio table for each PGFC.

Please refer to the [end-to-end metagenomics](end-to-end_metagenomics.md) or [recovering viruses from metatranscriptomics](recovering_viruses_from_metatranscriptomics.md) workflows for details on binning, clustering, and annotation.

_____________________________________________________

#### Steps:
0. Install `ensemble_networkx` (if not available)
1. Compile PGFCs
2. Loading `CategoricalEngineeredFeature` object into Python [Optional]

**Conda Environment:** `conda activate VEBA-mapping_env`
_______________________________________________________

#### 0. Install `ensemble_networkx` (if not available)

```
pip install "ensemble_networkx>=2024.2.5"
```

#### 1. Compile PGFCs

At this point, it's assumed you have the following: 

* Clustering results from the `cluster.py` module
* Annotations at the SSPC level from the `annotate.py` module
* Counts tables at the SSPC level from the `mapping.py` module (or `bowtie2|star_wrapper.py`)


Here we are going to use the results from previous workflows to compile the PGFCs which can be performed using the `compile_phylogenomic_functional_categories.py`.  This is not multi-threaded and the amount of memory will be dependent on your dataset but typical datasets will require between 32GB-64GB of memory. 


```
NAME="TestDataset"
ANNOTATION_RESULTS=veba_output/annotation/output/annotations.protein_clusters.tsv.gz
COUNTS=veba_output/counts/X_orthogroups.tsv.gz
FASTA=veba_output/cluster/output/global/representative_sequences.faa
OUTPUT_DIRECTORY=veba_output/phylogenomic_functional_categories/

# Compile PGFCs
compile_phylogenomic_functional_categories.py -i ${ANNOTATION_RESULTS} -X ${COUNTS} -f ${FASTA} -o ${OUTPUT_DIRECTORY} --level genome_cluster -n ${NAME}
```

This generates the following output files:

* categorical_engineered_feature.pgfc.genome_clusters.pkl - Serialized pickle file of `CategoricalEngineeredFeature` object from `ensemble_networkx`  
* categorical_engineered_feature.pgfc.genome_clusters.synopsis.tsv.gz - Metadata and summary statistics for PGFCs
* categorical_engineered_feature.pgfc.genome_clusters.transformed_counts.tsv.gz - Transformed counts tables for PGFCs (Rows=Samples, Columns=PGFCs)
* categorical_engineered_feature.pgfc.genome_clusters.module_completion_ratios.tsv.gz - Module completion ratios for PGFCS (Rows=Samples, Columns=PGFCs)


The main output files look like this: 

* categorical_engineered_feature.pgfc.genome_clusters.transformed_counts.tsv.gz

```
python -c 'import sys, pandas as pd; pd.read_csv("veba_output/phylogenomic_functional_categories/categorical_engineered_feature.pgfc.genome_clusters.transformed_counts.tsv.gz", sep="\t", index_col=0, header=[0,1]).iloc[:4,:6].to_csv(sys.stdout, sep="\t")'

Taxonomy	ESLC-1	ESLC-1	ESLC-1	ESLC-1	ESLC-1	ESLC-1
Functional	M00007	M00165	M00167	M00004	M00580	M00005
id_sample
S4	5	7	7	15	0	1
S3	0	0	0	0	0	0
S1	278	847	724	1086	38	146
S2	337	982	851	1243	30	150
```

* categorical_engineered_feature.pgfc.genome_clusters.module_completion_ratios.tsv.gz

```
 python -c 'import sys, pandas as pd; pd.read_csv("veba_output/phylogenomic_functional_categories/categorical_engineered_feature.pgfc.genome_clusters.module_completion_ratios.tsv.gz", sep="\t", index_col=0, header=[0,1]).iloc[:4,:6].to_csv(sys.stdout, sep="\t")'

Taxonomy	ESLC-1	ESLC-1	ESLC-1	ESLC-1	ESLC-1	ESLC-1
Functional	M00007	M00165	M00167	M00004	M00580	M00005
id_sample
S4	0.25	0.2727272727272727	0.4285714285714285	0.625	0.0	1.0
S3	0.0	0.0	0.0	0.0	0.0	0.0
S1	1.0	0.9090909090909092	1.0	1.0	0.5	1.0
S2	1.0	0.9090909090909092	1.0	1.0	0.5	1.0
```

* categorical_engineered_feature.pgfc.genome_clusters.synopsis.tsv.gz

```
gzip -d -c veba_output/phylogenomic_functional_categories/categorical_engineered_feature.pgfc.genome_clusters.synopsis.tsv.gz | head -n 4

Taxonomy	Functional	initial_features	number_of_features	Taxonomy(level:0)	Functional(level:1)	scaling_factors	sum(scaling_factors)	mean(scaling_factors)	sem(scaling_factors)	std(scaling_factors)
ESLC-1	M00007	['ESLC-1_SSPC-1011', 'ESLC-1_SSPC-6518', 'ESLC-1_SSPC-6643', 'ESLC-1_SSPC-8249', 'ESLC-1_SSPC-8882', 'ESLC-1_SSPC-1835', 'ESLC-1_SSPC-4092']	7	ESLC-1	M00007	[840, 930, 1026, 663, 2073, 2034, 684]	8250	1178.5714285714287	231.03657277406032	565.92171521784
ESLC-1	M00165	['ESLC-1_SSPC-1011', 'ESLC-1_SSPC-66', 'ESLC-1_SSPC-8954', 'ESLC-1_SSPC-2255', 'ESLC-1_SSPC-7572', 'ESLC-1_SSPC-8707', 'ESLC-1_SSPC-8882', 'ESLC-1_SSPC-3604', 'ESLC-1_SSPC-7959', 'ESLC-1_SSPC-10727', 'ESLC-1_SSPC-2990', 'ESLC-1_SSPC-5016', 'ESLC-1_SSPC-11715', 'ESLC-1_SSPC-8777', 'ESLC-1_SSPC-168', 'ESLC-1_SSPC-5211', 'ESLC-1_SSPC-1835']	17	ESLC-1	M00165	[840, 2001, 1128, 1209, 1104, 1110, 2073, 1254, 786, 1005, 1269, 1002, 834, 1083, 1200, 1215, 2034]	21147	1243.9411764705883	98.19726060140057	392.7890424056023
ESLC-1	M00167	['ESLC-1_SSPC-1011', 'ESLC-1_SSPC-66', 'ESLC-1_SSPC-8954', 'ESLC-1_SSPC-8707', 'ESLC-1_SSPC-8882', 'ESLC-1_SSPC-8777', 'ESLC-1_SSPC-5016', 'ESLC-1_SSPC-11715', 'ESLC-1_SSPC-168', 'ESLC-1_SSPC-5211', 'ESLC-1_SSPC-1835']	11	ESLC-1M00167	[840, 2001, 1128, 1110, 2073, 1083, 1002, 834, 1200, 1215, 2034]	14520	1320.0	143.680326988897	454.35708824267687
```

If you want to convert these to `anndata` or `biom` you might need to collapse the `pd.MultiIndex` header.


#### 2. Loading `CategoricalEngineeredFeature` object into Python

```python
>>> import ensemble_networkx as enx
>>> import pickle
>>> with open("veba_output/phylogenomic_functional_categories/categorical_engineered_feature.pgfc.genome_clusters.pkl", "rb") as f:
...     pgfcs = pickle.load(f)
...
>>> pgfcs
=======================================
CategoricalEngineeredFeature(Name:None)
=======================================
    * Number of levels: 2
    * Memory: 48 B
    * Compiled: True
    -----------------------------------
    | Types
    -----------------------------------
    * Initial feature type: protein_cluster
    * Engineered feature type: pgfc
    * Observation feature type: id_sample
    * Unit type: normalized_counts
    -----------------------------------
    | Statistics
    -----------------------------------
        Scaling Factors: True
        Summary: ['sum', 'mean', 'sem', 'std']
        Tests: []
    -----------------------------------
    | Categories
    -----------------------------------
    * Level 0 - Taxonomy:
        Number of initial features: 4888
        Number of categories: 16
    * Level 1 - Functional:
        Number of initial features: 4888
        Number of categories: 228
    -----------------------------------
    | Features
    -----------------------------------
      Number of initial features (Intersection): 4888
      Number of initial features (Union): 4888
      Number of engineered features: 2358

# You can also view the synopsis here
>>> pgfcs.synopsis_
                                                      initial_features number_of_features Taxonomy(level:0) Functional(level:1)  ... sum(scaling_factors) mean(scaling_factors) sem(scaling_factors) std(scaling_factors)
Taxonomy Functional                                                                                                              ...
ESLC-1   M00007      [ESLC-1_SSPC-1011, ESLC-1_SSPC-6518, ESLC-1_SS...                  7            ESLC-1              M00007  ...                 8250           1178.571429           231.036573           565.921715
         M00165      [ESLC-1_SSPC-1011, ESLC-1_SSPC-66, ESLC-1_SSPC...                 17            ESLC-1              M00165  ...                21147           1243.941176            98.197261           392.789042
         M00167      [ESLC-1_SSPC-1011, ESLC-1_SSPC-66, ESLC-1_SSPC...                 11            ESLC-1              M00167  ...                14520                1320.0           143.680327           454.357088
         M00004      [ESLC-1_SSPC-1011, ESLC-1_SSPC-4578, ESLC-1_SS...                 17            ESLC-1              M00004  ...                23766                1398.0           162.252051           649.008203
         M00580                                     [ESLC-1_SSPC-1011]                  1            ESLC-1              M00580  ...                  840                 840.0                  NaN                  0.0
...                                                                ...                ...               ...                 ...  ...                  ...                   ...                  ...                  ...
[2358 rows x 9 columns]

```

Check out the [documentation](https://github.com/jolespin/ensemble_networkx?tab=readme-ov-file#feature-engineering-using-categories) for this object to see how you can transform more datasets (e.g., metatranscriptomics that are paired with your metagenomics)
_____________________________________________________

#### Next steps:

Now it's time to analyze the data using compositional data analysis (CoDA).  Check out some of our other packages: 

* https://github.com/jolespin/compositional
* https://github.com/jolespin/ensemble_networkx
* https://github.com/jolespin/hive_networkx


##### Recommended reading for CoDA:

* Thomas P Quinn, Ionas Erb, Mark F Richardson, Tamsyn M Crowley, Understanding sequencing data as compositions: an outlook and review, Bioinformatics, Volume 34, Issue 16, 15 August 2018, Pages 2870–2878, https://doi.org/10.1093/bioinformatics/bty175
* Thomas P Quinn, Ionas Erb, Greg Gloor, Cedric Notredame, Mark F Richardson, Tamsyn M Crowley, A field guide for the compositional analysis of any-omics data, GigaScience, Volume 8, Issue 9, September 2019, giz107, https://doi.org/10.1093/gigascience/giz107
* Espinoza, J.L., Shah, N., Singh, S., Nelson, K.E. and Dupont, C.L. (2020), Applications of weighted association networks applied to compositional data in biology. Environ Microbiol, 22: 3020-3038. https://doi.org/10.1111/1462-2920.15091
