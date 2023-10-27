### Bioprospecting for biosynthetic gene clusters
This walkthrough goes through identification of biosynthetic gene clusters (BGC).

What you'll end up with at the end of this are tables with respect to genes/components in each BGC, synopsis of each BGC, homology hits to `MIBiG`, novelty scores for each BGC, and the fasta sequences for all the proteins involved.

It is assumed you've completed the [end-to-end metagenomics](end-to-end_metagenomics.md) walkthrough.

_____________________________________________________

#### Steps:

1. Compile table of genomes and gene models
2. Identify biosynthetic gene clusters and score novelty

_____________________________________________________


#### 1. Compile table of genomes and gene models
First we need to compile our genomes and gene models into a table.

```
compile_genomes_table.py -i veba_output/binning/ > veba_output/misc/genomes_table.tsv
```

This results in a table with the following columns: 

`[domain] [id_sample] [id_genome] [path/to/genome.fasta] [path/to/proteins.fasta] [path/to/cds.fasta] [path/to/gene_models.gff]`

If you [ran the clustering from the end-to-end metagenomics walkthrough](https://github.com/jolespin/veba/blob/main/walkthroughs/end-to-end_metagenomics.md#9-cluster-genomes-and-proteins) then you should already have this file: `veba_output/misc/genomes_table.tsv`. 

We only need the `[id_genome] [path/to/genome.fasta] [path/to/gene_models.gff]` columns so let's get some cut up for just the prokaryotes (note, this could be a one-shot with the `compile_genomes_table.py` command via piping):

`cat veba_output/misc/genomes_table.tsv | grep "^prokaryotic" | cut -f3,4,7 > veba_output/misc/genomes_gene-models.tsv`


#### 2. Identify biosynthetic gene clusters and score novelty

Now that we have our genome table formatted so it is `[id_genome] [path/to/genome.fasta] [path/to/gene_models.gff]` without headers, we can run the `biosynthetic.py` module to identify biosynthetic gene clusters via `antiSMASH` and detect homology of components to the `MIBiG` database.

**Conda Environment:** `conda activate VEBA-biosynthetic_env`


```
# Set the number of threads
N_JOBS=16

N=biosynthetic-prokaryotic

rm -f logs/${N}.*

# Proteome file paths
GENOMES=veba_output/misc/genomes_gene-models.tsv

# Output directory
OUT_DIR=veba_output/biosynthetic/prokaryotic

# Directory
CMD="source activate VEBA-biosynthetic_env && biosynthetic.py -i ${GENOMES} -o ${OUT_DIR} -p ${N_JOBS} -t bacteria"

# Either run this command or use SunGridEnginge/SLURM
```

The following output files will produced: 

* identifier\_mapping.bgcs.tsv.gz - Identifier mapping and synopsis for BGCs with identifiers, bgc_type, whether or not cluster is on contig edge, number of genes, number of hits to MIBiG, and novety score.
* identifier\_mapping.components.tsv.gz - Identifier mapping and synopsis for `antiSMASH` genbank files for BGCs compiled into a table format which includes all information included in the genbank files in addition to identifiers for BGC, component, genome, contig, region, gene, and BGC-type.  BGC identifiers are formatted as the following: `[id_genome]|[id_contig]|[id_region]` and component identifiers follow the format: `[id_genome]|[id_contig]|[id_region]_[position_of_gene_in_BGC]|[start_on_contig]-[end_on_contig]([strand])`. One-to-one mapping for component id and gene id called via `Prodigal` and augmented with `VEBA`.  Rows are with respect to gene/component.
* homology.tsv.gz - Diamond alignment results to `MIBiG` and `VFDB` as database.
* krona.html - Krona plot of BGC protocluster types
* genbanks/*.gbk - BGC genbank files from `antiSMASH`
* fasta/[id\_genome].faa - BGCs in protein space. Fasta header following: `[id_gene] [id_component]`.
* fasta/[id\_genome].fasta - BGCs in nucleotide space. Fasta header following: `[id_bgc] len={},gc={},n_genes={},edge={}`.
* prevalence_tables/bgcs.tsv.gz - Prevalence tables (row=genome, columns=BGC nucleotide cluster)
* prevalence_tables/components.tsv.gz - Prevalence tables (row=BGC, columns=BGC protein cluster)
* bgc_clusters.tsv - `MMSEQS2` clustering in nucleotide space
* component_clusters.tsv - `MMSEQS2` clustering in protein space

#### Next steps:

Cluster the prevalence tables, synthesize products, preserve the ecosystem, and save humanity.