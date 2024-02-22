<a name="readme-top"></a>

# Quick interpretation guide for *VEBA*

This guide serves as a source for interpreting the output of `VEBA` and the various modules. It is not meant as an [exhaustive explanation of files](../../bin/README.md) but instead just meant to showcase the main output files. 
___________________________________________________________________

### preprocess/preprocess-long
The main function of this module is to preprocess your reads for downstream analysis. 

There are a few different modes for the preprocessing modules.  Below are different scenarios and usages: 

* **Trimmed reads:**

    * Only reads are provided: 
        * `trimmed_1/2.fastq.gz` or `trimmed.fastq.gz` will contain the trimmed reads
    * If a contamination database (e.g., Human reference genome):
        * `cleaned_1/2.fastq.gz` or `cleaned.fastq.gz` will contain the decontaminated trimmed reads
        * `contaminated_1/2.fastq.gz` or `cleaned.fastq.gz` will contain the contaminated trimmed reads
    * If a k-mer database is provided (e.g., Ribokmers):
        * `kmer_hits_1/2.fastq.gz` or `kmer_hits.fastq.gz` will contain all the reads with hits to k-mer database
        * `non-kmer_hits_1/2.fastq.gz` or `non-kmer_hits.fastq.gz` will contain all the reads that did not have any hits to k-mer database

* **Read statistics:**

    * `seqkit_stats.concatenated.tsv` - Concatenated read statistics for all intermediate steps (e.g., fastp, bowtie2 removal of contaminated reads if provided, bbduk.sh removal of contaminated reads if provided)

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________


### assembly/assembly-long
The main function of this module is to assemble reads and align reads back to assembly. 

There are a few different modes for the assembly modules.  Below are different scenarios and usages: 

* **Assembly:**

    By default, sample identifier is preprended to each scaffold/contig/transcript name.  Each assembly is indexed with either `Bowtie2` or `Minimap2` depending on if you are using the `assembly` or `assembly-long` module.

    * If `SPAdes`, `metaSPAdes`, or `MEGAHIT` are used:
        * `scaffolds.fasta` is the main assembly.
    * If `rnaSPAdes` is used:
        * `transcripts.fasta` is the main assembly. 
    * If `Flye` or `MetaFlye` are used (i.e., long-read assembly):
        * `assembly.fasta` is the main assembly

* **Sorted alignment:**

    * `mapped.sorted.bam` is a sorted BAM file and `mapped.sorted.bam.bai` is the `samtools` index.

* **Assembly statistics:**

    * `seqkit_stats.tsv` - Assembly statistics such as N50, length, GC%, number of contigs, etc.

* **Counts tables:**

    * `featurecounts.tsv.gz` - Contig or transcript-level counts table for a single sample

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### coverage/coverage-long

The main function of this module is align reads to an assembly/co-assembly/pseudo-coassembly. 

There are a few different modes for the assembly modules.  Below are different scenarios and usages:

* **Reference:**

    * `reference.fasta` - Reference fasta (e.g., co-assembly or pseudo-coassembly)

* **Sorted alignments:**

    * `[sample_id]/mapped.sorted.bam` contains all the sorted BAM files from the samples provided. Also contains `[sample_id]/mapped.sorted.bam.bai` for the `samtools` index.

* **Sequence statistics:**

    * `seqkit_stats.tsv` - Reference fasta statistics such as N50, length, GC%, number of sequences, etc.

* **Counts tables:**

    * `featurecounts.tsv.gz` - Contig or transcript-level counts table for all samples

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### binning-prokaryotic

The main function of this module is to identify candidate prokaryotic genomes (bins), model genes, and quality assess to yield metagenome-assembled genomes (MAG).

* **Genomes:**

    * `genomes/[id_genome].fa` - MAG assembly fasta
    * `genomes/[id_genome].faa` - MAG protein fasta
    * `genomes/[id_genome].ffn` - MAG CDS fasta
    * `genomes/[id_genome].gff` - MAG gene models for assembly, CDS, rRNA, and tRNA
    * `genomes/[id_genome].rRNA` - MAG rRNA fasta
    * `genomes/[id_genome].tRNA` - MAG tRNA fasta
    * `genomes/identifier_mapping.tsv` - Identifier mapping between [id_orf, id_contig, id_mag]


* **Assembly statistics:**

    * `genome_statistics.tsv` - Genome assembly statistics
    * `gene_statistics.cds.tsv` - Gene sequence statistics (CDS)
    * `gene_statistics.rRNA.tsv` - Gene sequence statistics (rRNA)
    * `gene_statistics.tRNA.tsv` - Gene sequence statistics (tRNA)

* **Counts tables**:

    * `featurecounts.orf.tsv.gz` - ORF-level counts table for all samples provided as `mapped.sorted.bam` files.

* **Quality assessment:** 

    * `checkm2_results.filtered.tsv` - Filtered `CheckM2` output with completeness and contamination


<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### binning-eukaryotic

The main function of this module is to identify candidate eukaryotic genomes (bins), partition nuclear genome from organelles, model genes, and quality assess to yield metagenome-assembled genomes (MAG).

* **Genomes:**

    * `genomes/[id_genome].fa` - MAG assembly fasta
    * `genomes/[id_genome].faa` - MAG protein fasta
    * `genomes/[id_genome].ffn` - MAG CDS fasta
    * `genomes/[id_genome].gff` - MAG gene models for assembly, CDS, rRNA, and tRNA
    * `genomes/[id_genome].rRNA` - MAG rRNA fasta
    * `genomes/[id_genome].tRNA` - MAG tRNA fasta
    * `genomes/identifier_mapping.tsv` - Identifier mapping between [id_orf, id_contig, id_mag]
    * `genomes/[id_genome].seqtype.tsv` - Identifier mapping between [id_contig, sequence_type] {nuclear, mitochondrion, plastid}

* **Sequence statistics:**

    * `genome_statistics.tsv` - Genome assembly statistics
    * `gene_statistics.cds.tsv` - Gene sequence statistics (CDS)
    * `gene_statistics.rRNA.tsv` - Gene sequence statistics (rRNA)
    * `gene_statistics.tRNA.tsv` - Gene sequence statistics (tRNA)

* **Counts tables:**

    * `featurecounts.orf.tsv.gz` - ORF-level counts table for all samples provided as `mapped.sorted.bam` files.

* **Quality assessment:** 

    * `busco_results.filtered.tsv` - Filtered `BUSCO` output with completeness and contamination


<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### binning-viral

The main function of this module is to identify candidate viral genomes (bins), plasmids, model genes, and quality assess to yield metagenome-assembled genomes (MAG).

* **Genomes:**

    * `genomes/[id_genome].fa` - MAG assembly fasta
    * `genomes/[id_genome].faa` - MAG protein fasta
    * `genomes/[id_genome].ffn` - MAG CDS fasta
    * `genomes/[id_genome].gff` - MAG gene models for assembly, CDS, rRNA, and tRNA
    * `genomes/identifier_mapping.tsv` - Identifier mapping between [id_orf, id_contig, id_mag]


* **Sequence statistics**:

    * `genome_statistics.tsv` - Genome assembly statistics
    * `gene_statistics.cds.tsv` - Gene sequence statistics (CDS)

* **Counts tables**:

    * `featurecounts.orf.tsv.gz` - ORF-level counts table for all samples provided as `mapped.sorted.bam` files.

* **Quality assessment:** 

    * `checkv_results.filtered.tsv` - Filtered `CheckV` output with quality, completeness, and contamination


<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### classify-prokaryotic

The main function of this module is classify prokaryotic genomes using `GTDB-Tk`.  You can use either genomes recovered using `VEBA` or acquired elsewhere (e.g., NCBI, JGI, another tool).

* **Taxonomy:**

    * `taxonomy.tsv` - Prokaryotic genome classification based on `GTDB-Tk`
    * `taxonomy.clusters.tsv` - Prokaryotic cluster classification (If `--clusters` are provided)

* **Visualization:**

    * `krona.html` - Interactive `Krona` plot for various taxonomic levels

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### classify-eukaryotic

The main function of this module is classify eukaryotic genomes using `GTDB-Tk`.  You can use either genomes recovered using `VEBA` or acquired elsewhere (e.g., NCBI, JGI, another tool).

* **Taxonomy:**

    * `taxonomy.tsv` - Eukaryotic genome classification based on `VEBA`'s classification method.
    * `taxonomy.clusters.tsv` - Eukaryotic cluster classification (If `--clusters` are provided)
    * `gene-source_lineage.tsv` - Gene source lineage and scores for classifying MAGs [id_gene, id_scaffold, id_mag, id_target, id_source, lineage, bitscore]

* **Visualization:**

    * `krona.html` - Interactive `Krona` plot for various taxonomic levels


<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________


### classify-viral

The main function of this module is classify viral genomes using `geNomad`.  You can use either genomes recovered using `VEBA` or acquired elsewhere (e.g., NCBI, JGI, another tool).

* **Taxonomy:**

    * `taxonomy.tsv` - Viral genome classification based on `geNomad`
    * `taxonomy.clusters.tsv` - Viral cluster classification (If `--clusters` are provided)

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### cluster

The main function of this module is cluster genomes using `Skani|FastANI` and proteins within those pangenomes using `MMseqs2|DeepClust`.  You can use either genomes recovered using `VEBA` or acquired elsewhere (e.g., NCBI, JGI, another tool).  `global` refers to clustering across the entire dataset and `local` refers to clustering within each sample individually.

* **Identifier mapping:**

    * `identifier_mapping.genomes.tsv` - Identifier mapping for genomes [id_genome, organism_type, sample_of_origin, id_genome_cluster, number_of_proteins, number_of_singleton_protein_clusters, ratio_of_protein_cluster_are_singletons]
    * `identifier_mapping.proteins.tsv` - Identifier mapping for proteins [id_protein, organism_type, id_genome, sample_of_origin, id_genome_cluster, id_protein_cluster]
    * `identifier_mapping.scaffolds.tsv` - Identifier mapping for contigs [id_scaffold,	organism_type, id_genome, sample_of_origin, id_genome_cluster]

* **Metadata:**

    * `genome_clusters.tsv` - Machine-readable table for genome clusters [id_genome_cluster, number_of_components, number_of_samples_of_origin, components, samples_of_origin]
    * `protein_clusters.tsv` - Machine-readable table for protein clusters [id_protein_cluster, number_of_components, number_of_samples_of_origin, components, samples_of_origin]

* **Representatives:**

    * `representative_sequences.faa` - Protein sequences for cluster representatives. Header follows the following format: id_protein-cluster id_original_protein

* **Pangenome tables:**

    Prevalence tables for each pangenome showing genomes vs. protein clusters

    * `pangenome_tables/*.tsv.gz` - Pangenome tables for each SLC with prevalence values

* **ANI networks:**

    Networks of genomes with ANI-based edges

    * `serialization/*.networkx_graph.pkl` - `NetworkX` graphs for clusters

    You can open these in Python using the following code: 

    ```python
    import pickle
    import networkx as nx
    with open("serialization/prokaryotic.networkx_graph.pkl", "rb") as f: 
        graph = pickle.load(f)
    ```
    Each organism type has its own graph so you'll need to replace with `viral.networkx_graph.pkl` and `eukaryotic.networkx_graph.pkl` if you want to load those instead. 

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________


### index

The main function of this module is to build a `Bowtie2` index that is used by the `mapping` module. There are 2 modes, global and local.  For specifics on the differences between them, please refer to the following [walkthrough](read_mapping_and_counts_tables.md)

* **Reference sequences:**

    * `reference.fa.gz` - Concatenated reference fasta
    * `reference.saf` - SAF format for reference

* **Bowtie2 index:**

    * `reference.fa.gz.*.bt2` - Bowtie2 index of reference fasta

* **Gene models:**

    * `reference.gff` - Concatenated gene models

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### mapping

The main function of this module is map reads to `Bowtie2` index created by `index` module.

* **Sorted alignments:**

    * `mapped.sorted.bam` - Sorted BAM file
    * `mapped.sorted.bam.bai` - Sorted BAM file index

* **Coverage:**
    * `mapped.sorted.bam.coverage.tsv.gz` - Samtools coverage table
    * `genome_spatial_coverage.tsv.gz` - Spatial coverage for genome (i.e., ratio of bases covered) [Only if `--scaffolds_to_bins` is provided]

* **Counts tables**:

    * `counts.orfs.tsv.gz` - ORF-level counts table
    * `counts.scaffolds.tsv.gz` - Contig-level counts table
    * `counts.mags.tsv.gz` - MAG-level counts table [Only if `--scaffolds_to_bins` is provided]
    * `counts.clusters.tsv.gz` - SLC-level counts table [Only if `--scaffolds_to_clusters` is provided]
    * `counts.orthogroups.tsv.gz` - Orthogroup-level counts table [Only if `--orf_to_orthogroups` is provided]

* **Unmapped reads (Fastq)**:

    * `unmapped_1.fastq.gz` - Unmapped reads (forward)
    * `unmapped_2.fastq.gz` - Unmapped reads (reverse)

    To merge these counts tables together, run the following command from this [walkthrough](read_mapping_and_counts_tables.md): 

    ```bash
    # Merge contig-level counts
    merge_contig_mapping.py -m ${MAPPING_DIRECTORY} -c ${MAGS_TO_SLCS}  -i ${SCAFFOLDS_TO_MAGS} -o ${OUT_DIR}

    # Merge ORF-level counts
    merge_orf_mapping.py -m ${MAPPING_DIRECTORY} -c ${PROTEINS_TO_ORTHOGROUPS} -o ${OUT_DIR}
    ```

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### profile-taxonomy

The main function of this module is estimate taxonomic abundances of a custom genome database using `Sylph`.

* **Abundance tables**:

    * `sylph_profile.tsv.gz` - Output of `sylph profile`
    * `taxonomic_abundance.tsv.gz` - Genome-level taxonomic abundance (No header)
    * `taxonomic_abundance.clusters.tsv.gz` - SLC-level taxonomic abundance (No header, if `--genome_clusters` were provided)

    To merge these abundance tables together, run the following command from this [walkthrough](taxonomic_profiling_de-novo_genomes.md)

    ```bash
    merge_generalized_mapping.py -o veba_output/profiling/taxonomy/merged.taxonomic_abundance.tsv.gz veba_output/profiling/taxonomy/*/output/taxonomic_abundance.tsv.gz

    merge_generalized_mapping.py -o veba_output/profiling/taxonomy/merged.taxonomic_abundance.clusters.tsv.gz veba_output/profiling/taxonomy/*/output/taxonomic_abundance.clusters.tsv.gz
    ```

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### profile-pathway

The main function of this module is estimate genome-resolved pathway abundances of a custom genome database using `HUMAnN`.

* **Abundance tables**:

    * `humann_pathabundance.tsv` - Stratified abundance of taxa-specific metabolic pathways
    * `humann_pathcoverage.tsv` - Stratified pathway completion ratio of taxa-specific metabolic pathways
    * `humann_genefamilies.tsv` - Stratified abundance of taxa-specific gene families

    To merge these abundance tables together, run the following command from this [walkthrough](pathway_profiling_de-novo_genomes.md)

    ```bash
    merge_generalized_mapping.py -o veba_output/profiling/pathways/merged.humann_pathcoverage.tsv veba_output/profiling/pathways/*/output/humann_pathcoverage.tsv

    merge_generalized_mapping.py -o veba_output/profiling/pathways/merged.humann_pathabundance.tsv veba_output/profiling/pathways/*/output/humann_pathabundance.tsv
    
    merge_generalized_mapping.py -o veba_output/profiling/pathways/merged.humann_genefamilies.tsv veba_output/profiling/pathways/*/output/humann_genefamilies.tsv
    ```

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### phylogeny

The main function of this module is build phylogenetic trees from concatenated protein alignments. 

* **Phylogenetic trees:**

    * `concatenated_alignment.fasttree.nw` or `concatenated_alignment.veryfasttree.nw` newick format trees  based on concatenated alignment if `FastTree` or `VeryFastTree`, respectively.

    * `output.treefile` - `IQTREE2` newick format tree based on concatenated alignment (if `--no_iqtree` is not selected).  This tree takes a very long time to compute if there are lots genomes so only use this if you need extremely robust trees.


* **Multiple sequence alignment:**

    * `concatenated_alignment.fasta` - Concatenated protein alignment of all marker hits

* **Visualization:**

    * `[newick].pdf` newick trees above are also rendered into PDF files for quick visualization.


<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________

### biosynthetic

The main function of this module is run `antiSMASH`, reformat the results (tabular and fasta), then cluster the BGCs in both protein and nucleotide-space.

* **Identifier mapping:**

    * `identifier_mapping.bgcs.tsv.gz` - All of the BGCs in tabular format organized by genome, contig, region, and gene.
    * `identifier_mapping.components.tsv.gz` - All of the BGC components (i.e., genes in BGC) in tabular format organized by genome, contig, region, and gene.
    * `bgc_protocluster-types.tsv.gz` - Summary of BGCs detected organized by type.  Also includes summary of BGCs that are NOT on contig edge.

* **Sequences:**

    * `bgcs.representative_sequences.fasta.gz` - Full length BGC nucleotide cluster representatives
    * `components.representative_sequences.faa.gz` - BGC protein cluster representatives
    * `fasta/[id_genome].faa/fasta.gz` - BGC sequences in protein and nucleotide space
    * `genbanks/[id_genome]/*.gbk` - Genbank formatted antiSMASH results

* **Protein homology:**

    * `homology.tsv.gz` - Diamond results for MIBiG and VFDB

* **Visualization:**

    * `krona.html` - HTML showing Krona plot for number of BGCs per protocluster-type.

* **Prevalence tables:**

    Prevalence tables in protein and nucleotide-space.  

    * `prevalence_tables/bgcs.tsv.gz` - Genome vs. BGC nucleotide cluster prevalence table
    * `prevalence_tables/components.tsv.gz` - Genome vs. BGC protein cluster prevalence table

    These can be used for either clustering BGCs based on protein homology or samples based on the presence of BGCs. Since the prevalence tables are boolean, it is recommended to use [Jaccard distance](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.jaccard.html) with a clustering algorithm that allows for a [precomputed distance matrix](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.HDBSCAN.html#sklearn.cluster.HDBSCAN). 

    ```python
    import pandas as pd
    from scipy.spatial.distance import pdist, squareform
    from sklearn.cluster import HDBSCAN

    # Load prevalence table of boolean values
    df_prevalence = pd.read_csv("prevalence_tables/components.tsv.gz", sep="\t", index_col=0)

    # Calculate pairwise jaccard distances
    pairwise_distances = squareform(pdist(df_prevalence.values, metric="jaccard"))

    # Clustering using precomputed distances
    clustering = HDBSCAN(metric="precomputed", min_cluster_size=2)
    clustering.fit(pairwise_distances)

    # Get clusters and reassign names
    clusters = pd.Series(data=clustering.labels_, index=df_prevalence.index)
    ```


<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________