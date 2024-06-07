### What's next for *VEBA*?

*VEBA* is currently under active development. If you are interested in requesting features or wish to report a bug, please post a GitHub issue prefixed with the tag `[Feature Request]` and `[Bug]`, respectively.  If you want to contribute or have any other inquiries, contact me at `jol.espinoz[A|T]gmail[DOT]com`.

________________________________________________________________

#### Current Releases:

**Release v2.0.0 Highlights:**

* Added `-A/--from_antismash` in `biosynthetic.py` to use preexisting `antiSMASH` results.  Also changed `-i/--input` to `-i/--from_genomes`.
* Added `number_of_genomes`, `number_of_genome-clusters`, `number_of_proteins`, and `number_of_protein-clusters` to `feature_compression_ratios.tsv.gz` from `cluster.py`
* Added custom path for `conda` environments
* Added `busco_version` parameter to `merge_busco_json.py` with default set to `5.4.x` and additional support for `5.6.x`.
* Changed `antimash_genbanks_to_table.py` to `biosynthetic_genbanks_to_table.py` for future support of `DeepBGC` and `GECCO`


<details>
	<summary><b>Release v2.0.0 Details</b></summary>

* Changed default assembly algorithm to `metaflye` instead of `flye` in `assembly-long.py`
* Added `number_of_genomes`, `number_of_genome-clusters`, `number_of_proteins`, and `number_of_protein-clusters` to `feature_compression_ratios.tsv.gz` from `cluster.py`
* Added `-A/--from_antismash` in `biosynthetic.py` to use preexisting `antiSMASH` results.  Also changed `-i/--input` to `-i/--from_genomes`.
* Changed `antimash_genbanks_to_table.py` to `biosynthetic_genbanks_to_table.py` for future support of `DeepBGC` and `GECCO`
* Added `busco_version` parameter to `merge_busco_json.py` with default set to `5.4.x` and additional support for `5.6.x`.
* Added `CONDA_ENVS_PATH` to `update_environment_scripts.sh`, `update_environment_variables.sh`, and `check_installation.sh`
* Added `CONDA_ENVS_PATH` to `veba` to allow for custom environment locations
* Changed `install.sh` to support custom `CONDA_ENVS_PATH` argument `bash install.sh path/to/log path/to/envs/`
*  Added `merge_counts_with_taxonomy.py`


</details>

**Release v1.5.0 Highlights:**

* Added `VeryFastTree` to `phylogeny.py`
* Added `--blacklist` to `compile_eukaryotic_classifications.py`
* Added compatibility for `antismash_genbanks_to_table.py` to operate on `antiSMASH v7` genbanks
* Added `compile_phylogenomic_functional_categories.py` script which automates the methodology from [Espinoza et al. 2022 (doi:10.1093/pnasnexus/pgac239)](https://academic.oup.com/pnasnexus/article/1/5/pgac239/6762943)
* Fixed error in `annotations.protein_clusters.tsv` formatting from `annotate.py`
* Fixed situation where `unbinned.fasta` were not added in `binning-prokaryotic.py` and bad symlinks were created for GFF, rRNA, and tRNA when no genoems were detected.
* Fixed critical error where `classify_eukaryotic.py` was trying to access a deprecated database file from MicroEuk_v2.

<details>
	<summary><b>Release v1.5.0 Details</b></summary>

* Cleaned up installation files
* Changed `veba/src/` to `veba/bin/`
* Checked `SCRIPT_VERSIONS` to `VEBA_SCRIPT_VERSIONS` which are now in `bin/` of conda environment
* Fixed header being offset in `annotations.protein_clusters.tsv` where it could not be read with Pandas.
* Fixed `binning-prokaryotic.py` the creation of non-existing symlinks where "'*.gff'", "'*.rRNA'", and "'*.tRNA'" were created.
* Fixed .strip method on Pandas series in `antismash_genbanks_to_table.py` for compatibilty with `antiSMASH 6 and 7`
* Fixed situation where `unbinned.fasta` is empty in `binning-prokaryotic.py` when there are no bins that pass qc.
* Fixed minor error in `coverage.py` where `samtools sort --reference` was getting `reads_table.tsv` and not `reference.fasta`
* Changed default behavior from deterministic to not deterministic for increase in speed in `assembly-long.py`. (i.e., `--no_deterministic` --> `--deterministic`)
* Added `VeryFastTree` as an option to `phylogeny.py` with `FastTree` remaining as the default.
* Changed default `--leniency` parameter on `classify_eukaryotic.py` and `consensus_genome_classification_ranked.py` to `1.0` and added `--leniecy_genome_classification` as a separate option.
* Added `--blacklist` option to `compile_eukaryotic_classifications.py` with a default value of `species:uncultured eukaryote` in `classify_eukaryotic.py`
* Fixed critical error where `classify_eukaryotic.py` was trying to access a deprecated database file from MicrEuk_v2.
* Fixed minor error with `eukaryotic_gene_modeling_wrapper.py` not allowing for `Tiara` to run in backend.
* Added `compile_phylogenomic_functional_categories.py` script which automates the methodology from [Espinoza et al. 2022 (doi:10.1093/pnasnexus/pgac239)](https://academic.oup.com/pnasnexus/article/1/5/pgac239/6762943)
</details>


**Release v1.4.2 Highlights:**

* **`VEBA` Modules:**

	* Added `profile-taxonomic.py` module which uses `sylph` to build a sketch database for genomes and queries the genome database for taxonomic abundance.
	* Added long read support for `fastq_preprocessor`, `preprocess.py`, `assembly-long.py`, `coverage-long`, and all binning modules.
	* Redesign `binning-eukaryotic` module to handle custom `MetaEuk` databases
	* Added new usage syntax `veba --module preprocess --params “${PARAMS}”` where the Conda environment is abstracted and determined automatically in the backend.  Changed all the walkthroughs to reflect this change.
	* Added `skani` which is the new default for genome-level clustering based on ANI.
	* Added `Diamond DeepClust` as an alternative to `MMSEQS2` for protein clustering.

* **`VEBA` Database (`VDB_v6`)**:

	* Completely rebuilt `VEBA's Microeukaryotic Protein Database` to produce a clustered database `MicroEuk100/90/50` similar to `UniRef100/90/50`. Available on [doi:10.5281/zenodo.10139450](https://zenodo.org/records/10139451).

<details>
	<summary><b>Release v1.4.2 Details</b></summary>

* Fixed critical error where `classify_eukaryotic.py` was trying to access a deprecated database file from MicrEuk_v2.
* Added `profile-taxonomic.py` module which uses `sylph` to build a sketch database for genomes and queries the genome database similar to `Kraken` for taxonomic abundance.
* Removed requirement to have `--estimated_assembly_size` for Flye per [Flye Issue #652](https://github.com/fenderglass/Flye/issues/652).
* Added `sylph` to `VEBA-profile_env` for abundance profiling of genomes.
* Dereplicate duplicate contigs in `concatenate_fasta.py`.
* Added `--reference_gzipped` to `index.py` and `mapping.py` with new default being that the reference fasta is not gzipped.
* Added `skani` as new default for genome clustering in `cluster.py`, `global_clustering.py`, and `local_clustering.py`.
* Added support for long reads in `fastq_preprocessor`, `preprocess.py`, `assembly-long.py`, `coverage-long`, and all binning modules.
* Fixed `annotations.protein_clusters.tsv.gz` from `merge_annotations.py` added in patch update of `v1.3.1`.
* Added support for missing values in `compile_eukaryotic_classifications.py`.
* Added `--metaeuk_split_memory_limit` argument with (experimental) default set to `36G` in `binning-eukaryotic.py` and `eukaryotic_gene_modeling.py`.
* Added `--compressed 1` to `mmseqs createdb` in `download_databases.sh` installation script.
* Added a check to `check_fasta_duplicates.py` and `clean_fasta.py` to make sure there are no `>` characters in fasta sequence caused from concatenating fasta files that are missing linebreaks.
* Added `Diamond DeepClust` to `clustering_wrapper.py`, `global/local_clustering.py`, and `cluster.py`.  Changed `mmseqs2_wrapper.py` to `clustering_wrapper.py`.  Changed `easy-cluster` and `easy-linclust` to `mmseqs-cluster` and `mmseqs-linclust`.
* Fixed viral quality in `merge_genome_quality_assessments.py`
* Changed `consensus_genome_classification.py` to `consensus_genome_classification_ranked.py`.  Also, default behavior to allow for missing taxonomic levels.
* Fixed the `merge_annotations.py` resulting in a memory leak when creating the `annotations.protein_clusters.tsv.gz` output table.  However, still need to correct the formatting for empty sets and string lists.

</details>

**Release v1.3.0 Highlights:**

* **`VEBA` Modules:**
	* Added `profile-pathway.py` module and associated scripts for building `HUMAnN` databases from *de novo* genomes and annotations.  Essentially, a reads-based functional profiling method via `HUMAnN` using binned genomes as the database.
	* Added `marker_gene_clustering.py` script which identifies core marker proteins that are present in all genomes within a genome cluster (i.e., pangenome) and unique to only that genome cluster.  Clusters in either protein or nucleotide space.
	* Added `module_completion_ratios.py` script which calculates KEGG module completion ratios for genomes and pangenomes. Automatically run in backend of `annotate.py`.
	* Updated `annotate.py` and `merge_annotations.py` to provide better annotations for clustered proteins.
	* Added `merge_genome_quality.py` and `merge_taxonomy_classifications.py` which compiles genome quality and taxonomy, respectively, for all organisms.
	* Added BGC clustering in protein and nucleotide space to `biosynthetic.py`.  Also, produces prevalence tables that can be used for further clustering of BGCs.
	* Added `pangenome_core_sequences` in `cluster.py` writes both protein and CDS sequences for each genome cluster.
	* Added PDF visualization of newick trees in `phylogeny.py`.

	
* **`VEBA` Database (`VDB_v5.2`)**:
	* Added `CAZy`
	* Added `MicrobeAnnotator-KEGG`

<details>
	<summary><b>Release v1.3.0 Details</b></summary>
	
* Update `annotate.py` and `merge_annotations.py` to handle `CAZy`.  They also properly address clustered protein annotations now. 
* Added `module_completion_ratio.py` script which is a fork of `MicrobeAnnotator` [`ko_mapper.py`](https://github.com/cruizperez/MicrobeAnnotator/blob/master/microbeannotator/pipeline/ko_mapper.py).  Also included a database [Zenodo: 10020074](https://zenodo.org/records/10020074) which will be included in `VDB_v5.2`
* Added a checkpoint for `tRNAscan-SE` in `binning-prokaryotic.py` and `eukaryotic_gene_modeling_wrapper.py`.
* Added `profile-pathway.py` module and `VEBA-profile_env` environments which is a wrapper around `HUMAnN` for the custom database created from `annotate.py` and `compile_custom_humann_database_from_annotations.py`
* Added `GenoPype version` to log output
* Added `merge_genome_quality.py` which combines `CheckV`, `CheckM2`, and `BUSCO` results.
* Added `compile_custom_humann_database_from_annotations.py` which compiles a `HUMAnN` protein database table from the output of `annotate.py` and taxonomy classifications.
* Added functionality to `merge_taxonomy_classifications.py` to allow for `--no_domain` and `--no_header` which will serve as input to `compile_custom_humann_database_from_annotations.py`
* Added `marker_gene_clustering.py` script which gets core marker genes unique to each SLC (i.e., pangenome). `average_number_of_copies_per_genome` to protein clusters.
* Added `--minimum_core_prevalence` in `global_clustering.py`, `local_clustering.py`, and `cluster.py` which indicates prevalence ratio of protein clusters in a SLC will be considered core.  Also remove `--no_singletons` from `cluster.py` to avoid complications with marker genes.  Relabeled `--input` to `--genomes_table` in clustering scripts/module.
* Added a check in `coverage.py` to see if the `mapped.sorted.bam` files are created, if they are then skip them.  Not yet implemented for GNU parallel option.
* Changed default representative sequence format from table to fasta for `mmseqs2_wrapper.py`.
* Added `--nucleotide_fasta_output` to `antismash_genbank_to_table.py` which outputs the actual BGC DNA sequence. Changed `--fasta_output` to `--protein_fasta_output` and added output to `biosynthetic.py`. Changed BGC component identifiers to `[bgc_id]_[position_in_bgc]|[start]:[end]([strand])` to match with `MetaEuk` identifiers. Changed `bgc_type` to `protocluster_type`.  `biosynthetic.py` now supports GFF files from `MetaEuk` (exon and gene features not supported by `antiSMASH`).  Fixed error related to `antiSMASH` adding CDS (i.e., `allorf_[start]_[end]`) that are not in GFF so `antismash_genbank_to_table.py` failed in those cases. 
* Added `ete3` to `VEBA-phylogeny_env.yml` and automatically renders trees to PDF.
* Added presets for `MEGAHIT` using the `--megahit_preset` option. 
* The change for using `--mash_db` with `GTDB-Tk` violated the assumption that all prokaryotic classifications had a `msa_percent` field which caused the cluster-level taxonomy to fail.  `compile_prokaryotic_genome_cluster_classification_scores_table.py` fixes this by uses `fastani_ani` as the weight when genomes were classified using ANI and `msa_percent` for everything else.  Initial error caused unclassified prokaryotic for all cluster-level classifications.
* Fixed small error where empty gff files with an asterisk in the name were created for samples that didn't have any prokaryotic MAGs.
* Fixed critical error where descriptions in header were not being removed in `eukaryota.scaffolds.list` and did not remove eukaryotic scaffolds in `seqkit grep` so `DAS_Tool` output eukaryotic MAGs in `identifier_mapping.tsv` and `__DASTool_scaffolds2bin.no_eukaryota.txt`
* Fixed `krona.html` in `biosynthetic.py` which was being created incorrectly from `compile_krona.py` script.
* Create `pangenome_core_sequences` in `global_clustering.py` and `local_clustering.py` which writes both protein and CDS sequences for each SLC.  Also made default in `cluster.py` to NOT do local clustering switching `--no_local_clustering` to `--local_clustering`.
* `pandas.errors.InvalidIndexError: Reindexing only valid with uniquely valued Index objects` in `biosynthetic.py` when `Diamond` finds multiple regions in one hit that matches.  Added `--sort_by` and `--ascending` to `concatenate_dataframes.py` along with automatic detection and removal of duplicate indices.  Also added `--sort_by bitscore` in `biosynthetic.py`.
* Added core pangenome and singleton hits to clustering output
* Updated `--megahit_memory` default from 0.9 to 0.99
* Fixed error in `genomad_taxonomy_wrapper.py` where `viral_taxonomy.tsv` should have been `taxonomy.tsv`.
* Fixed minor error in `assembly.py` that was preventing users from using `SPAdes` programs that were not `spades.py`, `metaspades.py`, or `rnaspades.py` that was the result of using an incorrect string formatting.
* Updated `bowtie2` in preprocess, assembly, and mapping modules.  Updated `fastp` and `fastq_preprocessor` in preprocess module.

</details>


**Release v1.2.0 Highlights:**

* **`VEBA` Modules:**
	* Updated `GTDB-Tk` now uses `Mash` for ANI screening to speed up classification (now provided in `VDB_v5.1` database)
	* rRNA and tRNA are identified for prokaryotic and eukaryotic genomes via `BARRNAP` and `tRNAscan-SE`
	* Eukaryotic genes (CDS, rRNA, tRNA) are analyzed separately for nuclear, mitochondrion, and plastid sequences
	* Genome GFF files include contigs, CDS, rRNA, and tRNA with tags for mitochondrion and plastids when applicable
	* Clustering automatically generates pangenome protein prevalence tables for each genome cluster
	* Ratios of singletons in each genome are now calculated
	* [Virulence factor database](http://www.mgc.ac.cn/VFs/main.htm) (`VFDB`) is now included in annotations
	* [UniRef50/90](https://www.uniprot.org/help/uniref) is now included in annotations
	* `Krona` plots are generated for taxonomy classifications and biosynthetic gene cluster detection
	* Fixed a minor issue in `biosynthetic.py` where the fasta and genbank files were not properly symlinked.  Also added virulence factor results to synopsis.
	
	
* **`VEBA` Database (`VDB_v5.1`) **:
	* Added `VFDB`
	* Updated `GTDB v207_v2 → v214.1`
	* Changed `NR  → UniRef50/90` 
	* Deprecated [`RefSeq non-redundant`](https://www.ncbi.nlm.nih.gov/refseq/about/nonredundantproteins/) in place of `UniRef`

<details>
	<summary><b>Release v1.2.0 Details</b></summary>

* Fixed minor error in `binning-prokaryotic.py` where the `--veba_database` argument wasn't utilized and only the environment variable `VEBA_DATABASE` could be used.
* Updated the Docker images to have `/volumes/input`, `/volumes/output`, and `/volumes/database` directories to mount. 
* Replaced `prodigal` with `pyrodigal` as it is faster and under active development.
* Added support for missing classifications in `compile_krona.py` and `consensus_genome_classification.py`.
* Updated `GTDB-Tk` from version `2.1.3` → `2.3.0` and `GTDB` from version `r202_v2` → `r214`. Changed `${VEBA_DATABASE}/Classify/GTDBTk` → `${VEBA_DATABASE}/Classify/GTDB`.  Added `gtdb_r214.msh` to `GTDB` database for ANI screening.
* Added pangenome and singularity tables to `cluster.py` (and associated global/local clustering scripts) to output automatically.
* Added `compile_gff.py` to merge CDS, rRNA, and tRNA GFF files.  Used in `binning-prokaryotic.py` and `binning-viral.py`.  `binning-eukaryotic.py` uses the source of this in the backend of `filter_busco_results.py`. Includes GC content for contigs and various tags. 
* Updated `BUSCO v5.3.2 -> v5.4.3` which changes the json output structure and made the appropriate changes in `filter_busco_results.py`.
* Added `eukaryotic_gene_modeling_wrapper.py` which 1) splits nuclear, mitochondrial, and plastid genomes; 2) performs gene modeling via `MetaEuk` and `Pyrodigal`; 3) performs rRNA detection via `BARRNAP`; 4) performs tRNA detection via `tRNAscan-SE`; 5) merges processed GFF files; and 5) calculates sequences statistics. 
* Added `gene_biotype=protein_coding` to `P(y)rodigal(-GV)` GFF output. 
* Added `VFDB` to `annotate.py` and database.
* Compiled and pushed `gtdb_r214.msh` mash file to [Zenodo:8048187](https://zenodo.org/record/8048187) which is now used by default in `classify-prokaryotic.py`.  It is now included in `VDB_v5.1`.
* Cleaned up global and local clustering intermediate files.  Added pangenome tables and singelton information to outputs.

</details>


<details>
	<summary><b>Release v1.1.2 Details</b></summary>
	
* Created Docker images for all modules
* Replaced all absolute path symlinks with relative symlinks
* Changed `prokaryotic_taxonomy.tsv` and `prokaryotic_taxonomy.clusters.tsv` in `classify-prokaryotic.py` (along with eukaryotic and viral) files to `taxonomy.tsv` and `taxonomy.clusters.tsv` for uniformity.
* Updating all symlinks to relative links (also in `fastq_preprocessor`) to prepare for dockerization and updating all environments to use updated GenoPype 2023.4.13.
* Changed `nr` to `uniref` in `annotate.py` and added `propagate_annotations_from_representatives.py` script while simplifying `merge_annotations_and_taxonomy.py` to `merge_annotations.py` and excluding taxonomy operations.
* Changed `nr` to `UniRef90` and `UniRef50` in `VDB_v5` 
* Changed `orfs_to_orthogroups.tsv` to `proteins_to_orthogroups.tsv` for consistency with the `cluster.py` module.  Will eventually find some consitency with `scaffolds_to_bins/scaffolds_to_mags` but this will be later.
* Added a `scaffolds_to_mags.tsv` in the clustering output.
* Added `convert_counts_table.py` which converts a counts table (and metadata) to [Pandas pickle](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_pickle.html), [Anndata h5ad](https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.write_h5ad.html), or [Biom hdf5](https://biom-format.org/documentation/generated/biom.table.Table.to_hdf5.html#biom.table.Table.to_hdf5)
* Fixed output directory for `mapping.py` which now uses `output_directory/${NAME}` structure like `binning-*.py`.
* Removed "python" prefix for script calls and now uses shebang in script for executable. Also added single paranthesis around script filepath (e.g., `'[script_filepath]'`) to escape characters/spaces in filepath.
* Added support for `index.py` to accept individual `--references [file.fasta]` and `--gene_models [file.gff]`.
* Added `stdin` support for `scaffolds_to_bins.py` along with the ability to input genome tables [id_genome]<tab>[filepath].  Also added progress bars.
* As a result of [issues/22](https://github.com/jolespin/veba/issues/22), `assembly.py`, `assembly-sequential.py`, `binning-*.py`, and `mapping.py` will use `-p --countReadPairs` for `featureCounts` and updates `subread 2.0.1 → subread 2.0.3`.  For `binning-*.py`, long reads can be used with the `--long_reads` flag.
* Updated `cluster.py` and associated `global_clustering.py`/`local_clustering.py` scripts to use `mmseqs2_wrapper.py` which now automatically outputs representative sequences.  
* Added `check_fasta_duplicates.py` script that gives `0` and `1` exit codes for fasta without and with duplicates, respectively.  Added `reformat_representative_sequences.py` to reformat representative sequences from `MMSEQS2` into either a table or fasta file where the identifers are cluster labels.  Removed `--dbtype` from `[global/local]_clustering.py`.  Removed appended prefix for `.graph.pkl` and `dict.pkl` in `edgelist_to_clusters.py`.  Added `mmseqs2_wrapper.py` and `hmmer_wrapper.py` scripts.
* Added an option to `merge_generalized_mapping.py` to include the sample index in a filepath and also an option to remove empty features (useful for Salmon).  Added an `executable='/bin/bash'` option to the `subprocess.Popen` calls in `GenoPype` to address [issues/23](https://github.com/jolespin/veba/issues/23).
* Added `genbanks/[id_genome]/` to output directory of `biosynthetic.py` which has symlinks to all the BGC genbanks from `antiSMASH`.

</details>

<details>
	<summary><b>Release v1.1.1 Details</b></summary>

* Most important update includes fixing a broken VEBA-`binning-viral.yml` install recipe which had package conflicts for `aria2` 30e8b0a.
* Fixes on conda-related environment variables in the install scripts.
* Added `MIBiG` to database and `annotate.py`
* Added a composite label for annotations in `annotate.py`
* Added `--dastool_minimum_score` to `binning-prokaryotic.py` module
* Added a wrapper around `STAR` aligner
* Updated `merge_generalized_mapping.py` script to take in BAM files instead of being dependent on a specific directory.
* Added option to have no header in `subset_table.py`

</details>

<details>
	<summary><b>Release v1.1.0 Details</b></summary>
	
* **Modules**:

	* `annotate.py`
		* Added `NCBIfam-AMRFinder` AMR domain annotations
		* Added `AntiFam` contimination annotations
		* Uses `taxopy` instead of `ete3` in backend with `merge_annotations_and_score_taxonomy.py`
	
	* `assembly.py`
		* Added a `transcripts_to_genes.py` script which creates a `genes_to_transcripts.tsv` table that can be used with `TransDecoder`.

	* `binning-prokaryotic.py`
		* Updated `CheckM` → `CheckM2`.  This removes the dependency of `GTDB-Tk` and EXTREMELY REDUCES compute resource requirements (e.g., memory and time) as `CheckM2` automatically handles candidate phyla radiation.  With this, several backend scripts were deprecated.  This cleans up the binning pipeline and error messages SUBSTANTIALLY.
		* Uses `binning_wrapper.py` for all binning.  This makes it easier to add new binning algorithms in the future (e.g., `VAMB`).  Also, check out the new multi-split binning functionality described below.
		* Added `--skip_concoct` in addition to the already existing `--skip_maxbin2` option as `MaxBin2` takes very long when there's a lot of contigs and `CONCOCT` takes a long time when there are a lot of samples (i.e., BAM files).  `MetaBAT2` is not optional.  
		
	* `binning-viral.py`
		* Complete rewrite of this module which now uses `geNomad` as the default binning algorithm but still supports `VirFinder`.
		* If `VirFinder` is used, the `genomad annotate` is run via the `genomad_taxonomy_wrapper.py` script included in the update. 
		* Updated `Prodigal` → `Prodigal-GV` to handle additional viral genetic codes.
		
	* `biosynthetic.py`
		* Introduces `component_id` and `bgc_id` which are unique, pareseable, and informative.  For example, `component_id = SRR17458614__CONCOCT__P.2__9|NODE_3319_length_2682_cov_2.840502|region001_1|2-2681(+)` contains the unique `bgc_id` (i.e., `SRR17458614__CONCOCT__P.2__9|NODE_3319_length_2682_cov_2.840502|region001`), shows that it is the 1st gene in the cluster (the `_1` in `region001_1`), and the gene start/end/strand.  The `bgc_id` is composed of the `genome_id|contig_id|region_id`.
		
	* `classify-prokaryotic.py`
		* Updated `GTDB-Tk v2.1.1` →  `GTDB-Tk v2.2.3`.  For now, `--skip_ani_screen` is the only option because of [this thread](https://forum.gtdb.ecogenomic.org/t/how-can-i-use-the-mash-db-option-of-classify-wf/429/4).  However, `--mash_db` may be an option in the near future.
		* Added functionality to classify prokaryotic genomes that were not binned via `VEBA` which is available with the `--genomes` option (`--prokaryotic_binning_directory` is still available which can leverage existing intermediate files).		
	* `classify-eukaryotic.py`
		* Added functionality to classify eukaryotic genomes that were not binned via `VEBA` which is available with the `--genomes` option (`--eukaryotic_binning_directory` is still available which can leverage existing intermediate files). This is implemented by using the `eukaryota_odb10` markers from the `VEBA Microeukaryotic Database` to substantially improve performance and decrease resources required for gene models.
	
	* `classify-viral.py`
		* Complete rewrite of this module which does not rely on (deprecated) intermediate files from `CheckV`.
		* Uses taxonomy generated from `geNomad` and `consensus_genome_classification_unranked.py` (a wrapper around `taxopy`) that can handle the chaotic taxonomy of viruses.
		* Added functionality to classify viral genomes that were not binned via `VEBA` which is available with the `--genomes` option (`--viral_binning_directory` is still available which can leverage existing intermediate files).

	* `cluster.py`
		* Complete rewrite of this module which now uses `MMSEQS2` as the orthogroup detection algorithm instead of `OrthoFinder`.  `OrthoFinder` is overkill for creating protein clusters and it generates thousands of intermediate files (e.g., fasta, alignments, trees, etc.) which substantially increases the compute time.  `MMSEQS2` has very similar performance with a fraction of the resources and compute time.  Clustered the entire [*Plastisphere* dataset](https://figshare.com/articles/dataset/Genome_assemblies_gene_models_gene_annotations_taxonomy_classifications_clusters_and_counts_tables/20263974) on a local machine in ~30 minutes compared to several days on a HPC.
		* Now that the resources are minimal, clustering is performed at global level as before (i.e., all samples in the dataset) and now at the local level, optionally but ON by default, which clusters all genomes within a sample.  Accompanying wrapper scripts are `global_clustering.py` and `local_clustering.py`.
		* The genomic and functional feature compression ratios (FCR) (described [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04973-8)]) are now calculated automatically.  The calculation is `1 - number_of_clusters/number_of_features` which can easily be converted into an unsupervised biodiversity metric.  This is calculated at the global (original implementation) and local levels.
		* Input is now a table with the following columns: `[organism_type]<tab>[id_sample]<tab>[id_mag]<tab>[genome]<tab>[proteins]` and is generated easily with the `compile_genomes_table.py` script.  This allows clustering to be performed for prokaryotes, eukaryotes, and viruses all at the same time.
		* SLC-specific orthogroups (SSO) are now refered to as SLC-specific protein clusters (SSPC).
		* Support zfilling (e.g., `zfill=3, SLC7 → SLC007`) for genomic and protein clusters.
		* Deprecated `fastani_to_clusters.py` to now use the more generalizable `edgelist_to_clusters.py` which is used for both genomic and protein clusters.  This also outputs a `NetworkX` graph and a pickled dictionary `{"cluster_a":{"component_1", "component_2", ..., "component_n"}}`
		
	* `phylogeny.py`
		* Updated `MUSCLE` to `v5` which has `-align` and `-super5` algorithms which are now accessible with `--alignment_algorithm`.  Cannot use `stdin` so now the fasta files are not gzipped.  The `merge_msa.py` now output uncompressed fasta as default and can output gzipped with the `--gzip` flag.

		
* **`VEBA Database`**:
	* `VDB_v3.1` → `VDB_v4`
		* Updated `CheckV DB v1.0` → `CheckV DB v1.5`
		* Added `geNomad DB v1.2`
		* Added `CheckM2 DB`
		* Removed `CheckM DB`
		* Removed `taxa.sqlite` and `taxa.sqlite.traverse.pkl`
		* Added `reference.eukaryota_odb10.list` and corresponding `MMSEQS2` database (i.e., `microeukaryotic.eukaryota_odb10`)
		* Added `NCBIfam-AMRFinder` marker set for annotation
		* Added `AntiFam` marker set for contamination
		* Marker sets HMMs are now all gzipped (previously could not gzip because `CheckM` CPR workflow)

* **Scripts:**
	* Added:
		* `append_geneid_to_transdecoder_gff.py`
		* `bowtie2_wrapper.py`
		* `compile_genomes_table.py`
		* `consensus_genome_classification_unranked.py`
		* `cut_table.py`
		* `cut_table_by_column_labels.py`
		* `drop_missing_values.py`
		* `edgelist_to_clusters.py`
		* `filter_checkm2_results.py`
		* `genomad_taxonomy_wrapper.py`
		* `global_clustering.py`
		* `local_clustering.py`
		* `partition_multisplit_bins.py`
		* `scaffolds_to_clusters.py`
		* `scaffolds_to_samples.py`
		* `transcripts_to_genes.py`
		* `transdecoder_wrapper.py` (Note: Requires separate environment to run due to dependency conflicts)

	* Updated:
		* `antismash_genbanks_to_table.py` - Added option to output biosynthetic gene cluster (BGC) fasta. Adds unique (and parseable) BGC identifiers making the output much more useful.
		* `binning_wrapper.py` - This binning wrapper now includes functionality to use multi-split binning (i.e., concatenated contigs from different assemblies, map all reads to the contigs, bin all together, and then parition bins by sample).  This concept AFAIK was first introduced in the [`VAMB`](https://www.nature.com/articles/s41587-020-00777-4) paper.
		* `compile_reads_table.py` - Minimal change but now the extension excludes the `.` to make usage more consistent with other tools.
		* `consensus_genome_classification.py` - Changed the output to match that of `consensus_genome_classification_unranked.py`.
		* `filter_checkv_results.py` - Option to use taxonomy and viral summaries generated by `geNomad`.
		* `scaffolds_to_bins.py` - Support for getting scaffolds to bins for a list of genomes via `--genomes` argument while maintaining original support with `--binning_directory` argument.
		* `subset_table.py` - Added option to set index column and to drop duplicates.
		* `virfinder_wrapper.r` - Used to be `VirFinder_wrapper.R`.  This now has an option to use FDR values instead of P values.
		* `merge_annotations_and_score_taxonomy.py` - Completely rewritten.  Uses `taxopy` instead of `ete3`.
		* `merge_msa.py` - Output uncompressed protein fasta files by default and can compress with `--gzip` flag.

	* Deprecated:
		* `adjust_genomes_for_cpr.py`
		* `filter_checkm_results.py`
		* `fastani_to_clusters.py`
		* `partition_orthogroups.py`
		* `partition_clusters.py`
		* `compile_viral_classifications.py`
		* `build_taxa_sqlite.py`

* **Miscellaneous**:

	* Updated environments and now add versions to environments.
	* Added `mamba` to installation to speed up.
	* Added `transdecoder_wrapper.py` which is a wrapper around `TransDecoder` with direct support for `Diamond` and `HMMSearch` homology searches.  Also includes `append_geneid_to_transdecoder_gff.py` which is run in the backend to clean up the GFF file and make them compatible with what is output by `Prodigal` and `MetaEuk` runs of `VEBA`.
	* Added support for using `n_jobs -1` to use all available threads (similar to `scikit-learn` methodology).

</details>

<details>
	<summary><b>Release v1.0.4 Details</b></summary>
	
* Added `biopython` to `VEBA-assembly_env` which is needed when running `MEGAHIT` as the scaffolds are rewritten and [an error](https://github.com/jolespin/veba/issues/17) was raised. [aea51c3](https://github.com/jolespin/veba/commit/aea51c3e0b775aec90f7343f01cad6911f526f0a)
* Updated Microeukaryotic protein database to exclude a few higher eukaryotes that were present in database, changed naming scheme to hash identifiers (from `cat reference.faa | seqkit fx2tab -s -n > id_to_hash.tsv`).  Switching database from [FigShare](https://figshare.com/articles/dataset/Microeukaryotic_Protein_Database/19668855) to [Zenodo](https://zenodo.org/record/7485114#.Y6vZO-zMKDU).  Uses database version `VDB_v3` which has the updated microeukaryotic protein database (`VDB-Microeukaryotic_v2`) [0845ba6](https://github.com/jolespin/veba/commit/0845ba6be65f3486d61fe7ae21a2937efeb42ee9)

</details>

<details>
	<summary><b>Release v1.0.3e Details</b></summary>
	
* Patch fix for `install_veba.sh` where `install/environments/VEBA-assembly_env.yml` raised [a compatibilty error](https://github.com/jolespin/veba/issues/15) when creating the `VEBA-assembly_env` environment. [c2ab957](https://github.com/jolespin/veba/commit/c2ab957be132d34e6b99d6dea394be4572b83066)
* Patch fix for `VirFinder_wrapper.R` where `__version__ = ` variable was throwing [an R error](https://github.com/jolespin/veba/issues/13) when running `binning-viral.py` module. [19e8f38](https://github.com/jolespin/veba/commit/19e8f38a5050328b7ba88b2271f0221073748cbb)
* Patch fix for `filter_busco_results.py` where [an error](https://github.com/jolespin/veba/issues/12) arose that produced empty `identifier_mapping.metaeuk.tsv` subset tables. [359e4569](https://github.com/jolespin/veba/commit/359e45699fc6d6fdf739350263fd34c6e4a62f94)
* Patch fix for `compile_metaeuk_identifiers.py` where [a Python error](https://github.com/jolespin/veba/issues/11) arised when duplicate gene identifiers were present.  [c248527](https://github.com/jolespin/veba/commit/c248527da9edef5ba2ebee348d707d8ece29fbee)
* Patch fix for `install_veba.sh` where `install/environments/VEBA-preprocess_env.yml` raised [a compatibilty error](https://github.com/jolespin/veba/issues/10) when creating the `VEBA-preprocess_env` environment [8ed6eea](https://github.com/jolespin/veba/commit/8ed6eeaee1037694cf324d8fa4da6190578b9688)
* Added `biosynthetic.py` module which runs antiSMASH and converts genbank files to tabular format. [6c0ed82](https://github.com/jolespin/veba/commit/6c0ed82c804ad60a4f1ae51f3e5fecd14dba845f)
* Added `megahit` support for `assembly.py` module (not yet available in `assembly-sequential.py`). [6c0ed82](https://github.com/jolespin/veba/commit/6c0ed82c804ad60a4f1ae51f3e5fecd14dba845f) 
* Changed `-P/--spades_program` to `-P/--program` for `assembly.py`. [6c0ed82](https://github.com/jolespin/veba/commit/6c0ed82c804ad60a4f1ae51f3e5fecd14dba845f)
* Replaced penultimate step in `binning-prokaryotic.py` to use `adjust_genomes_for_cpr.py` instead of the extremely long series of bash commands.  This will make it easier to diagnose errors in this critical step.  [6c0ed82](https://github.com/jolespin/veba/commit/6c0ed82c804ad60a4f1ae51f3e5fecd14dba845f)
* Added support for contig descriptions and added MAG identifier in fasta files in `binning-eukaryotic.py`.  Now uses the `metaeuk_wrapper.py` script for the `MetaEuk` step.  [6c0ed82](https://github.com/jolespin/veba/commit/6c0ed82c804ad60a4f1ae51f3e5fecd14dba845f)
* Added separate option of `--run_metaplasmidspades` for `assembly-sequential.py` instead of making it mandatory (now it just runs `biosyntheticSPAdes` and `metaSPAdes` by default). [6c0ed82](https://github.com/jolespin/veba/commit/6c0ed82c804ad60a4f1ae51f3e5fecd14dba845f)
* Added `--use_mag_as_description` in `parition_gene_models.py` script to include the MAG identifier in the contig description of the fasta header which is default in `binning-prokaryotic.py`. [6c0ed82](https://github.com/jolespin/veba/commit/6c0ed82c804ad60a4f1ae51f3e5fecd14dba845f)
* Added `adjust_genomes_for_cpr.py` script to easier run and understand the CPR adjustment step of `binning-prokaryotic.py`. [6c0ed82](https://github.com/jolespin/veba/commit/6c0ed82c804ad60a4f1ae51f3e5fecd14dba845f)
* Added support for fasta header descriptions in `binning-prokaryotic.py`. [6c0ed82](https://github.com/jolespin/veba/commit/6c0ed82c804ad60a4f1ae51f3e5fecd14dba845f)
* Added functionality to `replace_fasta_descriptions.py` script to be able to use a string for replacing fasta headers in addition to the original functionality. [6c0ed82](https://github.com/jolespin/veba/commit/6c0ed82c804ad60a4f1ae51f3e5fecd14dba845f)

</details>

<details>
	<summary><b>Release v1.0.2a Details</b></summary>

* Updated *GTDB-Tk* in `VEBA-binning-prokaryotic_env` from `1.x` to `2.x` (this version uses much less memory): [f3507dd](https://github.com/jolespin/veba/commit/f3507dd13a42960e3671c9f8a106c9974fbfce21)
* Updated the *GTDB-Tk* database from `R202` to `R207_v2` to be compatible with *GTDB-Tk v2.x*: [f3507dd](https://github.com/jolespin/veba/commit/f3507dd13a42960e3671c9f8a106c9974fbfce21)
* Updated the [GRCh38 no-alt analysis set](https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip) to [T2T CHM13v2.0](https://genome-idx.s3.amazonaws.com/bt/chm13v2.0.zip) for the default human reference: [5ccb4e2](https://github.com/jolespin/veba/commit/5ccb4e20564513707fcc3420b18237974455e196)
* Added an experimental `amplicon.py` module for short-read ASV detection via the *DADA2* workflow of *QIIME2*: [cd4ed2b](https://github.com/jolespin/veba/commit/cd4ed2bfe35d5379a63dd3294c229f2c861f6f77)
* Added additional functionality to `compile_reads_table.py` to handle advanced parsing of samples from fastq directories while also maintaining support for parsing filenames from `veba_output/preprocess`: [cd4ed2b](https://github.com/jolespin/veba/commit/cd4ed2bfe35d5379a63dd3294c229f2c861f6f77)
* Added `sra-tools` to `VEBA-preprocess_env`: [f3507dd](https://github.com/jolespin/veba/commit/f3507dd13a42960e3671c9f8a106c9974fbfce21)
* Fixed symlinks to scripts for `install_veba.sh`: [d1fad03](https://github.com/jolespin/veba/commit/d1fad03b71537cc6cc0d47fee426b6610000752a)
* Added missing `CHECKM_DATA_PATH` environment variable to `VEBA-binning-prokaryotic_env` and `VEBA-classify_env`: [d1fad03](https://github.com/jolespin/veba/commit/d1fad03b71537cc6cc0d47fee426b6610000752a)

</details>



<details>
	<summary><b>Release v1.0.1 Details</b></summary>

* Fixed the fatal binning-eukaryotic.py error: [7c5addf](https://github.com/jolespin/veba/commit/7c5addf9ed6e8e45502274dd353f20b211838a41)
* Fixed the minor file naming in cluster.py: [5803845](https://github.com/jolespin/veba/commit/58038451dac0791899aa7fca3f9d79454cb9ed46)
* Removes left-over human genome tar.gz during database download/config: [5803845](https://github.com/jolespin/veba/commit/58038451dac0791899aa7fca3f9d79454cb9ed46)

</details>


<details>
	<summary><b>Release v1.0.0 Details</b></summary>
	
* Released with *BMC Bionformatics* publication (doi:10.1186/s12859-022-04973-8).

</details>

________________________________________________________________

#### Path to `v3.0.0`:

**Check:**

* Start/end positions on `MetaEuk` gene ID might be off.

**Critical:**

* Don't load all genomes, proteins, and cds into memory for clustering.
* Genome checkpoints in `tRNAscan-SE` aren't working properly.
* Dereplicate CDS sequences in GFF from `MetaEuk` for `antiSMASH` to work for eukaryotic genomes
* Error with `amplicon.py` that works when run manually... (Developmental module)

```
There was a problem importing veba_output/misc/reads_table.tsv:

  Missing one or more files for SingleLanePerSamplePairedEndFastqDirFmt: '.+_.+_L[0-9][0-9][0-9]_R[12]_001\\.fastq\\.gz'
```

**Definitely:**

* Add `BiNI` biosynthetic novelty index to `biosynthetic.py`
* `busco_wrapper.py` that relabels all the genes, runs analysis, then converts output to tsv.
* Script to update genome clusters
* Script to update protein clusters
* Script to add `Diamond` or `HMMSearch` annotations to `annotations.proteins.tsv.gz`
* Add `convert_reads_long_to_short.py` which will take windows of 150 bp for the long reads.
* Add option to `compile_custom_humann_database_from_annotations.py` to only output best hit of a UniRef identifier per genome.
* Use `pigz` instead of `gzip`
* Create a taxdump for `MicroEuk`
* Reimplement `compile_eukaryotic_classifications.py`
* Add representative to `identifier_mapping.proteins.tsv.gz`
* Use `aria2` in parallel instead of `wget`.
* Add support for `Salmon` in `mapping.py` and `index.py`.  This can be used instead of `STAR` which will require adding the `exon` field to `Prodigal` GFF file (`MetaEuk` modified GFF files already have exon ids). 
* [Optional] Number of plasmids (via `geNomad`) for each MAG.


**Eventually (Yes)?:**

* Support for `MAFFT` instead of `MUSCLE` as it performs especially well for multidomain protein with variable domain architectures.
* `NextFlow` support
* Install each module via `bioconda`
* Consistent usage of the following terms: 1) dataframe vs. table; 2) protein-cluster vs. orthogroup.  Dataframes should refer to generic tables while tables refer to specifics like "genomes table".
* Add coding density to GFF files
* Run `cmsearch` before `tRNAscan-SE`
* DN/DS from pangeome analysis
* Add a `metabolic.py` module	
* For viral binning, contigs that are not identified as viral via `geNomad -> CheckV` use with `vRhyme`.
* Add `vRhyme` to `binning_wrapper.py` and support `vRhyme` in `binning-viral.py`.


**...Maybe (Not)?**

* Swap [`TransDecoder`](https://github.com/TransDecoder/TransDecoder) for [`TransSuite`](https://github.com/anonconda/TranSuite)

________________________________________________________________


<details>
	<summary> <b>Daily Change Log:</b> </summary>

* [2024.6.7] - Adapted `phylogeny.py` and `partition_pyhmmsearch.py` to use `pyhmmsearch` instead of `hmmsearch` and `Kofam_Scan`.
* [2024.6.7] - Adapted `annotate.py`, `merge_annotations.py`, and `compile_ko_from_annotations.py` to use `pyhmmsearch` and `pykofamsearch` instead of `hmmsearch` and `Kofam_Scan`.
* [2024.6.6] - Changed `Diamond` output format from `-f 6 qseqid sseqid stitle pident length mismatch qlen qstart qend slen sstart send evalue bitscore qcovhsp scovhsp` to `-f 6 qseqid sseqid stitle pident evalue bitscore qcovhsp scovhsp`
* [2024.6.6] - Adapted `classify-eukaryotic.py` to use `pyhmmsearch` instead of `hmmsearch`.
* [2024.6.6] - Updating `GTDB-Tk` and `BUSCO` introduced conflicting dependencies.  To provide more flexibility for version updates, `VEBA-classify_env` has been split out into `VEBA-classify-eukaryotic_env`, `VEBA-classify-prokaryotic_env`, and `VEBA-classify-viral_env`.
* [2024.6.5] - Update `GTDB` version from `r214.1` to `r220` in VEBA database version `VDB_v7` and in `classify-prokaryotic.py`.  Corresponding mash database for `r220` is available here: 
* [2024.6.5] - Added `choose_fastest_mirror.py` to utility scripts which checks the speed of multiple urls and then outputs the fastest one.
* [2024.6.5] - Removing version name from `GTDB` .msh file.  Previous versions included `gtdb_r214.msh` but now they will be `gtdb.sh`.
* [2024.5.20] - Added `reformat_minpath_report.py` to reformat minpath reports.  `MinPath` isn't used directly by VEBA but it might be in the future.
* [2024.4.30] - Added `concatenate_files.py` which can concatenate files (and mixed compressed/decompressed files) using either arguments, list file, or glob.  Reason for this is that unix has a limit of arguments that can be used (e.g., `cat *.fasta > output.fasta` where *.fasta results in 50k files will crash)
* [2024.4.29] - Added `/volumes/workspace/` directory to Docker containers for situations when your input and output directories are the same. 
* [2024.4.29] - `featureCounts` can only handle 64 threads at a time so added `min(64, opts.n_jobs)` for all the modules/scripts that use `featureCounts` commands.
* [2024.4.23] - Added `uniprot_to_enzymes.py` which reformats tables and fasta from https://www.uniprot.org/uniprotkb?query=ec%3A* 
* [2024.4.18] - Developed a faster implementation of `KofamScan` called [`PyKofamSearch`](https://github.com/jolespin/pykofamsearch) which leverage `PyHmmer`.  This will be used in future versions of VEBA.
* [2024.3.26] - Added `--metaeuk_split_memory_limit` to `metaeuk_wrapper.py`.
* [2024.3.26] - Added `-d/--genome_identifier_directory_index` to `scaffolds_to_bins.py` for directories that are structured `path/to/genomes/bin_a/reference.fasta` where you would use `-d -2`.
* [2024.3.26] - Added `--minimum_af` to `edgelist_to_clusters.py` with an option to accept 4 column inputs `[id_1]<tab>[id_2]<tab>[weight]<tab>[alignment_fraction]`.  `global_clustering.py`, `local_clustering.py`, and `cluster.py` now use this by default `--af_threshold 30.0`.  If you want to retain previous behavior, just use `--af_threshold 0.0`.
* [2024.3.18] - `edgelist_to_clusters.py` only includes edges where both nodes are in identifiers set.  If `--identifiers` are provided, then only those identifiers are used.  If not, then it includes all nodes.
* [2024.3.18] - Added `--export_representatives` argument for `edgelist_to_clusters.py` to output table with `[id_node]<tab>[id_cluster]<tab>[intra-cluster_connectivity]<tab>[representative]`.  Also includes this information in `nx.Graph` objects.
* [2024.3.18] - Changed singleton weight to `np.nan` instead of `np.inf` for `edgelist_to_clusters.py` to allow for representative calculations. 
* [2024.3.8] - Changed default assembly algorithm to `metaflye` instead of `flye` in `assembly-long.py`
* [2024.3.8] - Added `number_of_genomes`, `number_of_genome-clusters`, `number_of_proteins`, and `number_of_protein-clusters` to `feature_compression_ratios.tsv.gz` from `cluster.py`
* [2024.3.5] - Added `-A/--from_antismash` in `biosynthetic.py` to use preexisting `antiSMASH` results.  Also changed `-i/--input` to `-i/--from_genomes`.
* [2024.3.4] - Changed `antimash_genbanks_to_table.py` to `biosynthetic_genbanks_to_table.py` for future support of `DeepBGC` and `GECCO`
* [2024.2.28] - Added `busco_version` parameter to `merge_busco_json.py` with default set to `5.4.x` and additional support for `5.6.x`.
* [2024.2.24] - Added `CONDA_ENVS_PATH` to `update_environment_scripts.sh`, `update_environment_variables.sh`, and `check_installation.sh`
* [2024.2.17] - Added `CONDA_ENVS_PATH` to `veba` to allow for custom environment locations
* [2024.2.16] - Changed `install.sh` to support custom `CONDA_ENVS_PATH` argument `bash install.sh path/to/log path/to/envs/`
* [2024.2.16] - Added `merge_counts_with_taxonomy.py`
* [2024.1.28] - Replaced `src/` with `bin/` and added `-V|--full_versions to show all VEBA versions`
* [2024.1.23] - Added `compile_phylogenomic_functional_categories.py` script which automates the methodology from [Espinoza et al. 2022 (doi:10.1093/pnasnexus/pgac239)](https://academic.oup.com/pnasnexus/article/1/5/pgac239/6762943)
* [2024.1.22] - Fixed header being offset in `annotations.protein_clusters.tsv` where it could not be read with Pandas.
* [2024.1.22] - Fixed `binning-prokaryotic.py` the creation of non-existing symlinks where "'*.gff'", "'*.rRNA'", and "'*.tRNA'" were created.
* [2024.1.16] - Fixed .strip method on Pandas series in `antismash_genbanks_to_table.py` for compatibilty with `antiSMASH 6 and 7`
* [2024.1.7] - Fixed situation where `unbinned.fasta` is empty in `binning-prokaryotic.py` when there are no bins that pass qc.
* [2024.1.7] - Fixed minor error in `coverage.py` where `samtools sort --reference` was getting `reads_table.tsv` and not `reference.fasta`
* [2023.1.4] - Changed default behavior from deterministic to not deterministic for increase in speed in `assembly-long.py`. (i.e., `--no_deterministic` --> `--deterministic`)
* [2024.1.2] - Added `VeryFastTree` as an option to `phylogeny.py` with `FastTree` remaining as the default.
* [2023.12.30] - Changed default `--leniency` parameter on `classify_eukaryotic.py` and `consensus_genome_classification_ranked.py` to `1.0` and added `--leniecy_genome_classification` as a separate option.
* [2023.12.28] - Added `--blacklist` option to `compile_eukaryotic_classifications.py` with a default value of `species:uncultured eukaryote` in `classify_eukaryotic.py`
* [2023.12.28] - Fixed critical error where `classify_eukaryotic.py` was trying to access a deprecated database file from MicrEuk_v2.
* [2023.12.22] - Fixed minor error with `eukaryotic_gene_modeling_wrapper.py` not allowing for `Tiara` to run in backend.
* [2023.12.21] - `GTDB-Tk` changed name of archaea summary file so VEBA was not adding this to final classification. Fixed this in `classify-prokaryotic.py`.
* [2023.12.20] - Fixed files not being closed in `compile_custom_humann_database_from_annotations.py` and added options to use different annotation file formats (i.e., multilevel, header, and no header).
* [2023.12.15] - Added `profile-taxonomic.py` module which uses `sylph` to build a sketch database for genomes and queries the genome database similar to `Kraken` for taxonomic abundance.
* [2023.12.14] - Removed requirement to have `--estimated_assembly_size` for Flye per [Flye Issue #652](https://github.com/fenderglass/Flye/issues/652).
* [2023.12.14] - Added `sylph` to `VEBA-profile_env` for abundance profiling of genomes.
* [2023.12.13] - Dereplicate duplicate contigs in `concatenate_fasta.py`.
* [2023.12.12] - Added `--reference_gzipped` to `index.py` and `mapping.py` with new default being that the reference fasta is not gzipped.
* [2023.12.11] - Added `skani` as new default for genome clustering in `cluster.py`, `global_clustering.py`, and `local_clustering.py`.
* [2023.12.11] - Added support for long reads in `fastq_preprocessor`, `preprocess.py`, `assembly-long.py`, and all binning modules.
* [2023.11.28] - Fixed `annotations.protein_clusters.tsv.gz` from `merge_annotations.py` added in patch update of `v1.3.1`.
* [2023.11.14] - Added support for missing values in `compile_eukaryotic_classifications.py`.
* [2023.11.13] - Added `--metaeuk_split_memory_limit` argument with (experimental) default set to `36G` in `binning-eukaryotic.py` and `eukaryotic_gene_modeling.py`.
* [2023.11.10] - Added `--compressed 1` to `mmseqs createdb` in `download_databases.sh` installation script.
* [2023.11.10] - Added a check to `check_fasta_duplicates.py` and `clean_fasta.py` to make sure there are no `>` characters in fasta sequence caused from concatenating fasta files that are missing linebreaks.
* [2023.11.10] - Added `Diamond DeepClust` to `clustering_wrapper.py`, `global/local_clustering.py`, and `cluster.py`.  Changed `mmseqs2_wrapper.py` to `clustering_wrapper.py`.  Changed `easy-cluster` and `easy-linclust` to `mmseqs-cluster` and `mmseqs-linclust`.
* [2023.11.9] - Fixed viral quality in `merge_genome_quality_assessments.py`
* [2023.11.3] - Changed `consensus_genome_classification.py` to `consensus_genome_classification_ranked.py`.  Also, default behavior to allow for missing taxonomic levels.
* [2023.11.2] - Fixed the `merge_annotations.py` resulting in a memory leak when creating the `annotations.protein_clusters.tsv.gz` output table.  However, still need to correct the formatting for empty sets and string lists.
* [2023.10.27] - Update `annotate.py` and `merge_annotations.py` to handle `CAZy`.  They also properly address clustered protein annotations now. 
* [2023.10.18] - Added `module_completion_ratio.py` script which is a fork of `MicrobeAnnotator` [`ko_mapper.py`](https://github.com/cruizperez/MicrobeAnnotator/blob/master/microbeannotator/pipeline/ko_mapper.py).  Also included a database [Zenodo: 10020074](https://zenodo.org/records/10020074) which will be included in `VDB_v5.2`
* [2023.10.16] - Added a checkpoint for `tRNAscan-SE` in `binning-prokaryotic.py` and `eukaryotic_gene_modeling_wrapper.py`.
* [2023.10.16] - Added `profile-pathway.py` module and `VEBA-profile_env` environments which is a wrapper around `HUMAnN` for the custom database created from `annotate.py` and `compile_custom_humann_database_from_annotations.py`
* [2023.10.16] - Added `GenoPype version` to log output
* [2023.10.16] - Added `merge_genome_quality.py` which combines `CheckV`, `CheckM2`, and `BUSCO` results.
* [2023.10.11] - Added `compile_custom_humann_database_from_annotations.py` which compiles a `HUMAnN` protein database table from the output of `annotate.py` and taxonomy classifications.
* [2023.10.11] - Added functionality to `merge_taxonomy_classifications.py` to allow for `--no_domain` and `--no_header` which will serve as input to `compile_custom_humann_database_from_annotations.py`
* [2023.10.5] - Added `marker_gene_clustering.py` script which gets core marker genes unique to each SLC (i.e., pangenome). `average_number_of_copies_per_genome` to protein clusters.
* [2023.10.5] - Added `--minimum_core_prevalence` in `global_clustering.py`, `local_clustering.py`, and `cluster.py` which indicates prevalence ratio of protein clusters in a SLC will be considered core.  Also remove `--no_singletons` from `cluster.py` to avoid complications with marker genes.  Relabeled `--input` to `--genomes_table` in clustering scripts/module.
* [2023.9.21] - Added a check in `coverage.py` to see if the `mapped.sorted.bam` files are created, if they are then skip them.  Not yet implemented for GNU parallel option.
* [2023.9.15] - Changed default representative sequence format from table to fasta for `mmseqs2_wrapper.py`.
* [2023.9.12] - Added `--nucleotide_fasta_output` to `antismash_genbank_to_table.py` which outputs the actual BGC DNA sequence. Changed `--fasta_output` to `--protein_fasta_output` and added output to `biosynthetic.py`. Changed BGC component identifiers to `[bgc_id]_[position_in_bgc]|[start]:[end]([strand])` to match with `MetaEuk` identifiers. Changed `bgc_type` to `protocluster_type`.  `biosynthetic.py` now supports GFF files from `MetaEuk` (exon and gene features not supported by `antiSMASH`).  Fixed error related to `antiSMASH` adding CDS (i.e., `allorf_[start]_[end]`) that are not in GFF so `antismash_genbank_to_table.py` failed in those cases. 
* [2023.9.12] - Added `ete3` to `VEBA-phylogeny_env.yml` and automatically renders trees to PDF. #! Need to test
* [2023.9.11] - Added presets for `MEGAHIT` using the `--megahit_preset` option. 
* [2023.9.11] - The change for using `--mash_db` with `GTDB-Tk` violated the assumption that all prokaryotic classifications had a `msa_percent` field which caused the cluster-level taxonomy to fail.  `compile_prokaryotic_genome_cluster_classification_scores_table.py` fixes this by uses `fastani_ani` as the weight when genomes were classified using ANI and `msa_percent` for everything else.  Initial error caused unclassified prokaryotic for all cluster-level classifications.
* [2023.9.8] - Fixed small error where empty gff files with an asterisk in the name were created for samples that didn't have any prokaryotic MAGs.
* [2023.9.8] - Fixed critical error where descriptions in header were not being removed in `eukaryota.scaffolds.list` and did not remove eukaryotic scaffolds in `seqkit grep` so `DAS_Tool` output eukaryotic MAGs in `identifier_mapping.tsv` and `__DASTool_scaffolds2bin.no_eukaryota.txt`
* [2023.9.5] - Fixed `krona.html` in `biosynthetic.py` which was being created incorrectly from `compile_krona.py` script.
* [2023.8.30] - Create `pangenome_core_sequences` in `global_clustering.py` and `local_clustering.py` which writes both protein and CDS sequences for each SLC.  Also made default in `cluster.py` to NOT do local clustering switching `--no_local_clustering` to `--local_clustering`.
* [2023.8.30] - `pandas.errors.InvalidIndexError: Reindexing only valid with uniquely valued Index objects` in `biosynthetic.py` when `Diamond` finds multiple regions in one hit that matches.  Added `--sort_by` and `--ascending` to `concatenate_dataframes.py` along with automatic detection and removal of duplicate indices.  Also added `--sort_by bitscore` in `biosynthetic.py`.
* [2023.8.28] - Added core pangenome and singleton hits to clustering output
* [2023.8.25] - Updated `--megahit_memory` default from 0.9 to 0.99
* [2023.8.16] - Fixed error in `genomad_taxonomy_wrapper.py` where `viral_taxonomy.tsv` should have been `taxonomy.tsv`.	
* [2023.7.26] - Fixed minor error in `assembly.py` that was preventing users from using `SPAdes` programs that were not `spades.py`, `metaspades.py`, or `rnaspades.py` that was the result of using an incorrect string formatting.
* [2023.7.25] - Updated `bowtie2` in preprocess, assembly, and mapping modules.  Updated `fastp` and `fastq_preprocessor` in preprocess module.
* [2023.7.7] - Added `compile_gff.py` to merge CDS, rRNA, and tRNA GFF files.  Used in `binning-prokaryotic.py` and `binning-viral.py`.  `binning-eukaryotic.py` uses the source of this in the backend of `filter_busco_results.py`. Includes GC content for contigs and various tags. 
* [2023.7.6] - Updated `BUSCO v5.3.2 -> v5.4.3` which changes the json output structure and made the appropriate changes in `filter_busco_results.py`.
* [2023.7.3] - Added `eukaryotic_gene_modeling_wrapper.py` which 1) splits nuclear, mitochondrial, and plastid genomes; 2) performs gene modeling via `MetaEuk` and `Pyrodigal`; 3) performs rRNA detection via `BARRNAP`; 4) performs tRNA detection via `tRNAscan-SE`; 5) merges processed GFF files; and 5) calculates sequences statistics. 
* [2023.6.29] - Added `gene_biotype=protein_coding` to `prodigal` GFF output. 
* [2023.6.20] - Added `VFDB` to `annotate.py` and database.
* [2023.6.16] - Compiled and pushed `gtdb_r214.msh` mash file to [Zenodo:8048187](https://zenodo.org/record/8048187) which is now used by default in `classify-prokaryotic.py`.  It is now included in `VDB_v5.1`.
* [2023.6.15] - Cleaned up global and local clustering intermediate files.  Added pangenome tables and singelton information to outputs.
* [2023.6.12] - Changed `${VEBA_DATABASE}/Classify/GTDBTk` → `${VEBA_DATABASE}/Classify/GTDB`.
* [2023.6.12] - Replace `prodigal` with `pyrodigal` in `binning-prokaryotic.py` (`prodigal` is still in environment b/c `DAS_Tool` dependency).
* [2023.6.12] - `consensus_genome_classification.py` now based missing classifications off of a missing weight value. Defaults for unclassified labels are `Unclassified prokaryote`, `Unclassified eukaryote`, and `Unclassified virus` for the various classification modules. Also changed "id_genome_cluster" to "id" and "genomes" to "components" to generalize for eukaryotic classification.
* [2023.6.12] - `global_clustering.py` and `local_clustering.py` (accessible through `cluster.py`) now outputs NetworkX graph and Python dictionary pickled objects.
* [2023.6.12] - Added support for missing values and unclassified taxa in `compile_krona.py` and `consensus_genome_classification.py`.  
* [2023.5.18] - Added `compile_protein_cluster_prevalence_table.py` script
* [2023.5.17] - Added `convert_table_to_fasta.py` script
* [2023.5.16] - Created Docker images for all modules
* [2023.5.16] - Replaced all absolute path symlinks with relative symlinks.
* [2023.5.15] - Changed `prokaryotic_taxonomy.tsv` and `prokaryotic_taxonomy.clusters.tsv` in `classify-prokaryotic.py` (along with eukaryotic and viral) files to `taxonomy.tsv` and `taxonomy.clusters.tsv` for uniformity.
* [2023.5.15] - Updating all symlinks to relative links (also in `fastq_preprocessor`) to prepare for dockerization and updating all environments to use updated GenoPype 2023.4.13.
* [2023.5.14] - Changed `nr` to `uniref` in `annotate.py` and added `propagate_annotations_from_representatives.py` script while simplifying `merge_annotations_and_taxonomy.py` to `merge_annotations.py` and excluding taxonomy operations.
* [2023.5.14] - Changed `nr` to `UniRef90` and `UniRef50` in `VDB_v5` 
* [2023.5.12] - Changed `orfs_to_orthogroups.tsv` to `proteins_to_orthogroups.tsv` for consistency with the `cluster.py` module.  Will eventually find some consitency with `scaffolds_to_bins/scaffolds_to_mags` but this will be later.
* [2023.5.12] - Added a `scaffolds_to_mags.tsv` in the clustering output.
* [2023.5.8] - Added `convert_counts_table.py` which converts a counts table (and metadata) to [Pandas pickle](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_pickle.html), [Anndata h5ad](https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.write_h5ad.html), or [Biom hdf5](https://biom-format.org/documentation/generated/biom.table.Table.to_hdf5.html#biom.table.Table.to_hdf5)
* [2023.5.8] - Fixed output directory for `mapping.py` which now uses `output_directory/${NAME}` structure like `binning-*.py`.
* [2023.5.8] - Removed "python" prefix for script calls and now uses shebang in script for executable. Also added single paranthesis around script filepath (e.g., `'[script_filepath]'`) to escape characters/spaces in filepath.
* [2023.5.8] - Added support for `index.py` to accept individual `--references [file.fasta]` and `--gene_models [file.gff]`.
* [2023.4.25] - Added `stdin` support for `scaffolds_to_bins.py` along with the ability to input genome tables [id_genome]<tab>[filepath].  Also added progress bars.
* [2023.4.23] - As a result of [issues/22](https://github.com/jolespin/veba/issues/22), `assembly.py`, `assembly-sequential.py`, `binning-*.py`, and `mapping.py` will use `-p --countReadPairs` for `featureCounts` and updates `subread 2.0.1 → subread 2.0.3`.  For `binning-*.py`, long reads can be used with the `--long_reads` flag.
* [2023.4.20] - Updated `cluster.py` and associated `global_clustering.py`/`local_clustering.py` scripts to use `mmseqs2_wrapper.py` which now automatically outputs representative sequences.  
* [2023.4.17] - Added `check_fasta_duplicates.py` script that gives `0` and `1` exit codes for fasta without and with duplicates, respectively.  Added `reformat_representative_sequences.py` to reformat representative sequences from `MMSEQS2` into either a table or fasta file where the identifers are cluster labels.  Removed `--dbtype` from `[global/local]_clustering.py`.  Removed appended prefix for `.graph.pkl` and `dict.pkl` in `edgelist_to_clusters.py`.  Added `mmseqs2_wrapper.py` and `hmmer_wrapper.py` scripts.
* [2023.4.13] - Added an option to `merge_generalized_mapping.py` to include the sample index in a filepath and also an option to remove empty features (useful for Salmon).  Added an `executable='/bin/bash'` option to the `subprocess.Popen` calls in `GenoPype` to address [issues/23](https://github.com/jolespin/veba/issues/23).
* [2023.3.23] - Added `genbanks/[id_genome]/` to output directory of `biosynthetic.py` which has symlinks to all the BGC genbanks from `antiSMASH`.
* [2023.3.20] - Added `database` field to `source_taxonomy.tsv.gz` in `VDB-Microeukaryotic_v2.1` as [an additional file](https://zenodo.org/record/7485114/files/source_taxonomy_with_database.tsv.gz?download=1) which wille eventually replace the default file.  Also changed `SourceID` to `id_source` in updated version.
* [2023.3.17] - Fixed rare bug when `antiSMASH` genbank files have a space appended to the contig.  Also fixed a typo in the BGC features fasta file name.
* [2023.3.13] - Fixed `--skip_maxbin2` and `--skip_concoct` arguments by adding missing `seed` parameters ([Issue #21](https://github.com/jolespin/veba/issues/21)).  Added a wrapper around `STAR` RNAseq-aligner (`star_wrapper.py`) in preperation to add as an option for `mapping.py`.  This also includes a helper script in compiling the summary log (`compile_star_statistics.py`).
* [2023.3.9] - Added `bgc_novelty_scorer.py` script to get novelty scores of biosynthetic gene clusters.
* [2023.3.7] - Added prefix and minimum contig length threshold to `assembly.py` by default.  Added `merge_generalized_mapping.py` which can be used for `bowtie2_wrapper.py` and (the future) `star_wrapper.py` helper scripts.
* [2023.3.6] - Added dereplicated `MIBiG` `Diamond` database to (`mibig_v3.1.dmnd`) `VDB_v4.1`.  Adds protein fasta files for genes in BGCs for `biosynthetic.py` which are used to run against the `mibig_v3.1.dmnd` database.  
* [2023.3.3] - Updated `binning-viral.py` module's `geNomad` run to use `--relaxed` settings by default since `CheckV` is used after with conservative settings (https://portal.nersc.gov/genomad/post_classification_filtering.html#default-parameters-and-presets)
* [2023.2.23] - The largest update to date.  Please refer to v1.1 for details on what has been changed.
* [2023.01.20] - Changed `-a --ani` to `-t --threshold` in `fastani_to_clusters.py` to match the usage in `edgelist_to_clusters.py` which is a generalization of `fastani_to_clusters.py` developed for `MMSEQS2` and `Diamond` implementations.
* [2023.01.12] - Updated `VDB-Microeukaryotic_v2` to `VDB-Microeukaryotic_v2.1` to include a `reference.eukaryota_odb10.list` containing all the eukaryotic core markers. To accomodate this, I've also updated `VDB_v3` to `VDB_v3.1`, the `download_databases.sh` script, and the `VEBA-database_env.yml` environment file.  Now a `microeukaryotic.eukaryota_odb10` will be available for streamlined eukaryotic classification.
* [2023.01.11] - `biosynthetic.py` automatically removes assembly.gbk and assembly.json files because they are big and unnecessary.
* [2023.01.08] - Added an internal checkpoint system for `biosynthetic.py` when re-running an incomplete `antiSMASH` step (useful when running large numbers of genomes).  Fixed the follow environment files: `VEBA-amplicon_env.yml`, `VEBA-binning-prokaryotic_env.yml`, and `VEBA-binning-eukaryotic_env.yml` as they had either PyPI or package conflict errors during installation.
* [2023.01.05] - Added start, end, and strand to antismash output table in `antismash_genbanks_to_table.py`. Output is sorted by `["genome_id", "contig_id", "start", "end"]` Fixed `VEBA-phylogney_env.yml` environment file.  Important fix in `update_environment_scripts.sh` for symlinking scripts in path. 
* [2023.01.03] - Moving `VEBA-biosynthetic_env` as a developmental environment so it won't be installed automatically.  The reasoning for this is that `antiSMASH` downloads and configures that `antiSMASH database` in the backend which uses a lot of compute resources and takes a long time.  Didn't want to slow up the installation more. 
* [2022.12.21] - Added `biopython` to `VEBA-assembly_env` which is needed when running `MEGAHIT` as the scaffolds are rewritten.
* [2022.12.12] - Fixed duplicate `step__step__program` labels for `classify-prokaryotic.py` module.  Added support for prepending index/column levels and `index_col` selection in `concatenate_dataframes.py`.
* [2022.12.07] - Fixed the compatibility issues for `VEBA-preprocess_env.yml` and issues with the following scripts: `compile_metaeuk_identifiers.py`, `filter_busco_results.py`, and `VirFinder_wrapper.R`.
* [2022.11.14] - Added `megahit` support for `assembly.py` module (not yet available in `assembly-sequential.py`).  Changed `-P/--spades_program` to `-P/--program` for `assembly.py`. Added `biosynthetic` module which runs antiSMASH and converts genbank files to tabular format.  `binning-prokaryotic.py` defaults to `TMPDIR` environment variable for CheckM step, if not available, then it uses `[PROJECT_DIRECTORY]/[ID]/tmp`.  See #12 of [FAQ](https://github.com/jolespin/veba/blob/main/FAQ.md). 
* [2022.11.8] - Replaced penultimate step in `binning-prokaryotic.py` to use `adjust_genomes_for_cpr.py` instead of the extremely long series of bash commands.  This will make it easier to diagnose errors in this critical step.  Also added support for contig descriptions and added MAG identifier in fasta files in `binning-eukaryotic.py`.  Now uses the `metaeuk_wrapper.py` script for the `MetaEuk` step.  Added separate option of `--run_metaplasmidspades` for `assembly-sequential.py` instead of making it mandatory (now it just runs `biosyntheticSPAdes` and `metaSPAdes` by default).
* [2022.11.7] - Added `--use_mag_as_description` in `parition_gene_models.py` script to include the MAG identifier in the contig description of the fasta header which is default in `binning-prokaryotic.py`. Added `adjust_genomes_for_cpr.py` script to easier run and understand the CPR adjustment step of `binning-prokaryotic.py`. Added support for fasta header descriptions in `binning-prokaryotic.py`.
* [2022.11.4] - Added functionality to `replace_fasta_descriptions.py` script to be able to use a string for replacing fasta headers in addition to the original functionality.
* [2022.10.26] - Fixed symlinks to scripts for `install_veba.sh` and added missing `CHECKM_DATA_PATH` environment variable.  Also added `uninstall_veba.sh`, added `update_environment_variables.sh` scripts, and cleaned up install/database scripts.
* [2022.10.25] - Updated default `GTDB-Tk` database from `R202` to `R207_v2` and along with this updated `GTDB-Tk` in `VEBA-binning-prokaryotic_env` and `VEBA-classify_env`.  Also, updated the `binning-prokaryotic.py` to include the `checkm_output.filtered.tsv` instead of unfiltered `output.tsv`.
* [2022.10.24] - Added new functionality to `compile_reads_table.py` by adding a method to compile reads tables from Fastq directories.  Compatible with `QIIME2` manifest. Defaults to absolute path with added option `--relative` for relative paths. Also added an experimental `amplicon.py` module for ASV detection/classification along with the appropriate environment recipe and README.md update.
* [2022.10.18] - Replace `GRCh38 alt analysis set` with the `CHM13v2.0 telomere-to-telomere build` for the included human reference.  Also updated the `VEBA-database_env` to include `unzip` and added a patch for users to update their human reference if desired.
* [2022.10.16] - Added `edgelist.tsv` and `graph.pkl` to output directory for cluster.py.  These files were already in intermediate 1__fastani but the file name was weird. (e.g., graph.pkl-ani_95.0.edgelist.tsv).  Fixed and also changed output of graph in `fastani_to_clusters.py`:nx.write_gpickle(graph, "{}-ani_{}.graph.pkl".format(opts.export_pickle, tol)) (Adding the .graph. part).
* [2022.08.27] - Added `metaeuk_wrapper.py` script
* [2022.08.17] - Added --scaffolds_to_bins option to mapping.py. Automated `samtools index` and `samtools coverage` steps for `mapped.sorted.bam` files for use spatial coverage calculation which is now automated as well producing a `genome_spatial_coverage.tsv.gz` file if `--scaffolds_to_bins` is provided.  Added `genome_spatial_coverage.py` script that uses the `samtools coverage` files from the `mapping.py` module or custom runs. Fixed --table_header error in `groupby_table.py` script.
* [2022.07.14] - Added `--absolute` argument to `compile_reads_table.py` to use absolute paths instead of relative paths.
* [2022.07.14] - Added `genome_coverage_from_spades.py` script to help with coverage calculations for NCBI submissions.
* [2022.07.13] - Added output files to documentation readme.  Also added subread to viral binning environment and changed the output filenames for viral classification.
* [2022.07.08] - Fixed database arguments to use --veba_database. Also added --bam support for binning-viral.py
* [2022.06.21] - Added `prefiltered_alignment_table.tsv.gz` to `phylogeny.py` and `merge_msa.py` for debugging and finding a balance between removing genomes and removing markers. Changed `samples` and `number_of_samples` to `genomes` and `number_of_genomes` in `merge_msa.py`. Also added `minimum_markers_aligned_ratio` to remove poor quality genomes.
* [2022.06.20] - Added `filter_hmmsearch_results.py` and score thresholding table input for `phylogeny.py` and `partition_hmmsearch.py`.  Changed `minimum_genomes_aligned_ratio` default to 0.95 instead of 0.5.
* [2022.06.03] - Changed `coassembly.py` to `coverage.py` which includes `veba_output/assembly/coassembly` to `veba_output/assembly/multisample` and `coassembly.fasta.*` to `reference.fasta.*`. Also change the GNU parallel default to an optional since it's much slower when the samples are different sizes. This can be selected using --one_task_per_cpu argument. Should probably do this for `cluster.py`. 
* [2022.05.25] - Added --skip_maxbin2 argument for binning-prokaryotic. The reason for this is that it takes an extremely long time. For a 1.5GB fasta file and ~50 or so samples in the coverage matrix it takes over 40 hours per MaxBin2 run (2 per run). This will be 30 days to run 10 iterations. 
* [2022.04.12] - Added coassembly module and support for multiple bam files in binning-prokaryotic, binning-eukaryotic, and binning-wrapper
* [2022.03.28] - Added GTDB-Tk to prokaryotic binning so check for CPR and then rerun CheckM using the proper parameters.
* [2022.03.14] - Created a `binning_wrapper.py` to normalize the binning process and add --minimum_genome_length capabilities. This is useful for eukaryotic binning but more complicated for prokaryotic binning because the current pipeline is hardcoded to handle errors on iterative binning. Also switched to CoverM for all coverage calculations because it's faster. Split out prokaryotic, eukaryotic, and viral binning environments. For eukaryotic binning, I've removed EukCC and use BUSCO v5 instead.
* [2022.03.01] - Added a domain classification script that is run during prokaryotic binning. I've created a hack that moves all of the eukaryotic genomes to another directory to allow for proper gene calls in a separate module. This hack will remain until DAS_Tool can handle custom gene sets because it cannot in the current version. The other option is to remove 
* [2022.02.24] - Added saf file to `assembly.py` and feature counts of scaffolds/transcripts
* [2022.02.22] - Made the original `preprocess.py` → `preprocess-kneaddata.py` and the new `preprocess.py` a wrapper around `fastq_preprocessor`
* [2022.02.22] - Made the `index.py` module
* [2022.02.22] - `concatenate_fasta.py` and `concatenate_gff.py`
* [2022.02.02] - `consensus_genome_classification.py`

</details>









