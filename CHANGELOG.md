### What's next for *VEBA*?

*VEBA* is currently under active development. If you are interested in requesting features or wish to report a bug, please post a GitHub issue prefixed with the tag `[Feature Request]` and `[Bug]`, respectively.  If you want to contribute or have any other inquiries, contact me at `jol.espinoz[A|T]gmail[DOT]com`.

________________________________________________________________

#### Current Releases:

# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.5.1] - 2025.04.11

### Added

*   Added `install-gpu.sh` which installs GPU accelerated environments when applicable (i.e., `VEBA-binning-prokaryotic_env` and `VEBA-binning-viral_env`)
*   Added `Dockerfile-GPU` which is experimental

### Changed

*   Changed `install.sh` so it only installs CPU-based environments [Issue #167](https://github.com/jolespin/veba/issues/167)
*   Changed `containerize_environments.sh` so it only installs CPU-based environments [Issue #167](https://github.com/jolespin/veba/issues/167)

### Deprecated

*   Deprecated `VirFinder` algorithm in `binning-viral.py` so now only `geNomad` is supported

## [2.5.0] - 2025.04.10

### Added

*   Added `VAMB` support to `binning-prokaryotic.py` (now a default binner) and `binning_wrapper.py`.
*   Added automatic gzipping of output files based on `.gz` extension in `edgelist_to_clusters.py` using `pyexeggutor.open_file_writer`.
*   Added `xxhash` dependency to `VEBA-binning-prokaryotic_env` for bin name reproducibility ([Issue #140](https://github.com/jolespin/veba/issues/140)).
*   Added `-e/--exclude` and `-d/--domain_predictions` options to `filter_binette_results.py` for removing eukaryotic genomes and setting up domain assignments ([Issue #153](https://github.com/jolespin/veba/issues/153)).
*   Added `semibin2-[biome]` option to `binning-prokaryotic.py` allowing specification of multiple biomes (e.g., `semibin2-global`, `semibin2-ocean`), replacing `--semibin2_biome` ([Issue #155](https://github.com/jolespin/veba/issues/155)).
*   Added `--semibin2_orf_finder` option to `binning_wrapper.py`.
*   Added `genome_statistics.tsv.gz`, `gene_statistics.cds.tsv.gz`, `gene_statistics.rRNA.tsv.gz`, and `gene_statistics.tRNA.tsv.gz` outputs to `essentials.py`.
*   Added `--identifiers`, `--index_name`, and `--no_header` options to `convert_metabat2_coverage.py` for broader applicability, including `VAMB`.
*   Added `-l eukaryota_odb12` as default but also allow `--auto-lineage-euk` for `BUSCO` in `binning-eukaryotic.py`

### Changed

*   Changed `binning-eukaryotic.py` behavior to provide a solution to [BUSCO Issue #447](https://gitlab.com/ezlab/busco/-/issues/447)
*   Changed `CHANGELOG.md` format to best practice [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
*   Changed `prodigal-gv` to `pyrodigal-gv` in multithreaded mode for `binning-viral.py` for performance.
*   Removed `metacoag` from the default set of binning algorithms in `binning-prokaryotic.py`.
*   Updated `geNomad` to `v1.11.0` and `geNomad database` to `v1.8` to resolve numpy import errors ([Issue #160](https://github.com/jolespin/veba/issues/160)).
*   Updated `Pyrodigal` usage in `binning-eukaryotic.py` for organelles to allow piping and threading.
*   Updated `BUSCO` to `v5.8.3` and associated databases.
*   Updated `Tiara` to `Tiara-NAL` in `VEBA-binning-prokaryotic_env` and `VEBA-binning-eukaryotic_env` to enable `stdin` usage.
*   Updated `biosynthetic.py` to use `antiSMASH v7` ([Issue #159](https://github.com/jolespin/veba/issues/159)).
*   Changed behavior when `--taxon fungi` is specified: precomputed genes are not used due to formatting issues.
*   Simplified the method for adding headers to `Diamond` outputs in `biosynthetic.py`.
*   Changed `Dockerfile` working directory from `/tmp/` to `/home/`.
*   Integrated `Tiara` and `consensus_domain_classification.py` into the `binette` step of `binning-prokaryotic.py`.
*   Renamed database identifier from `VDB` to `VEBA-DB`.
*   Updated `CheckM2` and `Binette` versions in `binning-prokaryotic.py`.
*   Updated `CheckM2 Diamond` database included in `VEBA-DB_v9` ([Issue #154](https://github.com/jolespin/veba/issues/154)).
*   Removed usage of precomputed genes in the `SemiBin2` wrapper due to [SemiBin2/issue-#185](https://github.com/BigDataBiology/SemiBin/issues/185).
*   Allowed faulty return codes in iterative mode for `binette` to permit convergence in genome recovery.

### Fixed

*   Fixed `CONDA_ENVS_PATH` detection in the `veba` controller executable to correctly handle environments outside the base Conda directory.
*   Fixed bug where `VFDB` hits were incorrectly counted as `MIBiG` in `biosynthetic.py` ([Issue #141](https://github.com/jolespin/veba/issues/141)).
*   Fixed `--tta_threshold` argument in `biosynthetic.py` which was previously defined but not connected to the command execution.
*   Removed capitalization from column headers in `filter_binette_results.py` output.
*   Fixed missing `--antismash_options` argument connection in `biosynthetic.py`.

### Removed

*   Removed `CONCOCT` support from `binning-eukaryotic.py`.

### Deprecated

*   Deprecated `amplicon.py` module in favor of external pipelines like `nf-core/ampliseq`.

## [2.4.2] - 2025.02.01

### Added

*   Added `--megahit_build_de_bruijn_graph` option to `assembly.py` to make de Bruijn graph construction optional for `MEGAHIT`.

### Fixed

*   Fixed an issue where the `MEGAHIT` de Bruijn graph was not included in the output directory even if constructed.

## [2.4.0] - 2025.01.29

### Added

*   Added `essentials.py` module.
*   Added `--serialized_annotations` option to `append_annotations_to_gff.py` to improve performance by avoiding re-parsing.
*   Added `--busco_options` and `--busco_offline` arguments to `binning-eukaryotic.py`.
*   Added `--semibin2_sequencing_type` argument to `binning_wrapper.py` and improved `--long_reads` handling.
*   Added support for precomputed coverage input for `metadecoder` in `binning_wrapper.py`.
*   Added support for `binette` and `tiara` to the `binning-prokaryotic.py` module.
*   Added `copy_attribute_in_gff.py` script to copy GFF attributes between fields.
*   Added `filter_binette_results.py` script.
*   Added support for `metacoag` binner and custom HMM support for `metadecoder` in `binning_wrapper.py`.
*   Added `prepend_de-bruijn_path.py` script and integrated into `assembly.py` and `assembly-long.py` to prefix SPAdes/Flye graph paths.
*   Added support for `SemiBin2` and `MetaDecoder` binners in `binning_wrapper.py`.
*   Added `update_genome_clusters.py` script to update genome clusters using `skani` against reference clusters (does not update protein clusters, graph, representatives, or proteins).
*   Added `Enzymes` annotation category support to `append_annotations_to_gff.py`.
*   **Database (`VDB_v8.1`)**: Added `kofam.enzymes.list` and `kofam.pathways.list` files for subsetting `pykofamsearch` queries. Added serialized KOfam database with enzyme support.
*   Added `Enzymes` annotation support to `annotate.py` and `merge_annotations.py`.
*   Added `--cluster_label_mode {"numeric", "random", "pseudo-random", "md5", "nodes"}` option to `edgelist_to_clusters.py` for flexible cluster labeling.
*   Added `--threshold2` option to `edgelist_to_clusters.py` for applying a secondary weight threshold.
*   Added `--wrap` option to `fasta_utility.py`.
*   Added `prepend_gff.py` script to prepend prefixes to GFF contig and attribute identifiers.
*   Added `append_annotations_to_gff.py` script.
*   Added `Initial_bins` information to `Binette` results in `filter_binette_results.py`.

### Changed

*   Changed default `--busco_completeness` threshold from `50` to `30` in `binning-eukaryotic.py`.
*   Moved `--long_reads` argument handling within `binning_wrapper.py`.
*   Changed default `--minimum_genome_size` from `150000` to `200000`.
*   Changed default `--cluster_label_mode` from `numeric` to `md5` in `cluster.py`, `global_clustering.py`, `local_clustering.py`, and `update_genome_clusters.py` for easier post hoc updates.
*   Changed `diamond` calls in `annotate.py` to use `--header simple` and updated `merge_annotations.py` accordingly.
*   **Database**: Updated VEBA database from `VDB_v8` to `VDB_v8.1`.
*   Updated `pyhmmsearch` and `pykofamsearch` versions in relevant environments (`VEBA-annotate_env`, `VEBA-classify-eukaryotic_env`, `VEBA-database_env`, `VEBA-phylogeny_env`) and updated executables/scripts (`annotate.py`, `classify-eukaryotic.py`, `phylogeny.py`, `download_databases-annotate.sh`).
*   Modified `fasta_utility.py` to split ID and description in FASTA headers, ensuring prefixes/suffixes apply only to the ID.
*   Changed default `--skani_minimum_af` from `15` to `50` in clustering scripts (`cluster.py`, `global_clustering.py`, `local_clustering.py`) to align with GTDB-Tk species-level clustering criteria.
*   Renamed `manual` mode in `binning-eukaryotic.py` to `metaeuk` mode for clarity when using preexisting `metaeuk` results.

### Fixed

*   Fixed bug in `binning_wrapper.py` where the script looked for `MetaCoAG` bins in the wrong directory.
*   Fixed bug in `merge_annotations.py` where `diamond` outputs were queried incorrectly.
*   Fixed issue in `consensus_domain_classification.py` where `softmax` returned a `numpy.array` instead of a `pandas.DataFrame`.
*   Fixed intermediate directory handling for `metacoag` in `binning_wrapper.py`.

## [2.3.0] - 2024.09.24

### Added

*   Added `${N_JOBS}` variable support to database download scripts, defaulting to maximum available threads.

### Changed

*   Replaced `MicrobeAnnotator-KEGG` with `KEGG Pathway Profiler` for KEGG module completion ratio calculations in `annotate.py`. Requires `VEBA-database_env` and `VEBA-annotate_env`.
*   Updated database structure: `${VEBA_DATABASE}/Annotate/MicrobeAnnotator-KEGG` replaced by `${VEBA_DATABASE}/Annotate/KEGG-Pathway-Profiler/`.
*   **Note**: The new module completion ratio output from `KEGG Pathway Profiler` does not include KEGG module class labels.

## [2.2.1] - 2024.09.24

### Added

*   Added `VERSION` file creation during `download_databases.sh` execution to track database version.
*   Added `--af_mode` option (`relaxed` or `strict`) to `edgelist_to_clusters.py`, `global_clustering.py`, `local_clustering.py`, and `cluster.py` to control how alignment fraction thresholds are applied to query and reference.
*   Added `-m/--include_mrna` option to `compile_metaeuk_identifiers.py` ([Issue #110](https://github.com/jolespin/veba/issues/110)).

### Changed

*   Renamed `choose_fastest_mirror.py` script to `determine_fastest_mirror.py`.

### Fixed

*   Fixed issue where genome clustering alignment fraction threshold (`--af_threshold`) was only applied to the reference genome; it now considers both reference and query based on the selected `--af_mode`.
*   Added missing `pigz` dependency to `VEBA-annotate_env`, required for Docker container execution.

## [2.2.0] - 2024.06.09

### Changed

*   Updated `GTDB-Tk` from version `2.3.0` to `2.4.0`.
*   Updated `GTDB` database requirement from `r214.1` to `r220`.
*   Split the `VEBA-classify_env` Conda environment into `VEBA-classify-eukaryotic_env`, `VEBA-classify-prokaryotic_env`, and `VEBA-prokaryotic_env`.
*   Rewrote `annotate.py`, `classify-eukaryotic.py`, `phylogeny.py`, and associated utility scripts to use `PyHMMER` (`pyhmmsearch`, `pykofamsearch`) for improved multithreaded performance over external `HMMER` calls.
*   Standardized database name `KOFAM` to `KOfam`.

## [2.1.0] - 2024.05.17

### Added

*   Added `concatenate_files.py` script for robust concatenation of many files (including mixed compressed/uncompressed), overcoming Unix argument limits.
*   Added `/volumes/workspace/` directory to Docker containers for use cases where input and output directories are the same.
*   Added `uniprot_to_enzymes.py` script to reformat data from UniProt enzyme queries.
*   Added `--metaeuk_split_memory_limit` parameter to `metaeuk_wrapper.py`.
*   Added `-d/--genome_identifier_directory_index` option to `scaffolds_to_bins.py` for handling specific directory structures (e.g., `path/to/genomes/[bin_id]/reference.fasta`).
*   Added `--minimum_af` option to `edgelist_to_clusters.py` supporting 4-column input (`id1`, `id2`, `weight`, `alignment_fraction`).
*   Added `--export_representatives` option to `edgelist_to_clusters.py` to output representative nodes per cluster along with connectivity information (also added to NetworkX graph object output).
*   Added official YouTube channel link: [VEBA Multiomics Channel](https://www.youtube.com/@VEBA-Multiomics).
*   (Internal Development) Developed faster CLI implementations `PyKofamSearch` and `PyHMMSearch` based on `PyHMMER` for future integration.

### Changed

*   Changed default behavior in `global_clustering.py`, `local_clustering.py`, and `cluster.py` to use an alignment fraction threshold (`--af_threshold 30.0`) by default via `edgelist_to_clusters.py`. Use `--af_threshold 0.0` to retain previous behavior.
*   Modified `edgelist_to_clusters.py` node filtering: if `--identifiers` are provided, only edges between those identifiers are considered; otherwise, all nodes are included.
*   Changed singleton weight representation in `edgelist_to_clusters.py` from `np.inf` to `np.nan` to facilitate representative calculations.

### Fixed

*   Limited `featureCounts` to use a maximum of 64 threads (`min(64, n_jobs)`) in relevant modules/scripts, as the tool cannot utilize more.

## [2.0.0] - 2024.03.10

### Added

*   Added `-A/--from_antismash` option in `biosynthetic.py` to use preexisting `antiSMASH` results.
*   Added `number_of_genomes`, `number_of_genome-clusters`, `number_of_proteins`, and `number_of_protein-clusters` columns to `feature_compression_ratios.tsv.gz` output from `cluster.py`.
*   Added support for custom paths for Conda environments via `CONDA_ENVS_PATH` variable in `update_environment_scripts.sh`, `update_environment_variables.sh`, `check_installation.sh`, `install.sh`, and core `veba` script.
*   Added `busco_version` parameter to `merge_busco_json.py` with default `5.4.x` and support for `5.6.x`.
*   Added `merge_counts_with_taxonomy.py` script.

### Changed

*   Changed `-i/--input` parameter to `-i/--from_genomes` in `biosynthetic.py`.
*   Renamed `antimash_genbanks_to_table.py` to `biosynthetic_genbanks_to_table.py` to prepare for future support of `DeepBGC` and `GECCO`.
*   Changed default assembly algorithm in `assembly-long.py` from `flye` to `metaflye`.
*   Modified `install.sh` to accept custom `CONDA_ENVS_PATH` argument: `bash install.sh path/to/log path/to/envs/`.

## [1.5.0] - 2024.01.30

### Added

*   Added `VeryFastTree` as an option to `phylogeny.py` (`FastTree` remains default).
*   Added `--blacklist` option to `compile_eukaryotic_classifications.py` (default: `species:uncultured eukaryote` in `classify_eukaryotic.py`).
*   Added `compile_phylogenomic_functional_categories.py` script automating methodology from [Espinoza et al. 2022](https://doi.org/10.1093/pnasnexus/pgac239).

### Changed

*   Updated `antismash_genbanks_to_table.py` for compatibility with `antiSMASH v6` and `v7` genbanks (fixed `.strip` method usage).
*   Changed default behavior in `assembly-long.py` to be non-deterministic for increased speed (changed `--no_deterministic` flag to `--deterministic`).
*   Changed default `--leniency` parameter on `classify_eukaryotic.py` and `consensus_genome_classification_ranked.py` to `1.0`. Added separate `--leniecy_genome_classification` option.
*   Refactored installation files and structure (`veba/src/` -> `veba/bin/`).
*   Updated environment variable checks from `SCRIPT_VERSIONS` to `VEBA_SCRIPT_VERSIONS` (now located in `bin/` of conda environment).

### Fixed

*   Fixed critical error where `classify_eukaryotic.py` accessed a deprecated database file from `MicroEuk_v2`.
*   Fixed incorrect header offset in `annotations.protein_clusters.tsv` preventing Pandas from reading it correctly.
*   Fixed issue in `binning-prokaryotic.py` where empty or non-existent `unbinned.fasta` files were not handled correctly, leading to bad symlinks for GFF, rRNA, and tRNA files when no valid bins were found.
*   Fixed creation of non-existing symlinks (`*.gff`, `*.rRNA`, `*.tRNA`) in `binning-prokaryotic.py`.
*   Fixed minor error in `coverage.py` where `samtools sort --reference` incorrectly received `reads_table.tsv` instead of `reference.fasta`.
*   Fixed minor error preventing `Tiara` from running in the backend of `eukaryotic_gene_modeling_wrapper.py`.

## [1.4.2] - 2023.12.21

### Added

*   Added `profile-taxonomic.py` module using `sylph` for sketch-based taxonomic abundance profiling.
*   Added support for long reads in `fastq_preprocessor`, `preprocess.py`, `assembly-long.py`, `coverage-long`, and all binning modules.
*   Added new usage syntax `veba --module <name> --params "${PARAMS}"` abstracting Conda environment activation. Updated walkthroughs accordingly.
*   Added `skani` as the new default algorithm for Average Nucleotide Identity (ANI) based genome clustering (`cluster.py`, `global_clustering.py`, `local_clustering.py`).
*   Added `Diamond DeepClust` as an alternative protein clustering method to `MMSEQS2` (`clustering_wrapper.py`, `global/local_clustering.py`, `cluster.py`).
*   Added `--reference_gzipped` option to `index.py` and `mapping.py` (default assumes uncompressed reference).
*   Added support for missing values in `compile_eukaryotic_classifications.py`.
*   Added `--metaeuk_split_memory_limit` argument (default `36G`, experimental) to `binning-eukaryotic.py` and `eukaryotic_gene_modeling.py`.
*   Added check in `check_fasta_duplicates.py` and `clean_fasta.py` for `>` characters within sequences.
*   Added `sylph` tool to `VEBA-profile_env`.

### Changed

*   **Database (`VDB_v6`)**: Rebuilt `VEBA Microeukaryotic Protein Database` into clustered versions (`MicroEuk100/90/50`), available at [doi:10.5281/zenodo.10139450](https://zenodo.org/records/10139451).
*   Redesigned `binning-eukaryotic` module to better handle custom `MetaEuk` databases.
*   Removed requirement for `--estimated_assembly_size` in `assembly-long.py` for Flye ([Flye Issue #652](https://github.com/fenderglass/Flye/issues/652)).
*   Enabled sequence compression (`--compressed 1`) during `mmseqs createdb` in `download_databases.sh`.
*   Renamed `mmseqs2_wrapper.py` to `clustering_wrapper.py`.
*   Renamed `MMSEQS2` commands `easy-cluster` and `easy-linclust` to `mmseqs-cluster` and `mmseqs-linclust` within wrappers.
*   Renamed `consensus_genome_classification.py` to `consensus_genome_classification_ranked.py` and changed default behavior to allow missing taxonomic levels.

### Fixed

*   Fixed critical error where `classify_eukaryotic.py` accessed a deprecated database file from `MicroEuk_v2`. (Note: Also fixed in 1.5.0, perhaps a regression or different aspect?)
*   Fixed dereplication step for duplicate contigs in `concatenate_fasta.py`.
*   Fixed formatting of `annotations.protein_clusters.tsv.gz` generated by `merge_annotations.py` (issue added in v1.3.1 patch).
*   Fixed calculation of viral quality scores in `merge_genome_quality_assessments.py`.
*   Fixed memory leak in `merge_annotations.py` when creating `annotations.protein_clusters.tsv.gz`. (Note: formatting for empty sets/lists still needed correction at time of release).

## [1.3.0] - 2023.10.27

### Added

*   **Modules & Scripts:**
    *   Added `profile-pathway.py` module for reads-based functional profiling using `HUMAnN` with custom databases built from VEBA outputs.
    *   Added `marker_gene_clustering.py` script to identify core and unique marker proteins within genome clusters (pangenomes).
    *   Added `module_completion_ratios.py` script to calculate KEGG module completion ratios (run automatically by `annotate.py`).
    *   Added `merge_genome_quality.py` to compile `CheckV`, `CheckM2`, and `BUSCO` results.
    *   Added `merge_taxonomy_classifications.py` to compile taxonomy classifications.
    *   Added Biosynthetic Gene Cluster (BGC) clustering (protein and nucleotide) to `biosynthetic.py`, including prevalence tables.
    *   Added `pangenome_core_sequences` output (protein and CDS) to `cluster.py`.
    *   Added PDF visualization rendering of Newick trees in `phylogeny.py` using `ete3`.
    *   Added `--minimum_core_prevalence` parameter to clustering scripts (`global_clustering.py`, `local_clustering.py`, `cluster.py`).
    *   Added checkpoint for `tRNAscan-SE` in `binning-prokaryotic.py` and `eukaryotic_gene_modeling_wrapper.py`.
    *   Added `compile_custom_humann_database_from_annotations.py` script.
    *   Added `--no_domain` and `--no_header` options to `merge_taxonomy_classifications.py`.
    *   Added check in `coverage.py` to skip processing if `mapped.sorted.bam` files already exist (GNU parallel option not yet implemented).
    *   Added `--nucleotide_fasta_output` to `biosynthetic_genbanks_to_table.py` (previously `antismash_genbanks_to_table.py`).
    *   Added support for `MEGAHIT` presets via `--megahit_preset` option.
    *   Added `--sort_by` and `--ascending` options to `concatenate_dataframes.py` with automatic duplicate index removal. Set `--sort_by bitscore` in `biosynthetic.py`.
    *   Added core pangenome and singleton hits information to clustering output tables.

*   **Database (`VDB_v5.2`)**:
    *   Added `CAZy` database support.
    *   Added `MicrobeAnnotator-KEGG` database support (including module definitions from [Zenodo:10020074](https://zenodo.org/records/10020074)).

### Changed

*   Updated `annotate.py` and `merge_annotations.py` to handle `CAZy` and improve clustered protein annotations.
*   Changed default representative sequence format from table to fasta in `clustering_wrapper.py` (previously `mmseqs2_wrapper.py`).
*   Renamed `--fasta_output` to `--protein_fasta_output` in `biosynthetic_genbanks_to_table.py`.
*   Updated BGC component identifiers in `biosynthetic_genbanks_to_table.py` to match `MetaEuk` format: `[bgc_id]_[position_in_bgc]|[start]:[end]([strand])`.
*   Renamed `bgc_type` to `protocluster_type` in `biosynthetic.py` output.
*   Enabled `biosynthetic.py` to process GFF files from `MetaEuk`.
*   Relabeled `--input` to `--genomes_table` in clustering scripts/module.
*   Changed default behavior in `cluster.py` to *not* perform local clustering (changed `--no_local_clustering` flag to `--local_clustering`).
*   Updated default `--megahit_memory` from `0.9` to `0.99`.
*   Updated `bowtie2`, `fastp`, and `fastq_preprocessor` versions in relevant module environments.

### Fixed

*   Fixed `pandas.errors.InvalidIndexError` in `biosynthetic.py` caused by multiple `Diamond` hits in one region by adding sorting and duplicate index removal via `concatenate_dataframes.py`.
*   Fixed `antismash_genbanks_to_table.py` failure when `antiSMASH` added CDS features (e.g., `allorf_[start]_[end]`) not present in the GFF.
*   Corrected cluster-level taxonomy weighting in `compile_prokaryotic_genome_cluster_classification_scores_table.py` when using `--mash_db` with `GTDB-Tk` by using `fastani_ani` or `msa_percent` appropriately.
*   Fixed creation of empty GFF files with asterisk names for samples without prokaryotic MAGs.
*   Fixed critical error in `binning-prokaryotic.py` where descriptions in `eukaryota.scaffolds.list` headers were not removed, causing `seqkit grep` to fail and `DAS_Tool` to incorrectly include eukaryotic MAGs.
*   Fixed incorrect `krona.html` generation in `biosynthetic.py` from `compile_krona.py`.
*   Fixed error in `genomad_taxonomy_wrapper.py` using incorrect filename (`viral_taxonomy.tsv` instead of `taxonomy.tsv`).
*   Fixed minor string formatting error in `assembly.py` preventing use of `SPAdes` programs other than `spades.py`, `metaspades.py`, or `rnaspades.py`.

### Removed

*   Removed `--no_singletons` option from `cluster.py`.

## [1.2.0] - 2023.07.11

### Added

*   Added rRNA (`BARRNAP`) and tRNA (`tRNAscan-SE`) identification for prokaryotic and eukaryotic genomes.
*   Added `compile_gff.py` script to merge CDS, rRNA, and tRNA GFF files, including GC content and organelle tags (used in binning modules).
*   Added `eukaryotic_gene_modeling_wrapper.py` for comprehensive eukaryotic gene prediction (split organelles, run `MetaEuk`/`Pyrodigal`, `BARRNAP`, `tRNAscan-SE`, merge GFFs, calculate stats).
*   Added automatic generation of pangenome protein prevalence tables and singularity ratio calculations in `cluster.py`.
*   Added [Virulence Factor Database (`VFDB`)](http://www.mgc.ac.cn/VFs/main.htm) annotation support in `annotate.py`.
*   Added [UniRef50/90](https://www.uniprot.org/help/uniref) annotation support in `annotate.py`.
*   Added `Krona` plot generation for taxonomy classifications and BGC detection.
*   Added `gene_biotype=protein_coding` attribute to `Pyrodigal(-GV)` GFF output.
*   Added `gtdb_r214.msh` (Mash sketch) for faster GTDB-Tk classification via ANI screening ([Zenodo:8048187](https://zenodo.org/record/8048187)).
*   Added support for missing classifications in `compile_krona.py` and `consensus_genome_classification_ranked.py`.

### Changed

*   **Database (`VDB_v5.1`)**:
    *   Added `VFDB`.
    *   Updated `GTDB r207_v2` to `r214.1`.
    *   Replaced `NR` annotations with `UniRef50/90`.
*   Updated `GTDB-Tk` classification in `classify-prokaryotic.py` to use Mash screening by default (`gtdb_r214.msh`).
*   Updated `GTDB-Tk` version `2.1.3` -> `2.3.0`.
*   Changed GTDB-Tk database path structure: `${VEBA_DATABASE}/Classify/GTDBTk` -> `${VEBA_DATABASE}/Classify/GTDB`.
*   Modified eukaryotic gene modeling to analyze nuclear, mitochondrial, and plastid sequences separately.
*   Updated genome GFF files to include contigs, CDS, rRNA, tRNA, with mitochondrial/plastid tags where applicable.
*   Updated `BUSCO` version `5.3.2` -> `5.4.3` and adapted `filter_busco_results.py` for new JSON output structure.
*   Replaced `prodigal` with `pyrodigal`.
*   Updated Docker images to define `/volumes/input`, `/volumes/output`, `/volumes/database` mount points.
*   Cleaned up intermediate files generated by global and local clustering scripts.

### Fixed

*   Fixed minor issue in `biosynthetic.py` where fasta and genbank files were not properly symlinked; added virulence factor results to synopsis output.
*   Fixed minor error in `binning-prokaryotic.py` where `--veba_database` argument was ignored (only `VEBA_DATABASE` environment variable worked).

### Deprecated

*   **Database (`VDB_v5.1`)**: Deprecated `RefSeq non-redundant protein` database in favor of `UniRef50/90`.

## [1.1.2] - 2023.05.16

### Added

*   Created Docker images for all VEBA modules.
*   Added `propagate_annotations_from_representatives.py` script.
*   Added `scaffolds_to_mags.tsv` file to clustering output.
*   Added `convert_counts_table.py` script to convert counts tables to Pandas pickle, Anndata h5ad, or Biom hdf5 formats.
*   Added support in `index.py` to accept individual `--references [file.fasta]` and `--gene_models [file.gff]` arguments.
*   Added `stdin` support and progress bars to `scaffolds_to_bins.py`, along with ability to process genome tables (`[id_genome]<tab>[filepath]`).
*   Added `check_fasta_duplicates.py` script (returns exit code 0/1).
*   Added `reformat_representative_sequences.py` script to format `MMSEQS2` representative sequences.
*   Added `mmseqs2_wrapper.py` and `hmmer_wrapper.py` helper scripts.
*   Added option to `merge_generalized_mapping.py` to include sample index in filepaths and remove empty features.
*   Added `genbanks/[id_genome]/` output directory in `biosynthetic.py` containing symlinks to BGC genbanks from `antiSMASH`.

### Changed

*   Replaced all absolute path symlinks with relative symlinks (also in `fastq_preprocessor`).
*   Standardized taxonomy output filenames: `prokaryotic_taxonomy.tsv` -> `taxonomy.tsv`, `prokaryotic_taxonomy.clusters.tsv` -> `taxonomy.clusters.tsv` (similarly for eukaryotic and viral).
*   Updated all Conda environments to use GenoPype `2023.4.13`.
*   Changed `nr` annotation source to `uniref` in `annotate.py`.
*   Simplified `merge_annotations_and_taxonomy.py` to `merge_annotations.py` (taxonomy operations removed).
*   Changed `orfs_to_orthogroups.tsv` to `proteins_to_orthogroups.tsv` for consistency.
*   Removed `python` prefix for script calls; now relies on shebang (`#!/usr/bin/env python`). Added single quotes around script filepaths in calls to handle spaces/special characters.
*   Updated `assembly.py`, `assembly-sequential.py`, `binning-*.py`, and `mapping.py` to use `featureCounts -p --countReadPairs` (addresses [issue #22](https://github.com/jolespin/veba/issues/22)). Updated `subread` version `2.0.1` -> `2.0.3`. Added `--long_reads` flag to `binning-*.py` for use with long read mapping.
*   Updated `cluster.py` (and associated `global/local_clustering.py`) to use `mmseqs2_wrapper.py` which automatically outputs representative sequences.
*   Added `executable='/bin/bash'` to `subprocess.Popen` calls in GenoPype (addresses [issue #23](https://github.com/jolespin/veba/issues/23)).

### Fixed

*   Fixed output directory structure for `mapping.py` to use `output_directory/${NAME}` like binning modules.

### Removed

*   Removed `--dbtype` argument from `global_clustering.py` and `local_clustering.py`.
*   Removed appended prefix for `.graph.pkl` and `dict.pkl` outputs in `edgelist_to_clusters.py`.

### Database (`VDB_v5`)

*   Changed `nr` database component to `UniRef90` and `UniRef50`.

## [1.1.1] - 2023.03.20

### Added

*   Added `MIBiG` database support to `annotate.py` and database downloads.
*   Added composite label generation for annotations in `annotate.py`.
*   Added `--dastool_minimum_score` parameter to `binning-prokaryotic.py`.
*   Added wrapper script for the `STAR` aligner.
*   Added option `--no_header` to `subset_table.py`.

### Changed

*   Updated `merge_generalized_mapping.py` script to accept BAM files directly instead of relying on a specific directory structure.

### Fixed

*   Fixed broken `VEBA-binning-viral.yml` Conda environment recipe due to `aria2` package conflicts ([commit 30e8b0a](https://github.com/jolespin/veba/commit/30e8b0a)).
*   Fixed issues with Conda environment variables in installation scripts.

## [1.1.0] - 2023.03.02

### Added

*   **Annotations (`annotate.py`):**
    *   Added `NCBIfam-AMRFinder` AMR domain annotations.
    *   Added `AntiFam` contamination annotations.
*   **Assembly (`assembly.py`):**
    *   Added `transcripts_to_genes.py` script to create `genes_to_transcripts.tsv` for `TransDecoder`.
*   **Binning (`binning-prokaryotic.py`, `binning-viral.py`):**
    *   Added `--skip_concoct` option to `binning-prokaryotic.py`.
    *   Added multi-split binning capability to `binning_wrapper.py` (concatenate contigs, map all reads, bin together, partition by sample).
*   **Classification (`classify-prokaryotic.py`, `classify-eukaryotic.py`, `classify-viral.py`):**
    *   Added `--genomes` option to classification modules (`classify-prokaryotic.py`, `classify-eukaryotic.py`, `classify-viral.py`) to classify external genomes not binned by VEBA.
    *   Implemented use of `eukaryota_odb10` markers for classifying external eukaryotic genomes, improving performance.
*   **Clustering (`cluster.py`):**
    *   Added optional (on by default) local-level clustering (within each sample) alongside global clustering.
    *   Added automatic calculation of genomic and functional Feature Compression Ratios (FCR) at global and local levels.
    *   Added support for zfilling cluster labels (e.g., `SLC7` -> `SLC007`).
*   **Scripts:**
    *   Added `append_geneid_to_transdecoder_gff.py`.
    *   Added `bowtie2_wrapper.py`.
    *   Added `compile_genomes_table.py`.
    *   Added `consensus_genome_classification_unranked.py`.
    *   Added `cut_table.py`.
    *   Added `cut_table_by_column_labels.py`.
    *   Added `drop_missing_values.py`.
    *   Added `edgelist_to_clusters.py` (generalizing `fastani_to_clusters.py`).
    *   Added `filter_checkm2_results.py`.
    *   Added `genomad_taxonomy_wrapper.py`.
    *   Added `global_clustering.py`.
    *   Added `local_clustering.py`.
    *   Added `partition_multisplit_bins.py`.
    *   Added `scaffolds_to_clusters.py`.
    *   Added `scaffolds_to_samples.py`.
    *   Added `transdecoder_wrapper.py` (requires separate environment).
*   **Miscellaneous:**
    *   Added versioning to Conda environment files.
    *   Added `mamba` to installation process for speed.
    *   Added support for `n_jobs = -1` to use all available threads.

### Changed

*   **Annotations (`annotate.py`):**
    *   Replaced `ete3` with `taxopy` in backend (`merge_annotations_and_score_taxonomy.py`).
*   **Binning (`binning-prokaryotic.py`, `binning-viral.py`):**
    *   Updated `CheckM` to `CheckM2` in `binning-prokaryotic.py`, significantly reducing resource requirements and simplifying the pipeline.
    *   Rewrote `binning-viral.py` module: uses `geNomad` as default (supports `VirFinder`), uses `Prodigal-GV` for viral genetic codes.
*   **Biosynthesis (`biosynthetic.py`):**
    *   Introduced unique and parseable `component_id` and `bgc_id` formats.
*   **Classification (`classify-prokaryotic.py`, `classify-viral.py`):**
    *   Updated `GTDB-Tk` from `v2.1.1` to `v2.2.3`. `--skip_ani_screen` is default for now.
    *   Rewrote `classify-viral.py` module: uses `geNomad` taxonomy and `consensus_genome_classification_unranked.py`.
*   **Clustering (`cluster.py`):**
    *   Rewrote module: replaced `OrthoFinder` with `MMSEQS2` for protein clustering, drastically reducing runtime and resource usage.
    *   Changed input format to a table: `[organism_type]<tab>[id_sample]<tab>[id_mag]<tab>[genome]<tab>[proteins]`.
    *   Renamed SLC-specific orthogroups (SSO) to SLC-specific protein clusters (SSPC).
*   **Phylogeny (`phylogeny.py`):**
    *   Updated `MUSCLE` to `v5` with new algorithms (`-align`, `-super5`) accessible via `--alignment_algorithm`. Input fasta files are no longer gzipped.
*   **Database (`VDB_v3.1` -> `VDB_v4`)**:
    *   Updated `CheckV DB v1.0` -> `v1.5`.
    *   Added `geNomad DB v1.2`.
    *   Added `CheckM2 DB`.
    *   Added `reference.eukaryota_odb10.list` and corresponding `MMSEQS2` database (`microeukaryotic.eukaryota_odb10`).
    *   Added `NCBIfam-AMRFinder` marker set.
    *   Added `AntiFam` marker set.
    *   Gzipped all HMM marker sets.
*   **Scripts:**
    *   Updated `antismash_genbanks_to_table.py`: added option for BGC fasta output, improved BGC identifiers.
    *   Updated `binning_wrapper.py` to handle multi-split binning.
    *   Updated `compile_reads_table.py`: file extension matching excludes the `.` for consistency.
    *   Updated `consensus_genome_classification.py`: output format matched to `consensus_genome_classification_unranked.py`.
    *   Updated `filter_checkv_results.py`: added option to use `geNomad` taxonomy/summaries.
    *   Updated `scaffolds_to_bins.py`: added `--genomes` argument support.
    *   Updated `subset_table.py`: added options to set index column and drop duplicates.
    *   Renamed `VirFinder_wrapper.R` to `virfinder_wrapper.r`, added option to use FDR instead of P-values.
    *   Rewrote `merge_annotations_and_score_taxonomy.py` using `taxopy`.
    *   Updated `merge_msa.py`: outputs uncompressed fasta by default (`--gzip` flag available).

### Removed

*   **Database (`VDB_v4`)**:
    *   Removed `CheckM DB`.
    *   Removed `taxa.sqlite` and `taxa.sqlite.traverse.pkl`.

### Deprecated

*   Deprecated backend scripts related to `CheckM` CPR handling in `binning-prokaryotic.py`.
*   Deprecated `adjust_genomes_for_cpr.py`.
*   Deprecated `filter_checkm_results.py`.
*   Deprecated `fastani_to_clusters.py` (replaced by `edgelist_to_clusters.py`).
*   Deprecated `partition_orthogroups.py`.
*   Deprecated `partition_clusters.py`.
*   Deprecated `compile_viral_classifications.py`.
*   Deprecated `build_taxa_sqlite.py`.

## [1.0.4] - 2022.12.27

### Changed

*   **Database (`VDB_v3`)**: Updated `Microeukaryotic Protein Database` (`VDB-Microeukaryotic_v2`), excluding some higher eukaryotes and changing identifiers to hashes. Moved hosting from FigShare to [Zenodo](https://zenodo.org/record/7485114) ([commit 0845ba6](https://github.com/jolespin/veba/commit/0845ba6be65f3486d61fe7ae21a2937efeb42ee9)).

### Fixed

*   Added missing `biopython` dependency to `VEBA-assembly_env` required when running `MEGAHIT` ([issue #17](https://github.com/jolespin/veba/issues/17), [commit aea51c3](https://github.com/jolespin/veba/commit/aea51c3e0b775aec90f7343f01cad6911f526f0a)).

## [1.0.3e] - 2022.12.14

### Added

*   Added `biosynthetic.py` module for running `antiSMASH` and converting genbanks to tables.
*   Added `megahit` assembler support to `assembly.py` module.
*   Added `--use_mag_as_description` option to `partition_gene_models.py` (default in `binning-prokaryotic.py`).
*   Added `adjust_genomes_for_cpr.py` script to simplify CPR adjustment step in `binning-prokaryotic.py`.
*   Added support for including MAG identifier in fasta header descriptions in `binning-prokaryotic.py` and `binning-eukaryotic.py`.
*   Added functionality to `replace_fasta_descriptions.py` to use a fixed string for replacement.

### Changed

*   Renamed `assembly.py` parameter `-P/--spades_program` to `-P/--program`.
*   Replaced complex bash commands in `binning-prokaryotic.py`'s penultimate step with `adjust_genomes_for_cpr.py`.
*   Updated `binning-eukaryotic.py` to use `metaeuk_wrapper.py` script for the `MetaEuk` step.
*   Changed `assembly-sequential.py` to run only `biosyntheticSPAdes` and `metaSPAdes` by default; `--run_metaplasmidspades` is now optional.

### Fixed

*   Fixed Conda environment creation error for `VEBA-assembly_env` in `install_veba.sh` ([issue #15](https://github.com/jolespin/veba/issues/15), [commit c2ab957](https://github.com/jolespin/veba/commit/c2ab957be132d34e6b99d6dea394be4572b83066)).
*   Fixed R error in `VirFinder_wrapper.R` caused by `__version__ = ` variable ([issue #13](https://github.com/jolespin/veba/issues/13), [commit 19e8f38](https://github.com/jolespin/veba/commit/19e8f38a5050328b7ba88b2271f0221073748cbb)).
*   Fixed error in `filter_busco_results.py` producing empty `identifier_mapping.metaeuk.tsv` tables ([issue #12](https://github.com/jolespin/veba/issues/12), [commit 359e456](https://github.com/jolespin/veba/commit/359e45699fc6d6fdf739350263fd34c6e4a62f94)).
*   Fixed Python error in `compile_metaeuk_identifiers.py` when duplicate gene identifiers were present ([issue #11](https://github.com/jolespin/veba/issues/11), [commit c248527](https://github.com/jolespin/veba/commit/c248527da9edef5ba2ebee348d707d8ece29fbee)).
*   Fixed Conda environment creation error for `VEBA-preprocess_env` in `install_veba.sh` ([issue #10](https://github.com/jolespin/veba/issues/10), [commit 8ed6eea](https://github.com/jolespin/veba/commit/8ed6eeaee1037694cf324d8fa4da6190578b9688)).

## [1.0.2a] - 2022.10.26

### Added

*   Added experimental `amplicon.py` module for short-read ASV detection using QIIME2/DADA2 workflow.
*   Added advanced sample parsing functionality to `compile_reads_table.py`.
*   Added `sra-tools` to `VEBA-preprocess_env`.
*   Added missing `CHECKM_DATA_PATH` environment variable definition to `VEBA-binning-prokaryotic_env` and `VEBA-classify_env`.

### Changed

*   Updated `GTDB-Tk` version from `1.x` to `2.x` in `VEBA-binning-prokaryotic_env`.
*   Updated `GTDB-Tk` database from `R202` to `R207_v2`.
*   Changed default human reference genome from GRCh38 no-alt to T2T CHM13v2.0.

### Fixed

*   Fixed incorrect script symlinks created by `install_veba.sh`.

## [1.0.1] - 2022.10.20

### Fixed

*   Fixed fatal error in `binning-eukaryotic.py` ([commit 7c5addf](https://github.com/jolespin/veba/commit/7c5addf9ed6e8e45502274dd353f20b211838a41)).
*   Fixed minor file naming issue in `cluster.py` ([commit 5803845](https://github.com/jolespin/veba/commit/58038451dac0791899aa7fca3f9d79454cb9ed46)).
*   Fixed issue where leftover human genome `tar.gz` file was not removed during database download/configuration ([commit 5803845](https://github.com/jolespin/veba/commit/58038451dac0791899aa7fca3f9d79454cb9ed46)).

## [1.0.0] - 2022.9.07

### Added

*   Initial public release accompanying the *BMC Bioinformatics* publication ([doi:10.1186/s12859-022-04973-8](https://doi.org/10.1186/s12859-022-04973-8)).

________________________________________________________________

#### Path to `v3.0.0`:

**Check:**

* Start/end positions on `MetaEuk` gene ID might be off.
* Why are some HMM annotations [number_of_hits, ids] while Enzyme annotations show the reverse order?

**Critical:**

* Symlink genomes in `classify-prokaryotic.py` instead of copying genomes [v2.6.0]
* Cluster module doesn't symlink global directory correctly on MacOS [v2.6.0]
* Return code for `cluster.py` when it fails during global and local clustering is 0 but should be 1. [v2.6.0]
* Don't load all genomes, proteins, and cds into memory for clustering.  Also, too many files are opened with large genomics datasets. [v2.6.0]
* Genome checkpoints in `tRNAscan-SE` aren't working properly.

**Definitely:**

* Update `eukaryota_odb10` -> `eukaryota_odb12` in `Markers` database
* Modify `COMEbin` so it can be used in the same environment as `Binette` and `CheckM2`
* When annotating proteins, create a hash representation and a dictionary of redundant sequences to decrease search space [v2.6.0]
* Replace `Bowtie2` with `strobealign` and `Fairy` when applicable (i.e., `coverage`, `assembly`, `binning`, `index`, and `mapping` modules) [v2.6.0] 
* Number of plasmids (via `geNomad`) for each MAG. [v2.6.0]
* Add a `isolate-prokaryotic.py` and `isolate-eukaryotic.py` modules which assembles, calls genes, checks for contamination and if there is, then bins genomes, and quality assesses.
* Add `metadecoder` with `eukaryota_odb12` marker set for `binning-eukaryotic.py`
* Script to add `PyHMMSearch` annotations to `annotations.proteins.tsv.gz`
* Develop `add_taxonomy_to_annotations.py` script
* Develop method for building and curating HMM cutoffs (e.g., comparing against false positives)
* Add `frozenset` for proteins that has all of the database identifiers
* Add `--proteins` option to `classify-eukaryotic.py` which aligns proteins to `MicroEuk100.eukaryota_odb10` via `MMseqs2` and then proceeds with the pipeline.
* Build SQL databases from all results
* Remove `p__Arthropoda` from `MicroEuk` database
* Add number of unique protein clusters to `identifier_mapping.genomes.tsv.gz` in `cluster.py` to assess most metabolicly diverse representative.
* Add `BiNI` biosynthetic novelty index to `biosynthetic.py`
* `busco_wrapper.py` that relabels all the genes, runs analysis, then converts output to tsv.
* Use `pigz` instead of `gzip`
* Add representative to `identifier_mapping.proteins.tsv.gz`
* Use `aria2` in parallel instead of `wget`.
* Add support for `Salmon` in `mapping.py` and `index.py`.  This can be used instead of `STAR` which will require adding the `exon` field to `Pyrodigal` GFF file (`MetaEuk` modified GFF files already have exon ids). 


**Eventually (Yes)?:**

* `NextFlow` support
* Install each module via `bioconda`
* Consistent usage of the following terms: 1) dataframe vs. table; 2) protein-cluster vs. orthogroup.  Dataframes should refer to generic tables while tables refer to specifics like "genomes table".
* Add coding density to GFF files
* Run `cmsearch` before `tRNAscan-SE`
* DN/DS from pangeome analysis
* For viral binning, contigs that are not identified as viral via `geNomad -> CheckV` use with `vRhyme`.
* Add `vRhyme` to `binning_wrapper.py` and support `vRhyme` in `binning-viral.py`.
* Different binning modes such as `ensemble-algorithm` (default) but also allow for `ensemble-seed` which allows for the same algorithm to be run multiple times different random states

**...Maybe (Not)?**
* Swap [`TransDecoder`](https://github.com/TransDecoder/TransDecoder) for [`TransSuite`](https://github.com/anonconda/TranSuite)


________________________________________________________________


<details>
	<summary> <b>Daily Change Log:</b> </summary>

* [2025.4.11] - Added `install-gpu.sh` which installs GPU accelerated environments when applicable (i.e., `VEBA-binning-prokaryotic_env` and `VEBA-binning-viral_env`)
* [2025.4.11] - Added `Dockerfile-GPU` which is experimental
* [2025.4.11] - Changed `install.sh` so it only installs CPU-based environments [Issue #167](https://github.com/jolespin/veba/issues/167)
* [2025.4.11] - Changed `containerize_environments.sh` so it only installs CPU-based environments [Issue #167](https://github.com/jolespin/veba/issues/167)
* [2025.4.11] - Deprecated `VirFinder` algorithm in `binning-viral.py` so now only `geNomad` is supported
* [2025.4.10] - Use `-l eukaryota_odb12` as default but also allow `--auto-lineage-euk` for `BUSCO` in `binning-eukaryotic.py`
* [2025.4.9] - Subset the coverage table for `VAMB` which is needed for iterations > 1 since the coverage index contigs must match the input contigs exactly
* [2025.4.8] - Added `VAMB` to `binning-prokaryotic.py` which is included as one of the default binners
* [2025.4.8] - Addd `VAMB` to `binning_wrapper.py`
* [2025.4.8] - Added `convert_metabat2_coverage.py` with `--identifiers`, `--index_name`, and `--no_header` so it can work in multiple scenarios including `VAMB`
* [2025.4.8] - Added `open_file_writer` from `pyexeggutor` to `edgelist_to_clusters.py` to handle automatic gzipping if `.gz` extension is detected.
* [2025.4.8] - Allow faulty return codes in iterative mode for `binette` since this allows genome recovery to converge
* [2025.4.4] - Deprecated `amplicon.py` because `nf-core/ampliseq` works great.
* [2025.4.4] - Changed `prodigal-gv` to `pyrodigal-gv` in multithreaded mode for `binning-viral.py`
* [2025.4.4] - Removed `metacoag`from default binning algorithms now that `SemiBin2` can be used with different biomes. 
* [2025.4.4] - Changed `CONDA_ENVS_PATH=${CONDA_ENVS_PATH:-"$(conda info --base)/envs/"}` to `CONDA_ENVS_PATH=$(dirname $CONDA_PREFIX)` in `veba` controller executable to handle new syntax when `VEBA` environments are in different location than the base conda.
* [2025.4.4] - Updated `geNomad` to `v1.11.0` and `geNomad database` to `v1.8` because of `ImportError: numpy.core.multiarray failed to import` error [Issue #160](https://github.com/jolespin/veba/issues/160)
* [2025.4.3] - Updated `Pyrodigal` in `binning-eukaryotic.py` for organelles to allow piping and threading
* [2025.4.3] - Remove `CONCOCT` support in `binning-eukaryotic.py`
* [2025.4.3] - Updated `BUSCO v5.8.3` and associated database
* [2025.4.3] - Updated `Tiara` to `Tiara-NAL` in `VEBA-binning-prokaryotic_env` and `VEBA-binning-eukaryotic_env` to handle `stdin`
* [2025.4.3] - Updated `biosynthetic.py` to use `antiSMASH v7` [Issue #159](https://github.com/jolespin/veba/issues/159)
* [2025.4.3] - If `--taxon fungi` then precomputed genes are not used because of formatting issues
* [2025.4.3] - Added `genome_statistics.tsv.gz`, `gene_statistics.cds.tsv.gz`, `gene_statistics.rRNA.tsv.gz`, and `gene_statistics.tRNA.tsv.gz` to `essentials.py`
* [2025.4.3] - Fixed bug where `VFDB` hits were counted as `MIBiG` in `biosynthetic.py` [Issue #141](https://github.com/jolespin/veba/issues/141)
* [2025.4.3] - Added `--tta_threshold` to `biosynthetic.py`.  This argument was already available but not actually connected to the command.
* ~~[2025.4.3] - Added `--clusterblast_database` to `biosynthetic.py` [Issue #143](https://github.com/jolespin/veba/issues/143)~~
* [2025.4.3] - Simplified the method used for adding the header to `Diamond` outputs in `biosynthetic.py`
* [2025.4.2] - Changed `Dockerfile` working directory from `/tmp/` to `/home/`
* [2025.4.2] - Added `xxhash` dependency to `VEBA-binning-prokaryotic_env` which is used for bin name reproducibility [Issue #140](https://github.com/jolespin/veba/issues/140)
* [2025.4.2] - Removed capitalization of column headers in `filter_binette_results.py`
* [2025.3.31] - Added `-e/--exclude` and `-d/--domain_predictions` to `filter_binette_results.py` to allow the removal of eukaryotic genomes and setup the domain assignments for `barrnap` and `tRNAscan-SE`.  [Issue #153](https://github.com/jolespin/veba/issues/153)
* [2025.3.31] - `Tiara` is now added into the `binette` step of `binning-prokaryotic.py` along with `consensus_domain_classification.py`.
* [2025.3.30] - Changing database name from `VDB` to `VEBA-DB`
* [2025.3.30] - Updated `CheckM2` and `Binette` versions in `binning-prokaryotic.py` which also includes a new `CheckM2` `Diamond` database in `VEBA-DB_v9` [Issue #154](https://github.com/jolespin/veba/issues/154)
* [2025.3.30] - Added `semibin2-[biome]` to `binning-prokaryotic.py` and remove `--semibin2_biome` which allows a user to specify multiple biomes (e.g., `semibin2-global` and `semibin2-ocean`) [Issue #155](https://github.com/jolespin/veba/issues/155)
* [2025.3.30] - Added `--semibin2_orf_finder` to `binning_wrapper.py` and remove precomputed genes from `SemiBin2` wrapper because [SemiBin2/issue-#185](https://github.com/BigDataBiology/SemiBin/issues/185).
* [2025.3.11] - Added `--antismash_options` missing bug in `biosynthetic.py`
* [2025.2.1] - Added `--megahit_build_de_bruijn_graph` to make de-Bruijn graph construction for `MEGAHIT` optional in `assembly.py`
* [2025.1.24] - Added `Initial_bins` to `Binette` results in `filter_binette_results.py`
* [2025.1.23] - Added `essentials.py` module
* [2025.1.16] - Added `--serialized_annotations` to `append_annotations_to_gff.py` to avoid overhead from reparsing the annotations
* [2025.1.15] - Fixed bug in `binning_wrapper.py` where script was looking for bins in the wrong directory for `MetaCoAG`
* [2025.1.14] - Fixed bug in `merge_annotations.py` where `diamond` outputs were queried incorrectly
* [2025.1.5] - Change default `--busco_completeness` from `50` to `30` in `binning-eukaryotic.py`
* [2025.1.5] - Added `--busco_options` and `--busco_offline` arguments for `binning-eukaryotic.py`
* [2024.12.28] - Added `--semibin2_sequencing_type` to `binning_wrapper.py` and added functionality for `--long_reads`.  Moved `--long_reads` argument to `parser_io` instead of `parser_featurecounts`
* [2024.12.27] - Fixed issue in `consensus_domain_classification.py` where `softmax` returns a `np.array` instead of a `pd.DataFrame`
* [2024.12.26] - Added support for precomputed coverage for `metadecoder` in `binning_wrapper.py`
* [2024.12.26] - Added support for `binette` and `tiara` in updated `binning_prokaryotic.py` module
* [2024.12.23] - Added `copy_attribute_in_gff.py` script which copies attributes to a source and destination attribute
* [2024.12.17] - Added `filter_binette_results.py` script
* [2024.12.16] - Added intermediate directory to `metacoag` in `binning_wrapper.py`
* [2024.12.12] - Added `metacoag` support and custom HMM support to `metadecoder` in `binning_wrapper.py` 
* [2024.12.11] - Added `prepend_de-bruijn_path.py` script and use this in `assembly.py` and `assembly-long.py` to prepend prefix to SPAdes/Flye de Bruijn graph paths.
* [2024.12.10] - Changed default `--minimum_genome_size` to `200000` from `150000`
* [2024.12.9] - Added support for `SemiBin2` and `MetaDecoder` in `binning_wrapper.py`
* [2024.11.21] - Updated `--cluster_label_mode` default to `md5` instead of `numeric` to allow for easier cluster updates post hoc. Change reflected in `cluster.py`, `global_clustering.py`, `local_clustering.py`, and `update_genome_clusters.py`
* [2024.11.18] - Added `update_genome_clusters.py` which runs `skani` against all reference genome clusters.  Does not do protein clustering nor does it update the graph, representatives, or proteins.
* [2024.11.15] - Added `--header simple` to `diamond` output in `annotate.py` and accounted for change in `merge_annotations.py`
* [2024.11.11] - Added `Enzymes` to `append_annotations_to_gff.py` script
* [2024.11.9] - Added `kofam.enzymes.list` and `kofam.pathways.list` in `VDB_v8.1` to provide subsets for `pykofamsearch`
* [2024.11.8] - Updating VEBA database `VDB_v8` to `VDB_v8.1` which adds serialized KOfam with enzyme support
* [2024.11.8] - Added `Enzymes` to `annotate.py` and `merge_annotations.py` [!untested]
* [2024.11.7] - Updated `pyhmmsearch` and `pykofamsearch` version in `VEBA-annotate_env.yml`, `VEBA-classify-eukaryotic_env.yml`,`VEBA-database_env`, and `VEBA-phylogeny_env`.  Also updated executables in `annotate.py`, `classify-eukaryotic.py`,  `phylogeny.py`, and `download_databases-annotate.sh`.
* [2024.11.7] - In `edgelist_to_clusters.py`, added `--cluster_label_mode {"numeric", "random", "pseudo-random", "md5", "nodes"}` to allow for different types of labels.  Added `--threshold2` option for a second weight.
* [2024.11.7] - Added `--wrap` to `fasta_utility.py` and split id and descriptions in header so prefix/suffix is only added to id.
* [2024.11.7] - Added `prepend_gff.py` to prepend a prefix to contig and attribute identifiers
* [2024.11.7] - Changed default `--skani_minimum_af` to `50` from `15` as this is used in GTDB-Tk for determining species-level clusters in `cluster.py`, `global_clustering.py`, and `local_clustering.py`
* [2024.11.6] - Added `append_annotations_to_gff.py` script
* [2024.10.29] - Changed `manual` mode to `metaeuk` mode for preexisting `metaeuk` results
* [2024.9.21] - Added `KEGG Pathway Profiler` to `VEBA-database_env` and `VEBA-annotate_env` which replaces `MicrobeAnnotator-KEGG` for module completion ratios.  Replacing `${VEBA_DATABASE}/Annotate/MicrobeAnnotator-KEGG` with `${VEBA_DATABASE}/Annotate/KEGG-Pathway-Profiler/` database files.  **Note: New module completion ratio output does not have classes labels for KEGG modules.**
* [2024.8.30] - Added ${N_JOBS} to download scripts with default set to maximum threads available
* [2024.8.29] - Added `VERSION` file created in `download_databases.sh`
* [2024.7.11] - Alignment fraction threshold for genome clustering only applied to reference but should also apply to query.  Added `--af_mode` with either `relaxed = max([Alignment_fraction_ref, Alignment_fraction_query]) > minimum_af` or `strict = (Alignment_fraction_ref > minimum_af) & (Alignment_fraction_query > minimum_af)` to `edgelist_to_clusters.py`, `global_clustering.py`, `local_clustering.py`, and `cluster.py`.
* [2024.7.3] - Added `pigz` to `VEBA-annotate_env` which isn't a problem with most `conda` installations but needed for `docker` containers.
* [2024.6.21] - Changed `choose_fastest_mirror.py` to `determine_fastest_mirror.py`
* [2024.6.20] - Added `-m/--include_mrna` to `compile_metaeuk_identifiers.py` for [Issue #110](https://github.com/jolespin/veba/issues/110)
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









