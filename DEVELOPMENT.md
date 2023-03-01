### What's next for *VEBA*?

*VEBA* is currently under active development. If you are interested in requesting features or wish to report a bug, please post a GitHub issue prefixed with the tag `[Feature Request]` and `[Bug]`, respectively.  If you want to contribute or have any other inquiries, contact me at `jespinoz[A|T]jcvi[DOT]org`.


________________________________________________________________

#### Current Releases:
##### Release v1.1 (Currently testing before official release)

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
		* Complete rewrite of this module which now uses `MMSEQS2` as the orthogroup detection algorithm instead of `OrthoFinder`.  `OrthoFinder` is overkill for creating protein clusters and it generates thousands of intermediate files (e.g., fasta, alignments, trees, etc.) which substantially increases the compute time.  `MMSEQS2` has very similar performance with a fraction of the resources and compute time.  Clustered the entire [Plastisphere dataset](https://figshare.com/articles/dataset/Genome_assemblies_gene_models_gene_annotations_taxonomy_classifications_clusters_and_counts_tables/20263974) on a local machine in ~30 minutes compared to several days on a HPC.
		* Now that the resources are minimal, clustering is performed at global level as before (i.e., all samples in the dataset) and now at the local level, optionally but ON by default, which clusters all genomes within a sample.  Accompanying wrapper scripts are `global_clustering.py` and `local_clustering.py`.
		* The genomic and functional feature compression ratios (FCR) (described [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04973-8)]) are now calculated automatically.  The calculation is `1 - number_of_clusters/number_of_features` which can easily be converted into an unsupervised biodiversity metric.  This is calculated at the global (original implementation) and local levels.
		* Input is now a table with the following columns: `[organism_type]<tab>[id_sample]<tab>[id_mag]<tab>[genome]<tab>[proteins]` and is generated easily with the `compile_genomes_table.py` script.  This allows clustering to be performed for prokaryotes, eukaryotes, and viruses all at the same time.
		* SLC-specific orthogroups (SSO) are now refered to as SLC-specific protein clusters (SSPC).
		* Support zfilling (e.g., `zfill=3, SLC7 → SLC007`) for genomic and protein clusters.
		* Deprecated `fastani_to_clusters.py` to now use the more generalizable `edgelist_to_clusters.py` which is used for both genomic and protein clusters.  This also outputs a `NetworkX` graph and a pickled dictionary `{"cluster_a":{"component_1", "component_2", ..., "component_n"}}`
		
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
		* Marker sets HMMs are now all gzipped (previously could not gzip because CheckM CPR workflow)

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


##### Release v1.0.4
* Added `biopython` to `VEBA-assembly_env` which is needed when running `MEGAHIT` as the scaffolds are rewritten and [an error](https://github.com/jolespin/veba/issues/17) was raised. [aea51c3](https://github.com/jolespin/veba/commit/aea51c3e0b775aec90f7343f01cad6911f526f0a)
* Updated Microeukaryotic protein database to exclude a few higher eukaryotes that were present in database, changed naming scheme to hash identifiers (from `cat reference.faa | seqkit fx2tab -s -n > id_to_hash.tsv`).  Switching database from [FigShare](https://figshare.com/articles/dataset/Microeukaryotic_Protein_Database/19668855) to [Zenodo](https://zenodo.org/record/7485114#.Y6vZO-zMKDU).  Uses database version `VDB_v3` which has the updated microeukaryotic protein database (`VDB-Microeukaryotic_v2`) [0845ba6](https://github.com/jolespin/veba/commit/0845ba6be65f3486d61fe7ae21a2937efeb42ee9)

___

##### Release v1.0.3e
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

___


##### Release v1.0.2a
* Updated *GTDB-Tk* in `VEBA-binning-prokaryotic_env` from `1.x` to `2.x` (this version uses much less memory): [f3507dd](https://github.com/jolespin/veba/commit/f3507dd13a42960e3671c9f8a106c9974fbfce21)
* Updated the *GTDB-Tk* database from `R202` to `R207_v2` to be compatible with *GTDB-Tk v2.x*: [f3507dd](https://github.com/jolespin/veba/commit/f3507dd13a42960e3671c9f8a106c9974fbfce21)
* Updated the [GRCh38 no-alt analysis set](https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip) to [T2T CHM13v2.0](https://genome-idx.s3.amazonaws.com/bt/chm13v2.0.zip) for the default human reference: [5ccb4e2](https://github.com/jolespin/veba/commit/5ccb4e20564513707fcc3420b18237974455e196)
* Added an experimental `amplicon.py` module for short-read ASV detection via the *DADA2* workflow of *QIIME2*: [cd4ed2b](https://github.com/jolespin/veba/commit/cd4ed2bfe35d5379a63dd3294c229f2c861f6f77)
* Added additional functionality to `compile_reads_table.py` to handle advanced parsing of samples from fastq directories while also maintaining support for parsing filenames from `veba_output/preprocess`: [cd4ed2b](https://github.com/jolespin/veba/commit/cd4ed2bfe35d5379a63dd3294c229f2c861f6f77)
* Added `sra-tools` to `VEBA-preprocess_env`: [f3507dd](https://github.com/jolespin/veba/commit/f3507dd13a42960e3671c9f8a106c9974fbfce21)
* Fixed symlinks to scripts for `install_veba.sh`: [d1fad03](https://github.com/jolespin/veba/commit/d1fad03b71537cc6cc0d47fee426b6610000752a)
* Added missing `CHECKM_DATA_PATH` environment variable to `VEBA-binning-prokaryotic_env` and `VEBA-classify_env`: [d1fad03](https://github.com/jolespin/veba/commit/d1fad03b71537cc6cc0d47fee426b6610000752a)

___


##### Release v1.0.1

Small patch fix:

* Fixed the fatal binning-eukaryotic.py error: [7c5addf](https://github.com/jolespin/veba/commit/7c5addf9ed6e8e45502274dd353f20b211838a41)
* Fixed the minor file naming in cluster.py: [5803845](https://github.com/jolespin/veba/commit/58038451dac0791899aa7fca3f9d79454cb9ed46)
* Removes left-over human genome tar.gz during database download/config: [5803845](https://github.com/jolespin/veba/commit/58038451dac0791899aa7fca3f9d79454cb9ed46)

___


##### Release v1.0.0
* Released with *BMC Bionformatics* publication (doi:10.1186/s12859-022-04973-8).

________________________________________________________________


#### Future Releases:

**Definitely:**

* Create a `assembly_longreads.py` module that uses `MetaFlye`
* Create a `noncoding.py` module that uses `tRNASCAN-SE` and other goodies.
* Expand Microeukaryotic Protein Database
* Automated consensus protein cluster annotations.  First need to create a hierarchical naming scheme that uses NR > KOFAM > Pfam.
* Create a wrapper around `hmmsearch` that takes in score cutoffs and outputs a useable table.  This will be used in place of `KOFAMSCAN` which creates thousands of intermediate files.

**Probably (Yes)?:**
* Add a `metabolic.py` module
* Swap [`TransDecoder`](https://github.com/TransDecoder/TransDecoder) for [`TransSuite`](https://github.com/anonconda/TranSuite)

**...Maybe (Not)?**

* Add `VAMB` as an option for `binning-prokaryotic.py` (requires python >= 3.7,<3.8)
* Add an option for sample name prefix in `assembly.py`




________________________________________________________________



#### Change Log:
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
* [2022.02.22] - Made the original `preprocess.py` -> `preprocess-kneaddata.py` and the new `preprocess.py` a wrapper around `fastq_preprocessor`
* [2022.02.22] - Made the `index.py` module
* [2022.02.22] - `concatenate_fasta.py` and `concatenate_gff.py`
* [2022.02.02] - `consensus_genome_classification.py``

________________________________________________________________

#### Next up:

* Add MAG-level counts to prokaryotic and eukaryotic. Add optional bam file for viral binning, if so then add MAG-level counts
* Add support for Anvi'o object export
* Add the [--name] prefix to all scaffolds to avoid rare situations where the contigs have the same name.
* Add spatial coverage to coverage.py script like in mapping.py script? Maybe just the samtools coverage output.


________________________________________________________________

#### Upcoming Modules:

* noncoding.py module: t-RNAscan-SE, BARRNAP, CORDON
* metabolism.py module: gapseq? Metage2Metabo?
* reassembly.py module: SPAdes

________________________________________________________________

#### Additional:
	
* Add a wrapper around `transDecoder` for metatranscriptomic assembly.


**GenoPype:**

* Get a mapping between {step:intermediate directory} and clear out intermediate directories for `--restart_from_checkpoint`.

**Classify:**

* Adapt `classify-eukaryotic.py` to handle actual genomes instead of just the results from `binning-eukaryotic.py`.  This should be straightforward, it will just need to run `metaeuk_wrapper.py` and use the output.

**Binning:**
* Make a `--skip_concoct` option because `CONCOCT` is slow when there's a lot of BAM files while `MaxBin2` is slow when there's a lot of contigs. 
* Make a wrapper that does coassembly binning but then splits out the bins into individual samples like in VAMB
* Need to replace the grep -v with non-eukaryota.list for final scaffolds_to_bins.tsv in binning-eukaryotic with subset_table.py using --inverse because sometimes it is returned as empty when list of contigs is large
* Add geNomad as main option: https://github.com/apcamargo/genomad/
* In very rare cases `identifier_mapping.tsv`from prokaryotic binning have entries where there's no MAG (i.e., columns 1 and 2 are not-null but columns 3 is null).  Identify and fix this error. 

**Index/Mapping:**

* Add STAR support. The limiting factor here is getting an analog to exon in the prodigal generated GTF for this option -sjdbGTFfeatureExon. Will need to make some adjustments to `append_geneid_to_prodigal_gff.py` but this will require new lines instead of modifying existing lines.

* Relevant GitHub issues:
https://github.com/alexdobin/STAR/issues/994
https://github.com/alexdobin/STAR/issues/867

**Cluster:**

* Add Anvi'o pangeome for each cluster



**Annotate:**

* If --identifier_mapping is `[id_orf]<tab>[id_contig]<tab>[id_mag]<tab>[id_orthogroup] then it does consensus annotations for orthogroups using UniFunc

**Scripts:**

* get_orthogroup\_consensus\_annotation.py [Not working w/ stdin]






