### What's next for *VEBA*?

*VEBA* is currently under active development. If you are interested in requesting features or wish to report a bug, please post a GitHub issue prefixed with the tag `[Feature Request]` and `[Bug]`, respectively.  If you want to contribute or have any other inquiries, contact me at `jespinoz[A|T]jcvi[DOT]org`.


________________________________________________________________

#### Current Releases:

##### Release v1.0.2
* Updated *GTDB-Tk* in `VEBA-binning-prokaryotic_env` from `1.x` to `2.x` (this version uses much less memory): [f3507dd](https://github.com/jolespin/veba/commit/f3507dd13a42960e3671c9f8a106c9974fbfce21)
* Updated the *GTDB-Tk* database from `R202` to `R207_v2` to be compatible with *GTDB-Tk v2.x*: [f3507dd](https://github.com/jolespin/veba/commit/f3507dd13a42960e3671c9f8a106c9974fbfce21)
* Updated the [GRCh38 no-alt analysis set](https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip) to [T2T CHM13v2.0](https://genome-idx.s3.amazonaws.com/bt/chm13v2.0.zip) for the default human reference: [5ccb4e2](https://github.com/jolespin/veba/commit/5ccb4e20564513707fcc3420b18237974455e196)
* Added an experimental `amplicon.py` module for short-read ASV detection via the *DADA2* workflow of *QIIME2*: [cd4ed2b](https://github.com/jolespin/veba/commit/cd4ed2bfe35d5379a63dd3294c229f2c861f6f77)
* Added additional functionality to `compile_reads_table.py` to handle advanced parsing of samples from fastq directories while also maintaining support for parsing filenames from `veba_output/preprocess`: [cd4ed2b](https://github.com/jolespin/veba/commit/cd4ed2bfe35d5379a63dd3294c229f2c861f6f77)
* Added `sra-tools` to `VEBA-preprocess_env`: [f3507dd](https://github.com/jolespin/veba/commit/f3507dd13a42960e3671c9f8a106c9974fbfce21)


##### Release v1.0.1

Small patch fix:

* Fixed the fatal binning-eukaryotic.py error: [7c5addf](https://github.com/jolespin/veba/commit/7c5addf9ed6e8e45502274dd353f20b211838a41)
* Fixed the minor file naming in cluster.py: [5803845](https://github.com/jolespin/veba/commit/58038451dac0791899aa7fca3f9d79454cb9ed46)
* Removes left-over human genome tar.gz during database download/config: [5803845](https://github.com/jolespin/veba/commit/58038451dac0791899aa7fca3f9d79454cb9ed46)

##### Release v1.0.0
* Released with *BMC Bionformatics* publication (doi:10.1186/s12859-022-04973-8).

________________________________________________________________


#### Future Releases:

##### Release v1.1.0 [In Development]

Completed:

* √ Change default human reference genome from [GRCh38 no-alt analysis set](https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip) to [T2T CHM13v2.0](https://genome-idx.s3.amazonaws.com/bt/chm13v2.0.zip)
* √ Update *GTDBTk* v1.x to v2.x and the database from [R202](https://data.gtdb.ecogenomic.org/releases/release202/202.0/auxillary_files/gtdbtk_r202_data.tar.gz) to [R207_v2](https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_v2_data.tar.gz)
* √ Add **experimental** `amplicon.py` module for ASV detection via DADA2 implemented in QIIME2. 

 
Pending:

* Add an option for sample name prefix in `assembly.py`
* Add support for *geNomad* for viral binning instead of *VirFinder*.
* [CONTINGENT] Once *CheckM2* is peer-reviewed and available on Conda, it will replace *CheckM* and the automated CPR workflow implemented by *VEBA*.


________________________________________________________________



#### Change Log:
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
* Automate feature compression ratios in cluster.py
* Switch CheckM2 [Contingent on publication and Conda release]
* Add `conda install -c bioconda sra-tools` to VEBA-preprocess_env
* Fix ClobberError in VEBA-binning-prokaryotic_env and VEBA-classify_env
* Add DADA2 pipeline as an amplicon.py module
* Add spatial coverage to coverage.py script like in mapping.py script? Maybe just the samtools coverage output.

________________________________________________________________

#### Upcoming Modules:

* Adapt DADA2_pipeline to amplicon.py module
* biosynthetic.py module: antiSMASH
* noncoding.py module that runs t-RNAscan-SE and BARRNAP (and CORDON?)
* metabolism.py module: gapseq

________________________________________________________________

#### Additional:

**Binning:**

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






