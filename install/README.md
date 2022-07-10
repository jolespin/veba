### Installation and Database Configuration Guide
____________________________________________________________
#### Software installation
One issue with having large-scale pipeline suites with open-source software is the issue of dependencies.  One solution for this is to have a modular software structure where each module has its own `conda` environment.  This allows for minimizing dependency constraints as this software suite uses an array of diverse packages from different developers. 

The basis for these environments is creating a separate environment for each module with the `VEBA-` prefix and `_env` as the suffix.  For example `VEBA-assembly_env` or `VEBA-binning-prokaryotic_env`.  Because of this, `VEBA` is currently not available as a `conda` package but each module will be in the near future.  In the meantime, please use the `install/install_veba.sh` script which installs each environment from the yaml files in `install/environments/` and then installs the database. To install the database separately, use the `install/download_databases` script.  To install anything manually, just read the scripts as they are well documented and refer to different URL and paths for specific installation options.

The majority of the time taken to build database is decompressing large archives, Diamond database creation of NR, and MMSEQS2 database creation of microeukaryotic protein database.

Total size is 369G but if you have certain databases installed already then you can just symlink them. 



```
# Download repository
git clone https://github.com/jolespin/veba/
cd veba/install
 
# Install VEBA environments 
bash install_veba.sh

# Activate that database conda environment
conda activate VEBA-database_env

# Download databases (This takes ~9.5 hrs using 8 threads with 128G memory)
bash download_databases.sh /path/to/veba_database

# Environment variable
# Check that `VEBA_DATABASE` environment variable is set. If not, then add it manually to ~/.bash_profile: `export VEBA_DATABASE=/path/to/veba_database`

```

#### Database Structure:
```
tree -L 3 .
.
.
├── ACCESS_DATE
├── Annotate
│   ├── KOFAM
│   │   ├── ko_list
│   │   └── profiles
│   ├── nr
│   │   └── nr.dmnd
│   └── Pfam
│       └── Pfam-A.full.gz
├── Classify
│   ├── CheckM
│   │   ├── distributions
│   │   ├── genome_tree
│   │   ├── hmms
│   │   ├── hmms_ssu
│   │   ├── img
│   │   ├── pfam
│   │   ├── selected_marker_sets.tsv
│   │   ├── taxon_marker_sets.tsv
│   │   └── test_data
│   ├── CheckV
│   │   ├── genome_db
│   │   ├── hmm_db
│   │   └── README.txt
│   ├── GTDBTk
│   │   ├── fastani
│   │   ├── manifest.tsv
│   │   ├── markers
│   │   ├── masks
│   │   ├── metadata
│   │   ├── mrca_red
│   │   ├── msa
│   │   ├── pplacer
│   │   ├── radii
│   │   └── taxonomy
│   ├── Microeukaryotic
│   │   ├── microeukaryotic
│   │   ├── microeukaryotic.dbtype
│   │   ├── microeukaryotic_h
│   │   ├── microeukaryotic_h.dbtype
│   │   ├── microeukaryotic_h.index
│   │   ├── microeukaryotic.index
│   │   ├── microeukaryotic.lookup
│   │   ├── microeukaryotic.source
│   │   ├── reference.rmdup.iupac.relabeled.no_deprecated.complete_lineage.faa.gz
│   │   ├── source_taxonomy.tsv.gz
│   │   ├── source_to_lineage.dict.pkl.gz
│   │   └── target_to_source.dict.pkl.gz
│   └── NCBITaxonomy
│       ├── citations.dmp
│       ├── delnodes.dmp
│       ├── division.dmp
│       ├── gc.prt
│       ├── gencode.dmp
│       ├── merged.dmp
│       ├── names.dmp
│       ├── nodes.dmp
│       ├── prot.accession2taxid.FULL.gz
│       ├── readme.txt
│       ├── taxa.sqlite
│       └── taxa.sqlite.traverse.pkl
├── Contamination
│   ├── grch38
│   │   ├── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.1.bt2
│   │   ├── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.2.bt2
│   │   ├── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.3.bt2
│   │   ├── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.4.bt2
│   │   ├── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.rev.1.bt2
│   │   └── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.rev.2.bt2
│   └── kmers
│       └── ribokmers.fa.gz
└── MarkerSets
    ├── Archaea_76.hmm
    ├── Bacteria_71.hmm
    ├── CPR_43.hmm
    ├── eukaryota_odb10.hmm
    ├── eukaryota_odb10.scores_cutoff.tsv.gz
    ├── Fungi_593.hmm
    ├── Protista_83.hmm
    └── README

33 directories, 47 files
```

#### Profile HMM Sources:
Please cite the following sources if these marker sets are used in any way:

```
* Archaea_76.hmm - (Anvi'o) Lee, https://doi.org/10.1093/bioinformatics/btz188 (https://github.com/merenlab/anvio/tree/master/anvio/data/hmm/Archaea_76)

* Bacteria_71.hmm - (Anvi'o) Lee modified, https://doi.org/10.1093/bioinformatics/btz188 (https://github.com/merenlab/anvio/tree/master/anvio/data/hmm/Bacteria_71)

* Protista_83.hmm - (Anvi'o) Delmont, http://merenlab.org/delmont-euk-scgs (https://github.com/merenlab/anvio/tree/master/anvio/data/hmm/Protista_83)

* Fungi_593.hmm - (FGMP) https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2782-9

* CPR_43.hmm - (CheckM) https://github.com/Ecogenomics/CheckM/tree/master/custom_marker_sets

* eukaryota_odb10 - (BUSCO) https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2020-09-10.tar.gz
```


#### Version Notes:
* The only mandatory datbase version is `R202` for `GTDBTk` because `VEBA` is built on `v1.5.0 ≤ GTDBTk ≤ v1.70`.  The next major `VEBA` update will upgrade to `GTDBTk v2.1.0` and the associated `R207_v2` database.
* The `nr.dmnd` should be built using NCBI's taxonomy info.  The `download_database.sh` takes care of this but if you are using a prebuilt `nr.dmnd` database then use following command for reference: `diamond makedb --in ${DATABASE_DIRECTORY}/nr.gz --db ${DATABASE_DIRECTORY}/Annotate/nr.dmnd --taxonmap ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/prot.accession2taxid.FULL.gz --taxonnodes ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/nodes.dmp --taxonnames ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/names.dmp`
* The NCBI taxonomy should be downloaded on the same date as NR to make sure the identifiers match up between datasets.  NCBI deprecates taxonomy identifiers and adds new ones between versions so downloading on the same day should minimize that discrepancy.  One caveat NCBI NR and taxonomy databases is the versioning is difficult to discern.
* For the human contamination, if you use `KneadData` and already have a `Bowtie2` index for human then you can use that instead.  The only module that uses this is `preprocess.py` and you have to specify this directly when running (i.e., it's optional) so it doesn't matter if it's in the database directory or not (same with ribokmers.fa.gz). 
* `CheckM` and `CheckV` only have 1 database version at this time so it isn't an issue. 
* `KOFAM` and `Pfam` just uses these as annotations so any version should work perfectly.
* Again, if you are low on disk space and already have these installed then just symlink them with the structure above. If so, them just comment out those sections of `download_databases.sh`.  

If you have any issues, please create a GitHub issue with the prefix `[DATABASE]` followed by the question.