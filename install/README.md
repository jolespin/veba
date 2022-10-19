### Installation and Database Configuration Guide
____________________________________________________________
#### Software installation
One issue with having large-scale pipeline suites with open-source software is the issue of dependencies.  One solution for this is to have a modular software structure where each module has its own `conda` environment.  This allows for minimizing dependency constraints as this software suite uses an array of diverse packages from different developers. 

The basis for these environments is creating a separate environment for each module with the `VEBA-` prefix and `_env` as the suffix.  For example `VEBA-assembly_env` or `VEBA-binning-prokaryotic_env`.  Because of this, `VEBA` is currently not available as a `conda` package but each module will be in the near future.  In the meantime, please use the `install/install_veba.sh` script which installs each environment from the yaml files in `install/environments/`. After installing the environments, use the `install/download_databases` script to download and configure the databases while also adding the environment variables to the activate/deactivate scripts in each environment.  To install anything manually, just read the scripts as they are well documented and refer to different URL and paths for specific installation options.

The majority of the time taken to build database is decompressing large archives, `Diamond` database creation of NR, and `MMSEQS2` database creation of microeukaryotic protein database.

Total size is 369G but if you have certain databases installed already then you can just symlink them. 

Each major version will be packaged as a [release](https://github.com/jolespin/veba/releases) which will include a log of module and script versions. 

____________________________________________________________

**There are 3 steps to install *VEBA*:**

* Download repository from GitHub

* Install conda environments

* Download/configure databases


**1. Download repository**

```
# For stable version, download and uncompress the tarball:
wget https://github.com/jolespin/veba/archive/refs/tags/v1.0.0.tar.gz
tar -xvf v1.0.0.tar.gz && mv veba-1.0.0 veba

# For developmental version, clone the repository:
# git clone https://github.com/jolespin/veba/

# Update the permissions
chmod 755 veba/src/*.py
chmod 755 veba/src/scripts/*

# Go into the install directory
cd veba/install
``` 

**2. Install VEBA environments**

This is not resource intensive and does not require grid access but took ~1.5 hours (~90 minutes) to download and configure the dependencies.  Advanced users may speed this up by replacing `conda` with [`mamba`](https://github.com/mamba-org/mamba) in the installation script assuming `mamba` is in the base environment.  Though, `mamba` only dropped the time down by 10 minutes.

```
bash install_veba.sh
```

**3. Activate the database conda environment, download, and configure databases**

⚠️ This step takes ~4.5 hrs using 8 threads with 128G memory and should be run using a compute grid via SLURM or SunGridEngine.  If this command is run on the head node it will likely fail or timeout if a connection is interrupted. The most computationally intensive steps are creating a `Diamond` database of NCBI's non-redundant reference and a `MMSEQS2` database of the microeukaryotic protein database.

If issues arise, please [submit a GitHub issue](https://github.com/jolespin/veba/issues) prefixed with [DATABASE]. We are here to help :)

**If you are running an interactive queue:**

```
conda activate VEBA-database_env

bash download_databases.sh /path/to/veba_database
```

**If you use job scheduling (e.g., sbatch or qsub):**

[If you're unfamiliar with SLURM or SunGridEnginer, we got you](https://github.com/jolespin/veba/blob/main/walkthroughs/README.md#basics).  

Running `conda activate` on a compute server might prompt you to run `conda init` even if you've already initilized on the head node.  To get around this you can use `source activate [environment]`.  Using the `source activate` command requires you to be in `base` conda environment.  You can do this via `conda deactivate` or `conda activate base` before you submit your job.  

Below is an example of how to do this: 

```
# Activate your base environment
conda activate base

# Set the number of threads you want to use.  
# Keep in mind that not all steps are parallelized (e.g., wget)
N_JOBS=8 

# Create a log directory
mkdir -p logs/

# Set name for log files when running on the grid
N="database_config"

# Adapt your command to a one-liner
CMD="source activate VEBA-database_env && bash download_databases.sh /path/to/veba_database"
	
# Note: You should either use SunGridEngine or SLURM not both. 
# You might need to adapt slightly but these worked on our systems.

# SunGridEngine:
qsub -o logs/${N}.o -e logs/${N}.e -cwd -N ${N} -j y -pe threaded ${N_JOBS} "${CMD}"
	
# SLURM:
# For SLURM you might need to specify which account and partition you are associated with for grid jobs
PARTITION=[partition name]
ACCOUNT=[account name]

sbatch -A ${ACCOUNT} -p ${PARTITION} -J ${N} -N 1 -c ${N_JOBS} --ntasks-per-node=1 -o logs/${N}.o -e logs/${N}.e --export=ALL -t 12:00:00 --mem=128G --wrap="${CMD}"
```

You should be done now. If you want to double check the installation and database configuration worked then activate a `VEBA` environment: 

```
VEBA-annotate_env
VEBA-assembly_env
VEBA-binning-eukaryotic_env
VEBA-binning-prokaryotic_env
VEBA-binning-viral_env
VEBA-classify_env
VEBA-cluster_env
VEBA-database_env
VEBA-mapping_env
VEBA-phylogeny_env
VEBA-preprocess_env
```
and check that `VEBA_DATABASE` environment variable is set. If not, then add it manually to ~/.bash_profile: `export VEBA_DATABASE=/path/to/veba_database`.

Future versions will have `bioconda` installation available.

#### Common installation errors that do not affect VEBA functionality:

You may get the following **non-fatal** errors but you can ignore these: 

* *SafetyError* from `CheckM v1.1.3`

```
SafetyError: The package for checkm-genome located at /path/to/anaconda3/pkgs/checkm-genome-1.1.3-py_1
appears to be corrupted. The path 'site-packages/checkm/DATA_CONFIG'
has an incorrect size.
  reported size: 215 bytes
  actual size: 251 bytes
```

* *ClobberError* from `Perl` packages
 
```
ClobberError: This transaction has incompatible packages due to a shared path.
  packages: bioconda/linux-64::perl-params-util-1.07-pl526h6bb024c_4, bioconda/linux-64::perl-package-stash-0.38-pl526hf484d3e_1
  path: 'lib/site_perl/5.26.2/x86_64-linux-thread-multi/auto/Params/Util/.packlist'
```

These errors have been reported during the creation of `VEBA-binning-prokaryotic_env` and `VEBA-classify_env` both of which are the only environments that install `CheckM`.  However, the SafetyError [is easily reproducible](https://github.com/Ecogenomics/CheckM/issues/349) and similar ClobberErrors have been [reported in other software suites](https://github.com/bxlab/metaWRAP/issues/72) which are likely due to Perl dependencies.  Working on addressing these minor errors but they are non-fatal so not the highest priority.  
____________________________________________________________


#### Database Structure:

**Current:**

```
tree -L 3 .
.
├── ACCESS_DATE
├── Annotate
│   ├── KOFAM
│   │   ├── ko_list
│   │   └── profiles
│   ├── nr
│   │   └── nr.dmnd
│   └── Pfam
│       └── Pfam-A.hmm.gz
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
│   ├── chm13v2.0
│   │   ├── chm13v2.0.1.bt2
│   │   ├── chm13v2.0.2.bt2
│   │   ├── chm13v2.0.3.bt2
│   │   ├── chm13v2.0.4.bt2
│   │   ├── chm13v2.0.rev.1.bt2
│   │   └── chm13v2.0.rev.2.bt2
│   └── kmers
│       └── ribokmers.fa.gz
├── MarkerSets
│   ├── Archaea_76.hmm
│   ├── Bacteria_71.hmm
│   ├── CPR_43.hmm
│   ├── eukaryota_odb10.hmm
│   ├── eukaryota_odb10.scores_cutoff.tsv.gz
│   ├── Fungi_593.hmm
│   ├── Protista_83.hmm
│   └── README
└── SIZE

33 directories, 48 files
```

**Previous:**
<details>
	<summary>v1.0.0 Database</summary>
	
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
│       └── Pfam-A.hmm.gz
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

</details>
____________________________________________________________

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

Espinoza, Josh (2022): Profile HMM marker sets. figshare. Dataset. https://doi.org/10.6084/m9.figshare.19616016.v1 

#### Microeukaryotic protein database:
A protein database is required not only for eukaryotic gene calls using MetaEuk but can also be used for MAG annotation.  Many eukaryotic protein databases exist such as MMETSP, EukZoo, and EukProt, yet these are limited to marine environments, include prokaryotic sequences, or include eukaryotic sequences for organisms that would not be expected to be binned out of metagenomes such as metazoans.  We combined and dereplicated MMETSP, EukZoo, EukProt, and NCBI non-redundant to include only microeukaryotes such as protists and fungi.  This optimized microeukaryotic database ensures that only eukaryotic exons expected to be represented in metagenomes are utilized for eukaryotic gene modeling and the resulting MetaEuk reference targets are used for eukaryotic MAG classification.  VEBA’s microeukaryotic protein database includes 48,006,918 proteins from 42,922 microeukaryotic strains.  

Espinoza, Josh (2022): Microeukaryotic Protein Database. figshare. Dataset. https://doi.org/10.6084/m9.figshare.19668855.v1 

____________________________________________________________

#### Version Notes:

* The `nr.dmnd` should be built using NCBI's taxonomy info.  The `download_database.sh` takes care of this but if you are using a prebuilt `nr.dmnd` database then use following command for reference: `diamond makedb --in ${DATABASE_DIRECTORY}/nr.gz --db ${DATABASE_DIRECTORY}/Annotate/nr.dmnd --taxonmap ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/prot.accession2taxid.FULL.gz --taxonnodes ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/nodes.dmp --taxonnames ${DATABASE_DIRECTORY}/Classify/NCBITaxonomy/names.dmp`
* The NCBI taxonomy should be downloaded on the same date as NR to make sure the identifiers match up between datasets.  NCBI deprecates taxonomy identifiers and adds new ones between versions so downloading on the same day should minimize that discrepancy.  One caveat NCBI NR and taxonomy databases is the versioning is difficult to discern.
* For the human contamination, if you use `KneadData` and already have a `Bowtie2` index for human then you can use that instead.  The only module that uses this is `preprocess.py` and you have to specify this directly when running (i.e., it's optional) so it doesn't matter if it's in the database directory or not (same with ribokmers.fa.gz). 
* `CheckM` and `CheckV` only have 1 database version at this time so it isn't an issue. 
* `KOFAM` and `Pfam` just uses these as annotations so any version should work perfectly.
* Again, if you are low on disk space and already have these installed then just symlink them with the structure above. If so, them just comment out those sections of `download_databases.sh`.  

**VEBA v1.0.0 Specific:**

* The only mandatory datbase version is `R202` for `GTDBTk` because `VEBA` is built on `v1.5.0 ≤ GTDBTk ≤ v1.70`.  The next major `VEBA` update will upgrade to `≥ GTDBTk v2.0.0` and the associated `R207_v2` database.

If you have any issues, please create a GitHub issue with the prefix `[DATABASE]` followed by the question.