### VEBA Database:

#### Profile HMM Sources:
Please cite the following sources if these marker sets are used in any way:

```
* Archaea_76.hmm.gz - (Anvi'o) Lee, https://doi.org/10.1093/bioinformatics/btz188 (https://github.com/merenlab/anvio/tree/master/anvio/data/hmm/Archaea_76)

* Bacteria_71.hmm.gz - (Anvi'o) Lee modified, https://doi.org/10.1093/bioinformatics/btz188 (https://github.com/merenlab/anvio/tree/master/anvio/data/hmm/Bacteria_71)

* Protista_83.hmm.gz - (Anvi'o) Delmont, http://merenlab.org/delmont-euk-scgs (https://github.com/merenlab/anvio/tree/master/anvio/data/hmm/Protista_83)

* Fungi_593.hmm.gz - (FGMP) https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2782-9

* CPR_43.hmm.gz - (CheckM) https://github.com/Ecogenomics/CheckM/tree/master/custom_marker_sets

* eukaryota_odb10.hmm.gz - (BUSCO) https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2020-09-10.tar.gz
```

Espinoza, Josh (2022): Profile HMM marker sets. figshare. Dataset. https://doi.org/10.6084/m9.figshare.19616016.v1 

#### Microeukaryotic protein database:
VEBA’s Microeukaryotic Protein Database has been completely redesigned using the logic of UniRef and their clustered database.  The previous microeukaryotic protein database contained 48,006,918 proteins from 44,647 source organisms while the updated database, MicroEuk, contains 79,920,430 proteins from 52,495 source organisms.  As in the prior major release, MicroEuk concentrates on microeukaryotic organisms while excluding higher eukaryotes as these organisms are the primary eukaryotes targeted by shotgun metagenomics and metatranscriptomics.  Source organisms in this context are defined as organisms in which the proteins were derived.  

**Number of sequences:**
 * MicroEuk100 = 79,920,431 (19 GB)
 * MicroEuk90  = 51,767,730 (13 GB)
 * MicroEuk50  = 29,898,853 (6.5 GB)

**Number of source organisms per dataset:**

* MycoCosm = 2503
* PhycoCosm = 174
* EnsemblProtists = 233
* MMETSP = 759
* TARA_SAGv1 = 8
* EukProt = 366
* EukZoo = 27
* TARA_SMAGv1 = 389
* NR_Protists-Fungi = 48217

**Current:**

* [MicroEuk\_v3](https://zenodo.org/records/10139451) available on Zenodo

**Deprecated:**

* [MicroEuk\_v2](https://zenodo.org/record/7485114) available on Zenodo

* [MicroEuk\_v1](https://figshare.com/articles/dataset/Microeukaryotic_Protein_Database/19668855) available on FigShare

#### Database Structure:

**Current:**

 *VEBA Database* version: `VDB_v8` (~300 GB)
 
* Replaced [MicrobeAnnotator-KEGG](https://zenodo.org/records/10020074) reimplementation with [KEGG Pathway Profiler](https://github.com/jolespin/kegg_pathway_profiler) database

```
├── ACCESS_DATE
├── Annotate
│   ├── CAZy
│   │   └── CAZyDB.07262023.dmnd
│   ├── KEGG-Pathway-Profiler
│   │   ├── database.pkl.gz
│   │   ├── database.tsv
│   │   ├── database.version
│   │   ├── pathway_names.tsv.gz
│   │   └── pathways
│   ├── KOfam
│   │   ├── ko_list
│   │   └── profiles
│   ├── MIBiG
│   │   └── mibig_v3.1.dmnd
│   ├── NCBIfam-AMRFinder
│   │   ├── NCBIfam-AMRFinder.changelog.txt
│   │   ├── NCBIfam-AMRFinder.changelog.txt.1
│   │   ├── NCBIfam-AMRFinder.hmm.gz
│   │   ├── NCBIfam-AMRFinder.tsv
│   │   └── NCBIfam-AMRFinder.tsv.1
│   ├── Pfam
│   │   ├── Pfam-A.clans.tsv.gz
│   │   ├── Pfam-A.clans.tsv.gz.1
│   │   ├── Pfam-A.hmm.gz
│   │   ├── Pfam-A.hmm.gz.1
│   │   ├── relnotes.txt
│   │   └── relnotes.txt.1
│   ├── UniRef
│   │   ├── uniref50.dmnd
│   │   ├── uniref50.release_note
│   │   ├── uniref90.dmnd
│   │   ├── uniref90.release_note
│   │   └── uniref90.release_note.1
│   └── VFDB
│       ├── VFDB_setA_pro.dmnd
│       └── VFs.xls.gz
├── Classify
│   ├── CheckM2
│   │   └── uniref100.KO.1.dmnd
│   ├── CheckV
│   │   ├── database_version
│   │   ├── genome_db
│   │   ├── hmm_db
│   │   └── README.txt
│   ├── geNomad
│   │   ├── genomad_db
│   │   ├── genomad_db.dbtype
│   │   ├── genomad_db_h
│   │   ├── genomad_db_h.dbtype
│   │   ├── genomad_db_h.index
│   │   ├── genomad_db.index
│   │   ├── genomad_db.lookup
│   │   ├── genomad_db_mapping
│   │   ├── genomad_db.source
│   │   ├── genomad_db_taxonomy
│   │   ├── genomad_integrase_db
│   │   ├── genomad_integrase_db.dbtype
│   │   ├── genomad_integrase_db_h
│   │   ├── genomad_integrase_db_h.dbtype
│   │   ├── genomad_integrase_db_h.index
│   │   ├── genomad_integrase_db.index
│   │   ├── genomad_integrase_db.lookup
│   │   ├── genomad_integrase_db.source
│   │   ├── genomad_marker_metadata.tsv
│   │   ├── genomad_mini_db
│   │   ├── genomad_mini_db.dbtype
│   │   ├── genomad_mini_db_h
│   │   ├── genomad_mini_db_h.dbtype
│   │   ├── genomad_mini_db_h.index
│   │   ├── genomad_mini_db.index
│   │   ├── genomad_mini_db.lookup
│   │   ├── genomad_mini_db_mapping
│   │   ├── genomad_mini_db.source
│   │   ├── genomad_mini_db_taxonomy
│   │   ├── mini_set_ids
│   │   ├── names.dmp
│   │   ├── nodes.dmp
│   │   ├── plasmid_hallmark_annotation.txt
│   │   ├── version.txt
│   │   └── virus_hallmark_annotation.txt
│   ├── GTDB
│   │   ├── database_version
│   │   ├── markers
│   │   ├── mash
│   │   ├── masks
│   │   ├── metadata
│   │   ├── mrca_red
│   │   ├── msa
│   │   ├── pplacer
│   │   ├── radii
│   │   ├── RELEASE_NOTES.txt
│   │   ├── skani
│   │   ├── split
│   │   └── taxonomy
│   ├── MicroEuk
│   │   ├── MicroEuk100
│   │   ├── MicroEuk100.dbtype
│   │   ├── MicroEuk100.eukaryota_odb10
│   │   ├── MicroEuk100.eukaryota_odb10.dbtype
│   │   ├── MicroEuk100.eukaryota_odb10_h
│   │   ├── MicroEuk100.eukaryota_odb10_h.dbtype
│   │   ├── MicroEuk100.eukaryota_odb10_h.index
│   │   ├── MicroEuk100.eukaryota_odb10.index
│   │   ├── MicroEuk100.eukaryota_odb10.lookup
│   │   ├── MicroEuk100.eukaryota_odb10.source
│   │   ├── MicroEuk100_h
│   │   ├── MicroEuk100_h.dbtype
│   │   ├── MicroEuk100_h.index
│   │   ├── MicroEuk100.index
│   │   ├── MicroEuk100.lookup
│   │   ├── MicroEuk100.source
│   │   ├── MicroEuk50
│   │   ├── MicroEuk50.dbtype
│   │   ├── MicroEuk50_h
│   │   ├── MicroEuk50_h.dbtype
│   │   ├── MicroEuk50_h.index
│   │   ├── MicroEuk50.index
│   │   ├── MicroEuk50.lookup
│   │   ├── MicroEuk50.source
│   │   ├── MicroEuk90
│   │   ├── MicroEuk90.dbtype
│   │   ├── MicroEuk90_h
│   │   ├── MicroEuk90_h.dbtype
│   │   ├── MicroEuk90_h.index
│   │   ├── MicroEuk90.index
│   │   ├── MicroEuk90.lookup
│   │   ├── MicroEuk90.source
│   │   ├── source_taxonomy.tsv.gz
│   │   ├── source_to_lineage.dict.pkl.gz
│   │   └── target_to_source.dict.pkl.gz
│   └── NCBITaxonomy
│       ├── citations.dmp
│       ├── delnodes.dmp
│       ├── division.dmp
│       ├── gc.prt
│       ├── gencode.dmp
│       ├── images.dmp
│       ├── merged.dmp
│       ├── names.dmp
│       ├── nodes.dmp
│       └── readme.txt
├── Contamination
│   ├── AntiFam
│   │   ├── AntiFam.hmm.gz
│   │   ├── relnotes
│   │   └── version
│   ├── chm13v2.0
│   │   ├── chm13v2.0.1.bt2
│   │   ├── chm13v2.0.2.bt2
│   │   ├── chm13v2.0.3.bt2
│   │   ├── chm13v2.0.4.bt2
│   │   ├── chm13v2.0.rev.1.bt2
│   │   └── chm13v2.0.rev.2.bt2
│   └── kmers
│       └── ribokmers.fa.gz
├── MarkerSets
│   ├── Archaea_76.hmm.gz
│   ├── Bacteria_71.hmm.gz
│   ├── CPR_43.hmm.gz
│   ├── eukaryota_odb10.hmm.gz
│   ├── eukaryota_odb10.scores_cutoff.tsv.gz
│   ├── Fungi_593.hmm.gz
│   ├── Protista_83.hmm.gz
│   └── README
├── VDB_v8.db_size
└── VERSION
```

**Current:**

 *VEBA Database* version: `VDB_v7` (301 GB)
 
* Updated GTDB r214.1 -> r220
* Added Pfam clans
* Changed KOFAM -> KOfam

```
tree -L 3 .
.
|-- ACCESS_DATE
|-- Annotate
|   |-- CAZy
|   |   `-- CAZyDB.07262023.dmnd
|   |-- KOfam
|   |   |-- ko_list
|   |   `-- profiles
|   |-- MIBiG
|   |   `-- mibig_v3.1.dmnd
|   |-- MicrobeAnnotator-KEGG
|   |   |-- KEGG_Bifurcating_Module_Information.pkl
|   |   |-- KEGG_Bifurcating_Module_Information.pkl.md5
|   |   |-- KEGG_Module_Information.txt
|   |   |-- KEGG_Module_Information.txt.md5
|   |   |-- KEGG_Regular_Module_Information.pkl
|   |   |-- KEGG_Regular_Module_Information.pkl.md5
|   |   |-- KEGG_Structural_Module_Information.pkl
|   |   `-- KEGG_Structural_Module_Information.pkl.md5
|   |-- NCBIfam-AMRFinder
|   |   |-- NCBIfam-AMRFinder.changelog.txt
|   |   |-- NCBIfam-AMRFinder.hmm.gz
|   |   `-- NCBIfam-AMRFinder.tsv
|   |-- Pfam
|   |   |-- Pfam-A.clans.tsv.gz
|   |   |-- Pfam-A.hmm.gz
|   |   `-- relnotes.txt
|   |-- UniRef
|   |   |-- uniref50.dmnd
|   |   |-- uniref50.release_note
|   |   |-- uniref90.dmnd
|   |   `-- uniref90.release_note
|   `-- VFDB
|       |-- VFDB_setA_pro.dmnd
|       `-- VFs.xls.gz
|-- Classify
|   |-- CheckM2
|   |   `-- uniref100.KO.1.dmnd
|   |-- CheckV
|   |   |-- README.txt
|   |   |-- database_version
|   |   |-- genome_db
|   |   `-- hmm_db
|   |-- GTDB
|   |   |-- RELEASE_NOTES.txt
|   |   |-- database_version
|   |   |-- markers
|   |   |-- mash
|   |   |-- masks
|   |   |-- metadata
|   |   |-- mrca_red
|   |   |-- msa
|   |   |-- pplacer
|   |   |-- radii
|   |   |-- skani
|   |   |-- split
|   |   `-- taxonomy
|   |-- MicroEuk
|   |   |-- MicroEuk100
|   |   |-- MicroEuk100.dbtype
|   |   |-- MicroEuk100.eukaryota_odb10
|   |   |-- MicroEuk100.eukaryota_odb10.dbtype
|   |   |-- MicroEuk100.eukaryota_odb10.index
|   |   |-- MicroEuk100.eukaryota_odb10.lookup
|   |   |-- MicroEuk100.eukaryota_odb10.source
|   |   |-- MicroEuk100.eukaryota_odb10_h
|   |   |-- MicroEuk100.eukaryota_odb10_h.dbtype
|   |   |-- MicroEuk100.eukaryota_odb10_h.index
|   |   |-- MicroEuk100.index
|   |   |-- MicroEuk100.lookup
|   |   |-- MicroEuk100.source
|   |   |-- MicroEuk100_h
|   |   |-- MicroEuk100_h.dbtype
|   |   |-- MicroEuk100_h.index
|   |   |-- MicroEuk50
|   |   |-- MicroEuk50.dbtype
|   |   |-- MicroEuk50.index
|   |   |-- MicroEuk50.lookup
|   |   |-- MicroEuk50.source
|   |   |-- MicroEuk50_h
|   |   |-- MicroEuk50_h.dbtype
|   |   |-- MicroEuk50_h.index
|   |   |-- MicroEuk90
|   |   |-- MicroEuk90.dbtype
|   |   |-- MicroEuk90.index
|   |   |-- MicroEuk90.lookup
|   |   |-- MicroEuk90.source
|   |   |-- MicroEuk90_h
|   |   |-- MicroEuk90_h.dbtype
|   |   |-- MicroEuk90_h.index
|   |   |-- source_taxonomy.tsv.gz
|   |   |-- source_to_lineage.dict.pkl.gz
|   |   `-- target_to_source.dict.pkl.gz
|   |-- NCBITaxonomy
|   |   |-- citations.dmp
|   |   |-- delnodes.dmp
|   |   |-- division.dmp
|   |   |-- gc.prt
|   |   |-- gencode.dmp
|   |   |-- images.dmp
|   |   |-- merged.dmp
|   |   |-- names.dmp
|   |   |-- nodes.dmp
|   |   `-- readme.txt
|   `-- geNomad
|       |-- genomad_db
|       |-- genomad_db.dbtype
|       |-- genomad_db.index
|       |-- genomad_db.lookup
|       |-- genomad_db.source
|       |-- genomad_db_h
|       |-- genomad_db_h.dbtype
|       |-- genomad_db_h.index
|       |-- genomad_db_mapping
|       |-- genomad_db_taxonomy
|       |-- genomad_integrase_db
|       |-- genomad_integrase_db.dbtype
|       |-- genomad_integrase_db.index
|       |-- genomad_integrase_db.lookup
|       |-- genomad_integrase_db.source
|       |-- genomad_integrase_db_h
|       |-- genomad_integrase_db_h.dbtype
|       |-- genomad_integrase_db_h.index
|       |-- genomad_marker_metadata.tsv
|       |-- genomad_mini_db -> genomad_db
|       |-- genomad_mini_db.dbtype
|       |-- genomad_mini_db.index
|       |-- genomad_mini_db.lookup -> genomad_db.lookup
|       |-- genomad_mini_db.source -> genomad_db.source
|       |-- genomad_mini_db_h -> genomad_db_h
|       |-- genomad_mini_db_h.dbtype -> genomad_db_h.dbtype
|       |-- genomad_mini_db_h.index -> genomad_db_h.index
|       |-- genomad_mini_db_mapping -> genomad_db_mapping
|       |-- genomad_mini_db_taxonomy -> genomad_db_taxonomy
|       |-- mini_set_ids
|       |-- names.dmp
|       |-- nodes.dmp
|       |-- plasmid_hallmark_annotation.txt
|       |-- version.txt
|       `-- virus_hallmark_annotation.txt
|-- Contamination
|   |-- AntiFam
|   |   |-- AntiFam.hmm.gz
|   |   |-- relnotes
|   |   `-- version
|   |-- chm13v2.0
|   |   |-- chm13v2.0.1.bt2
|   |   |-- chm13v2.0.2.bt2
|   |   |-- chm13v2.0.3.bt2
|   |   |-- chm13v2.0.4.bt2
|   |   |-- chm13v2.0.rev.1.bt2
|   |   `-- chm13v2.0.rev.2.bt2
|   `-- kmers
|       `-- ribokmers.fa.gz
`-- MarkerSets
    |-- Archaea_76.hmm.gz
    |-- Bacteria_71.hmm.gz
    |-- CPR_43.hmm.gz
    |-- Fungi_593.hmm.gz
    |-- Protista_83.hmm.gz
    |-- README
    |-- eukaryota_odb10.hmm.gz
    `-- eukaryota_odb10.scores_cutoff.tsv.gz

35 directories, 127 files
```


**Deprecated:**

<details>
	<summary> VEBA Database version: VDB_v6 (272 GB) </summary>

 * Added `MicroEuk_v3`

```
 tree -L 3 .
.
├── ACCESS_DATE
├── Annotate
│   ├── CAZy
│   │   └── CAZyDB.07262023.dmnd
│   ├── KOFAM
│   │   ├── ko_list
│   │   └── profiles
│   ├── MIBiG
│   │   └── mibig_v3.1.dmnd
│   ├── MicrobeAnnotator-KEGG
│   │   ├── KEGG_Bifurcating_Module_Information.pkl
│   │   ├── KEGG_Bifurcating_Module_Information.pkl.md5
│   │   ├── KEGG_Module_Information.txt
│   │   ├── KEGG_Module_Information.txt.md5
│   │   ├── KEGG_Regular_Module_Information.pkl
│   │   ├── KEGG_Regular_Module_Information.pkl.md5
│   │   ├── KEGG_Structural_Module_Information.pkl
│   │   └── KEGG_Structural_Module_Information.pkl.md5
│   ├── MicrobeAnnotator-KEGG.tar.gz
│   ├── NCBIfam-AMRFinder
│   │   ├── NCBIfam-AMRFinder.changelog.txt
│   │   ├── NCBIfam-AMRFinder.hmm.gz
│   │   └── NCBIfam-AMRFinder.tsv
│   ├── Pfam
│   │   ├── Pfam-A.hmm.gz
│   │   └── relnotes.txt
│   ├── UniRef
│   │   ├── uniref50.dmnd
│   │   ├── uniref50.release_note
│   │   ├── uniref90.dmnd
│   │   └── uniref90.release_note
│   └── VFDB
│       └── VFDB_setA_pro.dmnd
├── Classify
│   ├── CheckM2
│   │   └── uniref100.KO.1.dmnd
│   ├── CheckV
│   │   ├── genome_db
│   │   ├── hmm_db
│   │   └── README.txt
│   ├── geNomad
│   │   ├── genomad_db
│   │   ├── genomad_db.dbtype
│   │   ├── genomad_db_h
│   │   ├── genomad_db_h.dbtype
│   │   ├── genomad_db_h.index
│   │   ├── genomad_db.index
│   │   ├── genomad_db.lookup
│   │   ├── genomad_db_mapping
│   │   ├── genomad_db.source
│   │   ├── genomad_db_taxonomy
│   │   ├── genomad_integrase_db
│   │   ├── genomad_integrase_db.dbtype
│   │   ├── genomad_integrase_db_h
│   │   ├── genomad_integrase_db_h.dbtype
│   │   ├── genomad_integrase_db_h.index
│   │   ├── genomad_integrase_db.index
│   │   ├── genomad_integrase_db.lookup
│   │   ├── genomad_integrase_db.source
│   │   ├── genomad_marker_metadata.tsv
│   │   ├── genomad_mini_db -> genomad_db
│   │   ├── genomad_mini_db.dbtype
│   │   ├── genomad_mini_db_h -> genomad_db_h
│   │   ├── genomad_mini_db_h.dbtype -> genomad_db_h.dbtype
│   │   ├── genomad_mini_db_h.index -> genomad_db_h.index
│   │   ├── genomad_mini_db.index
│   │   ├── genomad_mini_db.lookup -> genomad_db.lookup
│   │   ├── genomad_mini_db_mapping -> genomad_db_mapping
│   │   ├── genomad_mini_db.source -> genomad_db.source
│   │   ├── genomad_mini_db_taxonomy -> genomad_db_taxonomy
│   │   ├── mini_set_ids
│   │   ├── names.dmp
│   │   ├── nodes.dmp
│   │   ├── plasmid_hallmark_annotation.txt
│   │   ├── version.txt
│   │   └── virus_hallmark_annotation.txt
│   ├── GTDB
│   │   ├── fastani
│   │   ├── markers
│   │   ├── mash
│   │   ├── masks
│   │   ├── metadata
│   │   ├── mrca_red
│   │   ├── msa
│   │   ├── pplacer
│   │   ├── radii
│   │   ├── split
│   │   ├── taxonomy
│   │   └── temp
│   ├── MicroEuk
│   │   ├── MicroEuk100
│   │   ├── MicroEuk100.dbtype
│   │   ├── MicroEuk100.eukaryota_odb10
│   │   ├── MicroEuk100.eukaryota_odb10.dbtype
│   │   ├── MicroEuk100.eukaryota_odb10_h
│   │   ├── MicroEuk100.eukaryota_odb10_h.dbtype
│   │   ├── MicroEuk100.eukaryota_odb10_h.index
│   │   ├── MicroEuk100.eukaryota_odb10.index
│   │   ├── MicroEuk100.eukaryota_odb10.lookup
│   │   ├── MicroEuk100.eukaryota_odb10.source
│   │   ├── MicroEuk100_h
│   │   ├── MicroEuk100_h.dbtype
│   │   ├── MicroEuk100_h.index
│   │   ├── MicroEuk100.index
│   │   ├── MicroEuk100.lookup
│   │   ├── MicroEuk100_mapping
│   │   ├── MicroEuk100.source
│   │   ├── MicroEuk100_taxonomy
│   │   ├── MicroEuk50
│   │   ├── MicroEuk50.dbtype
│   │   ├── MicroEuk50_h
│   │   ├── MicroEuk50_h.dbtype
│   │   ├── MicroEuk50_h.index
│   │   ├── MicroEuk50.index
│   │   ├── MicroEuk50.lookup
│   │   ├── MicroEuk50.source
│   │   ├── MicroEuk90
│   │   ├── MicroEuk90.dbtype
│   │   ├── MicroEuk90_h
│   │   ├── MicroEuk90_h.dbtype
│   │   ├── MicroEuk90_h.index
│   │   ├── MicroEuk90.index
│   │   ├── MicroEuk90.lookup
│   │   ├── MicroEuk90.source
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
│       └── readme.txt
├── Contamination
│   ├── AntiFam
│   │   ├── AntiFam.hmm.gz
│   │   ├── relnotes
│   │   └── version
│   ├── chm13v2.0
│   │   ├── chm13v2.0.1.bt2
│   │   ├── chm13v2.0.2.bt2
│   │   ├── chm13v2.0.3.bt2
│   │   ├── chm13v2.0.4.bt2
│   │   ├── chm13v2.0.rev.1.bt2
│   │   └── chm13v2.0.rev.2.bt2
│   └── kmers
│       └── ribokmers.fa.gz
└── MarkerSets
    ├── Archaea_76.hmm.gz
    ├── Bacteria_71.hmm.gz
    ├── CPR_43.hmm.gz
    ├── eukaryota_odb10.hmm.gz
    ├── eukaryota_odb10.scores_cutoff.tsv.gz
    ├── Fungi_593.hmm.gz
    ├── Protista_83.hmm.gz
    └── README

36 directories, 124 files

```

</summary>

<details>
	<summary> VEBA Database version: VDB_v5.2 (243 GB) </summary>

*  Added `MicrobeAnnotator-KEGG` [Zenodo: 10020074](https://zenodo.org/records/10020074) which includes KEGG module pathway information from [`MicrobeAnnotator`](https://doi.org/10.1186/s12859-020-03940-5).
*  Added `CAZy` protein sequences from [`dbCAN2`](https://academic.oup.com/nar/article/46/W1/W95/4996582)

```
tree -L 3 .
.
├── ACCESS_DATE
├── Annotate
│   ├── CAZy
│   │   └── CAZyDB.07262023.dmnd
│   ├── KOFAM
│   │   ├── ko_list
│   │   └── profiles
│   ├── MIBiG
│   │   └── mibig_v3.1.dmnd
│   ├── MicrobeAnnotator-KEGG
│   │   ├── KEGG_Bifurcating_Module_Information.pkl
│   │   ├── KEGG_Bifurcating_Module_Information.pkl.md5
│   │   ├── KEGG_Module_Information.txt
│   │   ├── KEGG_Module_Information.txt.md5
│   │   ├── KEGG_Regular_Module_Information.pkl
│   │   ├── KEGG_Regular_Module_Information.pkl.md5
│   │   ├── KEGG_Structural_Module_Information.pkl
│   │   └── KEGG_Structural_Module_Information.pkl.md5
│   ├── MicrobeAnnotator-KEGG.tar.gz
│   ├── NCBIfam-AMRFinder
│   │   ├── NCBIfam-AMRFinder.changelog.txt
│   │   ├── NCBIfam-AMRFinder.hmm.gz
│   │   └── NCBIfam-AMRFinder.tsv
│   ├── Pfam
│   │   ├── Pfam-A.hmm.gz
│   │   └── relnotes.txt
│   ├── UniRef
│   │   ├── uniref50.dmnd
│   │   ├── uniref50.release_note
│   │   ├── uniref90.dmnd
│   │   └── uniref90.release_note
│   └── VFDB
│       └── VFDB_setA_pro.dmnd
├── Classify
│   ├── CheckM2
│   │   └── uniref100.KO.1.dmnd
│   ├── CheckV
│   │   ├── genome_db
│   │   ├── hmm_db
│   │   └── README.txt
│   ├── geNomad
│   │   ├── genomad_db
│   │   ├── genomad_db.dbtype
│   │   ├── genomad_db_h
│   │   ├── genomad_db_h.dbtype
│   │   ├── genomad_db_h.index
│   │   ├── genomad_db.index
│   │   ├── genomad_db.lookup
│   │   ├── genomad_db_mapping
│   │   ├── genomad_db.source
│   │   ├── genomad_db_taxonomy
│   │   ├── genomad_integrase_db
│   │   ├── genomad_integrase_db.dbtype
│   │   ├── genomad_integrase_db_h
│   │   ├── genomad_integrase_db_h.dbtype
│   │   ├── genomad_integrase_db_h.index
│   │   ├── genomad_integrase_db.index
│   │   ├── genomad_integrase_db.lookup
│   │   ├── genomad_integrase_db.source
│   │   ├── genomad_marker_metadata.tsv
│   │   ├── genomad_mini_db -> genomad_db
│   │   ├── genomad_mini_db.dbtype
│   │   ├── genomad_mini_db_h -> genomad_db_h
│   │   ├── genomad_mini_db_h.dbtype -> genomad_db_h.dbtype
│   │   ├── genomad_mini_db_h.index -> genomad_db_h.index
│   │   ├── genomad_mini_db.index
│   │   ├── genomad_mini_db.lookup -> genomad_db.lookup
│   │   ├── genomad_mini_db_mapping -> genomad_db_mapping
│   │   ├── genomad_mini_db.source -> genomad_db.source
│   │   ├── genomad_mini_db_taxonomy -> genomad_db_taxonomy
│   │   ├── mini_set_ids
│   │   ├── names.dmp
│   │   ├── nodes.dmp
│   │   ├── plasmid_hallmark_annotation.txt
│   │   ├── version.txt
│   │   └── virus_hallmark_annotation.txt
│   ├── GTDB
│   │   ├── fastani
│   │   ├── markers
│   │   ├── mash
│   │   ├── masks
│   │   ├── metadata
│   │   ├── mrca_red
│   │   ├── msa
│   │   ├── pplacer
│   │   ├── radii
│   │   ├── split
│   │   ├── taxonomy
│   │   └── temp
│   ├── Microeukaryotic
│   │   ├── clean_ftp.sh
│   │   ├── humann_uniref50_annotations.tsv.gz
│   │   ├── md5_checksums
│   │   ├── microeukaryotic
│   │   ├── microeukaryotic.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10
│   │   ├── microeukaryotic.eukaryota_odb10.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10_h
│   │   ├── microeukaryotic.eukaryota_odb10_h.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10_h.index
│   │   ├── microeukaryotic.eukaryota_odb10.index
│   │   ├── microeukaryotic.eukaryota_odb10.lookup
│   │   ├── microeukaryotic.eukaryota_odb10.source
│   │   ├── microeukaryotic_h
│   │   ├── microeukaryotic_h.dbtype
│   │   ├── microeukaryotic_h.index
│   │   ├── microeukaryotic.index
│   │   ├── microeukaryotic.lookup
│   │   ├── microeukaryotic.source
│   │   ├── reference.eukaryota_odb10.list
│   │   ├── RELEASE_NOTES
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
│       └── readme.txt
├── Contamination
│   ├── AntiFam
│   │   ├── AntiFam.hmm.gz
│   │   ├── relnotes
│   │   └── version
│   ├── chm13v2.0
│   │   ├── chm13v2.0.1.bt2
│   │   ├── chm13v2.0.2.bt2
│   │   ├── chm13v2.0.3.bt2
│   │   ├── chm13v2.0.4.bt2
│   │   ├── chm13v2.0.rev.1.bt2
│   │   └── chm13v2.0.rev.2.bt2
│   └── kmers
│       └── ribokmers.fa.gz
└── MarkerSets
    ├── Archaea_76.hmm.gz
    ├── Bacteria_71.hmm.gz
    ├── CPR_43.hmm.gz
    ├── eukaryota_odb10.hmm.gz
    ├── eukaryota_odb10.scores_cutoff.tsv.gz
    ├── Fungi_593.hmm.gz
    ├── Protista_83.hmm.gz
    └── README

37 directories, 112 files
```

</details>

<details>
	<summary> VEBA Database version: VDB_v5.1 </summary>

* `VDB_v5` → `VDB_v5.1` updates `GTDB` database from `r207_v2` → `r214`.  
* Changes `${VEBA_DATABASE}/Classify/GTDBTk` → `${VEBA_DATABASE}/Classify/GTDB`.
* Adds `gtdb_r214.msh` to `${VEBA_DATABASE}/Classify/GTDB/mash/` for ANI screens. [Available on Zenodo:8048187](https://zenodo.org/record/8048187)
* Adds `VFDB` to `${VEBA_DATABASE}/Annotate/VFDB/VFDB_setA_pro.dmnd` for virulence factor annotation.

```
tree -L 3 .
├── ACCESS_DATE
├── Annotate
│   ├── KOFAM
│   │   ├── ko_list
│   │   └── profiles
│   ├── MIBiG
│   │   └── mibig_v3.1.dmnd
│   ├── NCBIfam-AMRFinder
│   │   ├── NCBIfam-AMRFinder.changelog.txt
│   │   ├── NCBIfam-AMRFinder.hmm.gz
│   │   └── NCBIfam-AMRFinder.tsv
│   ├── Pfam
│   │   ├── Pfam-A.hmm.gz
│   │   └── relnotes.txt
│   └── UniRef
│       ├── uniref50.dmnd
│       ├── uniref50.release_note
│       ├── uniref90.dmnd
│       └── uniref90.release_note
├── Classify
│   ├── CheckM2
│   │   └── uniref100.KO.1.dmnd
│   ├── CheckV
│   │   ├── genome_db
│   │   ├── hmm_db
│   │   └── README.txt
│   ├── geNomad
│   │   ├── genomad_db
│   │   ├── genomad_db.dbtype
│   │   ├── genomad_db_h
│   │   ├── genomad_db_h.dbtype
│   │   ├── genomad_db_h.index
│   │   ├── genomad_db.index
│   │   ├── genomad_db.lookup
│   │   ├── genomad_db_mapping
│   │   ├── genomad_db.source
│   │   ├── genomad_db_taxonomy
│   │   ├── genomad_integrase_db
│   │   ├── genomad_integrase_db.dbtype
│   │   ├── genomad_integrase_db_h
│   │   ├── genomad_integrase_db_h.dbtype
│   │   ├── genomad_integrase_db_h.index
│   │   ├── genomad_integrase_db.index
│   │   ├── genomad_integrase_db.lookup
│   │   ├── genomad_integrase_db.source
│   │   ├── genomad_marker_metadata.tsv
│   │   ├── genomad_mini_db -> genomad_db
│   │   ├── genomad_mini_db.dbtype
│   │   ├── genomad_mini_db_h -> genomad_db_h
│   │   ├── genomad_mini_db_h.dbtype -> genomad_db_h.dbtype
│   │   ├── genomad_mini_db_h.index -> genomad_db_h.index
│   │   ├── genomad_mini_db.index
│   │   ├── genomad_mini_db.lookup -> genomad_db.lookup
│   │   ├── genomad_mini_db_mapping -> genomad_db_mapping
│   │   ├── genomad_mini_db.source -> genomad_db.source
│   │   ├── genomad_mini_db_taxonomy -> genomad_db_taxonomy
│   │   ├── mini_set_ids
│   │   ├── names.dmp
│   │   ├── nodes.dmp
│   │   ├── plasmid_hallmark_annotation.txt
│   │   ├── version.txt
│   │   └── virus_hallmark_annotation.txt
│   ├── GTDB
│   │   ├── fastani
│   │   ├── markers
│   │   ├── mash
│   │   ├── masks
│   │   ├── metadata
│   │   ├── mrca_red
│   │   ├── msa
│   │   ├── pplacer
│   │   ├── radii
│   │   ├── split
│   │   ├── taxonomy
│   │   └── temp
│   ├── Microeukaryotic
│   │   ├── humann_uniref50_annotations.tsv.gz
│   │   ├── md5_checksums
│   │   ├── microeukaryotic
│   │   ├── microeukaryotic.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10
│   │   ├── microeukaryotic.eukaryota_odb10.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10_h
│   │   ├── microeukaryotic.eukaryota_odb10_h.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10_h.index
│   │   ├── microeukaryotic.eukaryota_odb10.index
│   │   ├── microeukaryotic.eukaryota_odb10.lookup
│   │   ├── microeukaryotic.eukaryota_odb10.source
│   │   ├── microeukaryotic_h
│   │   ├── microeukaryotic_h.dbtype
│   │   ├── microeukaryotic_h.index
│   │   ├── microeukaryotic.index
│   │   ├── microeukaryotic.lookup
│   │   ├── microeukaryotic.source
│   │   ├── reference.eukaryota_odb10.list
│   │   ├── RELEASE_NOTES
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
│       └── readme.txt
├── Contamination
│   ├── AntiFam
│   │   ├── AntiFam.hmm.gz
│   │   ├── relnotes
│   │   └── version
│   ├── chm13v2.0
│   │   ├── chm13v2.0.1.bt2
│   │   ├── chm13v2.0.2.bt2
│   │   ├── chm13v2.0.3.bt2
│   │   ├── chm13v2.0.4.bt2
│   │   ├── chm13v2.0.rev.1.bt2
│   │   └── chm13v2.0.rev.2.bt2
│   └── kmers
│       └── ribokmers.fa.gz
└── MarkerSets
    ├── Archaea_76.hmm.gz
    ├── Bacteria_71.hmm.gz
    ├── CPR_43.hmm.gz
    ├── eukaryota_odb10.hmm.gz
    ├── eukaryota_odb10.scores_cutoff.tsv.gz
    ├── Fungi_593.hmm.gz
    ├── Protista_83.hmm.gz
    └── README
```

</details>

<details>
	<summary> VEBA Database version: VDB_v5 </summary>
	
`VDB_v4` → `VDB_v5` replaces `nr` with `UniRef90` and `UniRef50`.  Also includes `MiBIG` database.

```
tree -L 3 .
├── ACCESS_DATE
├── Annotate
│   ├── KOFAM
│   │   ├── ko_list
│   │   └── profiles
│   ├── MIBiG
│   │   └── mibig_v3.1.dmnd
│   ├── NCBIfam-AMRFinder
│   │   ├── NCBIfam-AMRFinder.changelog.txt
│   │   ├── NCBIfam-AMRFinder.hmm.gz
│   │   └── NCBIfam-AMRFinder.tsv
│   ├── Pfam
│   │   ├── Pfam-A.hmm.gz
│   │   └── relnotes.txt
│   └── UniRef
│       ├── uniref50.dmnd
│       └── uniref90.dmnd
├── Classify
│   ├── CheckM2
│   │   └── uniref100.KO.1.dmnd
│   ├── CheckV
│   │   ├── genome_db
│   │   ├── hmm_db
│   │   └── README.txt
│   ├── geNomad
│   │   ├── genomad_db
│   │   ├── genomad_db.dbtype
│   │   ├── genomad_db_h
│   │   ├── genomad_db_h.dbtype
│   │   ├── genomad_db_h.index
│   │   ├── genomad_db.index
│   │   ├── genomad_db.lookup
│   │   ├── genomad_db_mapping
│   │   ├── genomad_db.source
│   │   ├── genomad_db_taxonomy
│   │   ├── genomad_integrase_db
│   │   ├── genomad_integrase_db.dbtype
│   │   ├── genomad_integrase_db_h
│   │   ├── genomad_integrase_db_h.dbtype
│   │   ├── genomad_integrase_db_h.index
│   │   ├── genomad_integrase_db.index
│   │   ├── genomad_integrase_db.lookup
│   │   ├── genomad_integrase_db.source
│   │   ├── genomad_marker_metadata.tsv
│   │   ├── genomad_mini_db -> genomad_db
│   │   ├── genomad_mini_db.dbtype
│   │   ├── genomad_mini_db_h -> genomad_db_h
│   │   ├── genomad_mini_db_h.dbtype -> genomad_db_h.dbtype
│   │   ├── genomad_mini_db_h.index -> genomad_db_h.index
│   │   ├── genomad_mini_db.index
│   │   ├── genomad_mini_db.lookup -> genomad_db.lookup
│   │   ├── genomad_mini_db_mapping -> genomad_db_mapping
│   │   ├── genomad_mini_db.source -> genomad_db.source
│   │   ├── genomad_mini_db_taxonomy -> genomad_db_taxonomy
│   │   ├── mini_set_ids
│   │   ├── names.dmp
│   │   ├── nodes.dmp
│   │   ├── plasmid_hallmark_annotation.txt
│   │   ├── version.txt
│   │   └── virus_hallmark_annotation.txt
│   ├── GTDBTk
│   │   ├── fastani
│   │   ├── markers
│   │   ├── masks
│   │   ├── metadata
│   │   ├── mrca_red
│   │   ├── msa
│   │   ├── pplacer
│   │   ├── radii
│   │   ├── split
│   │   ├── taxonomy
│   │   └── temp
│   ├── Microeukaryotic
│   │   ├── humann_uniref50_annotations.tsv.gz
│   │   ├── md5_checksums
│   │   ├── microeukaryotic
│   │   ├── microeukaryotic.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10
│   │   ├── microeukaryotic.eukaryota_odb10.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10_h
│   │   ├── microeukaryotic.eukaryota_odb10_h.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10_h.index
│   │   ├── microeukaryotic.eukaryota_odb10.index
│   │   ├── microeukaryotic.eukaryota_odb10.lookup
│   │   ├── microeukaryotic.eukaryota_odb10.source
│   │   ├── microeukaryotic_h
│   │   ├── microeukaryotic_h.dbtype
│   │   ├── microeukaryotic_h.index
│   │   ├── microeukaryotic.index
│   │   ├── microeukaryotic.lookup
│   │   ├── microeukaryotic.source
│   │   ├── reference.eukaryota_odb10.list
│   │   ├── RELEASE_NOTES
│   │   ├── source_taxonomy.tsv.gz
│   │   ├── source_to_lineage.dict.pkl.gz
│   │   └── target_to_source.dict.pkl.gz
│   └── NCBITaxonomy
│       ├── citations.dmp
│       ├── delnodes.dmp
│       ├── division.dmp
│       ├── gc.prt
│       ├── gencode.dmp
│       ├── images.dmp
│       ├── merged.dmp
│       ├── names.dmp
│       ├── nodes.dmp
│       ├── prot.accession2taxid.FULL.gz
│       └── readme.txt
├── Contamination
│   ├── AntiFam
│   │   ├── AntiFam.hmm.gz
│   │   ├── relnotes
│   │   └── version
│   ├── chm13v2.0
│   │   ├── chm13v2.0.1.bt2
│   │   ├── chm13v2.0.2.bt2
│   │   ├── chm13v2.0.3.bt2
│   │   ├── chm13v2.0.4.bt2
│   │   ├── chm13v2.0.rev.1.bt2
│   │   └── chm13v2.0.rev.2.bt2
│   └── kmers
│       └── ribokmers.fa.gz
└── MarkerSets
    ├── Archaea_76.hmm.gz
    ├── Bacteria_71.hmm.gz
    ├── CPR_43.hmm.gz
    ├── eukaryota_odb10.hmm.gz
    ├── eukaryota_odb10.scores_cutoff.tsv.gz
    ├── Fungi_593.hmm.gz
    ├── Protista_83.hmm.gz
    └── README
```

</details>

<details>
	<summary> VEBA Database version: VDB_v4 </summary>
	

`VDB_v4` is `VDB_v3.1` with the following changes: 1) `CheckM1` database swapped for `CheckM2` database;  includes `geNomad` database; and 3) updates `CheckV` database.  Refer to [development log](https://github.com/jolespin/veba/blob/main/DEVELOPMENT.md#release-v11-currently-testing-before-official-release) for specifics.

```
tree -L 3 .
.
├── ACCESS_DATE
├── Annotate
│   ├── KOFAM
│   │   ├── ko_list
│   │   └── profiles
│   ├── NCBIfam-AMRFinder
│   │   ├── NCBIfam-AMRFinder.changelog.txt
│   │   ├── NCBIfam-AMRFinder.hmm.gz
│   │   └── NCBIfam-AMRFinder.tsv
│   ├── nr
│   │   └── nr.dmnd
│   └── Pfam
│       ├── Pfam-A.hmm.gz
│       └── relnotes.txt
├── Classify
│   ├── CheckM2
│   │   └── uniref100.KO.1.dmnd
│   ├── CheckV
│   │   ├── genome_db
│   │   ├── hmm_db
│   │   └── README.txt
│   ├── geNomad
│   │   ├── genomad_db
│   │   ├── genomad_db.dbtype
│   │   ├── genomad_db_h
│   │   ├── genomad_db_h.dbtype
│   │   ├── genomad_db_h.index
│   │   ├── genomad_db.index
│   │   ├── genomad_db.lookup
│   │   ├── genomad_db_mapping
│   │   ├── genomad_db.source
│   │   ├── genomad_db_taxonomy
│   │   ├── genomad_integrase_db
│   │   ├── genomad_integrase_db.dbtype
│   │   ├── genomad_integrase_db_h
│   │   ├── genomad_integrase_db_h.dbtype
│   │   ├── genomad_integrase_db_h.index
│   │   ├── genomad_integrase_db.index
│   │   ├── genomad_integrase_db.lookup
│   │   ├── genomad_integrase_db.source
│   │   ├── genomad_marker_metadata.tsv
│   │   ├── genomad_mini_db -> genomad_db
│   │   ├── genomad_mini_db.dbtype
│   │   ├── genomad_mini_db_h -> genomad_db_h
│   │   ├── genomad_mini_db_h.dbtype -> genomad_db_h.dbtype
│   │   ├── genomad_mini_db_h.index -> genomad_db_h.index
│   │   ├── genomad_mini_db.index
│   │   ├── genomad_mini_db.lookup -> genomad_db.lookup
│   │   ├── genomad_mini_db_mapping -> genomad_db_mapping
│   │   ├── genomad_mini_db.source -> genomad_db.source
│   │   ├── genomad_mini_db_taxonomy -> genomad_db_taxonomy
│   │   ├── mini_set_ids
│   │   ├── names.dmp
│   │   ├── nodes.dmp
│   │   ├── plasmid_hallmark_annotation.txt
│   │   ├── version.txt
│   │   └── virus_hallmark_annotation.txt
│   ├── GTDBTk
│   │   ├── fastani
│   │   ├── markers
│   │   ├── masks
│   │   ├── metadata
│   │   ├── mrca_red
│   │   ├── msa
│   │   ├── pplacer
│   │   ├── radii
│   │   ├── split
│   │   ├── taxonomy
│   │   └── temp
│   ├── Microeukaryotic
│   │   ├── humann_uniref50_annotations.tsv.gz
│   │   ├── md5_checksums
│   │   ├── microeukaryotic
│   │   ├── microeukaryotic.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10
│   │   ├── microeukaryotic.eukaryota_odb10.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10_h
│   │   ├── microeukaryotic.eukaryota_odb10_h.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10_h.index
│   │   ├── microeukaryotic.eukaryota_odb10.index
│   │   ├── microeukaryotic.eukaryota_odb10.lookup
│   │   ├── microeukaryotic.eukaryota_odb10.source
│   │   ├── microeukaryotic_h
│   │   ├── microeukaryotic_h.dbtype
│   │   ├── microeukaryotic_h.index
│   │   ├── microeukaryotic.index
│   │   ├── microeukaryotic.lookup
│   │   ├── microeukaryotic.source
│   │   ├── reference.eukaryota_odb10.list
│   │   ├── RELEASE_NOTES
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
│       └── readme.txt
├── Contamination
│   ├── AntiFam
│   │   ├── AntiFam.hmm.gz
│   │   ├── relnotes
│   │   └── version
│   ├── chm13v2.0
│   │   ├── chm13v2.0.1.bt2
│   │   ├── chm13v2.0.2.bt2
│   │   ├── chm13v2.0.3.bt2
│   │   ├── chm13v2.0.4.bt2
│   │   ├── chm13v2.0.rev.1.bt2
│   │   └── chm13v2.0.rev.2.bt2
│   └── kmers
│       └── ribokmers.fa.gz
└── MarkerSets
    ├── Archaea_76.hmm.gz
    ├── Bacteria_71.hmm.gz
    ├── CPR_43.hmm.gz
    ├── eukaryota_odb10.hmm.gz
    ├── eukaryota_odb10.scores_cutoff.tsv.gz
    ├── Fungi_593.hmm.gz
    ├── Protista_83.hmm.gz
    └── README

31 directories, 96 files
```

</details>


<details>
	<summary>VEBA Database version: VDB_v3.1</summary>
	
The same as `VDB_v3` but updates `VDB-Microeukaryotic_v2` to `VDB-Microeukaryotic_v2.1` which has a `reference.eukaryota_odb10.list` containing only the subset of identifiers that core eukaryotic markers (useful for classification).
	
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
│   │   ├── markers
│   │   ├── masks
│   │   ├── metadata
│   │   ├── mrca_red
│   │   ├── msa
│   │   ├── pplacer
│   │   ├── radii
│   │   ├── split
│   │   ├── taxonomy
│   │   └── temp
│   ├── Microeukaryotic
│   │   ├── humann_uniref50_annotations.tsv.gz
│   │   ├── md5_checksums
│   │   ├── microeukaryotic
│   │   ├── microeukaryotic.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10
│   │   ├── microeukaryotic.eukaryota_odb10.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10_h
│   │   ├── microeukaryotic.eukaryota_odb10_h.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10_h.index
│   │   ├── microeukaryotic.eukaryota_odb10.index
│   │   ├── microeukaryotic.eukaryota_odb10.lookup
│   │   ├── microeukaryotic.eukaryota_odb10.source
│   │   ├── microeukaryotic_h
│   │   ├── microeukaryotic_h.dbtype
│   │   ├── microeukaryotic_h.index
│   │   ├── microeukaryotic.index
│   │   ├── microeukaryotic.lookup
│   │   ├── microeukaryotic.source
│   │   ├── reference.eukaryota_odb10.list
│   │   ├── reference.eukaryota_odb10.list.md5
│   │   ├── reference.faa.gz
│   │   ├── RELEASE_NOTES
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

35 directories, 60 files
```

</details>


<details>
	<summary>VEBA Database version: VDB_v3</summary>

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
│   │   ├── markers
│   │   ├── masks
│   │   ├── metadata
│   │   ├── mrca_red
│   │   ├── msa
│   │   ├── pplacer
│   │   ├── radii
│   │   ├── split
│   │   ├── taxonomy
│   │   └── temp
│   ├── Microeukaryotic
│   │   ├── humann_uniref50_annotations.tsv.gz
│   │   ├── md5_checksums
│   │   ├── microeukaryotic
│   │   ├── microeukaryotic.dbtype
│   │   ├── microeukaryotic_h
│   │   ├── microeukaryotic_h.dbtype
│   │   ├── microeukaryotic_h.index
│   │   ├── microeukaryotic.index
│   │   ├── microeukaryotic.lookup
│   │   ├── microeukaryotic.source
│   │   ├── reference.faa.gz
│   │   ├── RELEASE_NOTES
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

35 directories, 50 files
```

</details>


<details>
	<summary>VEBA Database version: VDB_v2</summary>
	
* Compatible with *VEBA* version: `v1.0.2a+`
	

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
│   │   ├── markers
│   │   ├── masks
│   │   ├── metadata
│   │   ├── mrca_red
│   │   ├── msa
│   │   ├── pplacer
│   │   ├── radii
│   │   ├── split
│   │   ├── taxonomy
│   │   └── temp
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

35 directories, 47 files

```
</details>


<details>
	<summary>VEBA Database version: VDB_v1</summary>
	

* Compatible with *VEBA* version: `v1.0.0`, `v1.0.1`
	
	
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
│   │   ├── markers
│   │   ├── masks
│   │   ├── metadata
│   │   ├── mrca_red
│   │   ├── msa
│   │   ├── pplacer
│   │   ├── radii
│   │   ├── split
│   │   ├── taxonomy
│   │   └── temp
│   ├── Microeukaryotic
│   │   ├── humann_uniref50_annotations.tsv.gz
│   │   ├── md5_checksums
│   │   ├── microeukaryotic
│   │   ├── microeukaryotic.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10
│   │   ├── microeukaryotic.eukaryota_odb10.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10_h
│   │   ├── microeukaryotic.eukaryota_odb10_h.dbtype
│   │   ├── microeukaryotic.eukaryota_odb10_h.index
│   │   ├── microeukaryotic.eukaryota_odb10.index
│   │   ├── microeukaryotic.eukaryota_odb10.lookup
│   │   ├── microeukaryotic.eukaryota_odb10.source
│   │   ├── microeukaryotic_h
│   │   ├── microeukaryotic_h.dbtype
│   │   ├── microeukaryotic_h.index
│   │   ├── microeukaryotic.index
│   │   ├── microeukaryotic.lookup
│   │   ├── microeukaryotic.source
│   │   ├── reference.eukaryota_odb10.list
│   │   ├── reference.eukaryota_odb10.list.md5
│   │   ├── reference.faa.gz
│   │   ├── RELEASE_NOTES
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

35 directories, 60 files
```
	
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

___________________________________________________________

#### Version Notes:

* For the human contamination, if you use `KneadData` and already have a `Bowtie2` index for human then you can use that instead.  The only module that uses this is `preprocess.py` and you have to specify this directly when running (i.e., it's optional) so it doesn't matter if it's in the database directory or not (same with ribokmers.fa.gz). 
* `CheckV`, `GTDB-Tk`, and `geNomad` require specific versions so don't mess with those unless you know what you're doing.
* If you are low on disk space and already have these installed then just symlink them with the structure above. If so, them just comment out those sections of `download_databases.sh`.  

_______________________________________________________

If you have any issues, please create a GitHub issue with the prefix `[Database]` followed by the question.