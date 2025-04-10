### VEBA Database:

# VEBA Database Version Changes

This document summarizes the main changes introduced between different versions of the VEBA Database (VDB).

## Summary of Changes

### VDB_v8.1 (from VDB_v8)
*   Added `kofam.enzymes.list` and `kofam.pathways.list` in `Annotate/KOfam/` to provide subsets for `pykofamsearch`.
*   Updated `Annotate/KOfam/` with serialized KOfam data including enzyme support (replacing older `ko_list` and `profiles` structure).
*   Minor size increase (300 GB -> 306 GB).

### VDB_v8 (from VDB_v7)
*   Replaced the `Annotate/MicrobeAnnotator-KEGG` database implementation with the `Annotate/KEGG-Pathway-Profiler` database.
*   Slight size decrease (301 GB -> 300 GB).

### VDB_v7 (from VDB_v6)
*   Updated GTDB database from r214.1 to r220 (Note: Specific GTDB version files not explicitly listed in the provided tree, change based on summary text).
*   Added Pfam clan information (`Annotate/Pfam/Pfam-A.clans.tsv.gz`).
*   Renamed database directory `Annotate/KOFAM` to `Annotate/KOfam`.
*   Size increase (272 GB -> 301 GB).

### VDB_v6 (from VDB_v5.2)
*   Added/Updated `MicroEuk` database (likely `MicroEuk_v3` based on summary text). The directory structure under `Classify/MicroEuk/` appears more comprehensive than `Classify/Microeukaryotic/` in v5.2.
*   Size increase (243 GB -> 272 GB).

### VDB_v5.2 (from VDB_v5.1)
*   Added `Annotate/MicrobeAnnotator-KEGG` database based on [Zenodo: 10020074](https://zenodo.org/records/10020074).
*   Added `Annotate/CAZy` database (DIAMOND format) from dbCAN2.
*   Size increased (No size for v5.1 provided -> 243 GB).

### VDB_v5.1 (from VDB_v5)
*   Updated GTDB database from r207_v2 to r214.
*   Renamed database directory `Classify/GTDBTk` to `Classify/GTDB`.
*   Added `gtdb_r214.msh` to `Classify/GTDB/mash/` for ANI screening.
*   Added `Annotate/VFDB/VFDB_setA_pro.dmnd` for virulence factor annotation.

### VDB_v5 (from VDB_v4)
*   Replaced `Annotate/nr` (NCBI nr DIAMOND database) with `Annotate/UniRef/uniref90.dmnd` and `Annotate/UniRef/uniref50.dmnd`.
*   Added `Annotate/MIBiG/mibig_v3.1.dmnd` database.

### VDB_v4 (from VDB_v3.1)
*   Replaced `Classify/CheckM` (CheckM1 database) with `Classify/CheckM2` database.
*   Added `Classify/geNomad` database.
*   Updated `Classify/CheckV` database (specific file changes not detailed in summary).

### VDB_v3.1 (from VDB_v3)
*   Updated `Classify/Microeukaryotic` database (`VDB-Microeukaryotic_v2` to `VDB-Microeukaryotic_v2.1`).
*   Specifically added `Classify/Microeukaryotic/reference.eukaryota_odb10.list`.

### VDB_v3 (from VDB_v2)
*   Added `eukaryota_odb10` related files to `Classify/Microeukaryotic/`.
*   Added `eukaryota_odb10` marker HMMs and score cutoffs to `MarkerSets/`.
*   Added `taxa.sqlite` and `taxa.sqlite.traverse.pkl` to `Classify/NCBITaxonomy/` (likely related to Microeukaryotic updates).

### VDB_v2 (from VDB_v1)
*   Updated `Classify/Microeukaryotic` database structure and reference files.
*   Replaced `Contamination/grch38` human genome index with `Contamination/chm13v2.0` index.

## Compatibility Table

| VEBA Database Version | Compatible VEBA Software Version(s)       | Legacy Database Name |
|-----------------------|-------------------------------------------|----------------------|
| VEBA-DB_v9            | v2.5.0                                    | -                    |
| VEBA-DB_v8.1          | v2.4.0, v2.4.1, v2.4.2                    | VDB_v8.1             |
| VEBA-DB_v8            | v2.3.0                                    | VDB_v8               |
| VEBA-DB_v7            | v2.2.0, v2.2.1                            | VDB_v7               |
| VEBA-DB_v6            | v1.4.1, v1.4.2, v1.5.0, v2.0.0, v2.1.0    | VDB_v6               |
| VEBA-DB_v5.2          | v1.3.0                                    | VDB_v5.2             |
| VEBA-DB_v5.1          | v1.2.0                                    | VDB_v5.1             |
| VEBA-DB_v5            | v1.1.2                                    | VDB_v5               |
| VEBA-DB_v4            | v1.1.0                                    | VDB_v4               |
| VEBA-DB_v3.1          | *Likely v1.1.0*¹                          | VDB_v3.1             |
| VEBA-DB_v3            | v1.0.4                                    | VDB_v3               |
| VEBA-DB_v2            | v1.0.2a+, v1.0.3e                         | VDB_v2               |
| VEBA-DB_v1            | v1.0.0, v1.0.1                            | VDB_v1               |

*¹ Note: Version `1.1.0` explicitly mentions updating *from* `VDB_v3.1` to `VDB_v4`, implying compatibility with `VDB_v3.1` before the update.*

*Note: Please update the "Information needed" fields as compatibility information becomes available.*

## Profile HMM Sources:
Please cite the following sources if these marker sets are used in any way:

```
* Archaea_76.hmm.gz - (Anvi'o) Lee, https://doi.org/10.1093/bioinformatics/btz188 (https://github.com/merenlab/anvio/tree/master/anvio/data/hmm/Archaea_76)

* Bacteria_71.hmm.gz - (Anvi'o) Lee modified, https://doi.org/10.1093/bioinformatics/btz188 (https://github.com/merenlab/anvio/tree/master/anvio/data/hmm/Bacteria_71)

* Protista_83.hmm.gz - (Anvi'o) Delmont, http://merenlab.org/delmont-euk-scgs (https://github.com/merenlab/anvio/tree/master/anvio/data/hmm/Protista_83)

* Fungi_593.hmm.gz - (FGMP) https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2782-9

* CPR_43.hmm.gz - (CheckM) https://github.com/Ecogenomics/CheckM/tree/master/custom_marker_sets

* eukaryota_odb12.hmm.gz - (BUSCO) https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb12.2025-02-25.tar.gz
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


##

___________________________________________________________

#### Version Notes:

* For the human contamination, if you use `KneadData` and already have a `Bowtie2` index for human then you can use that instead.  The only module that uses this is `preprocess.py` and you have to specify this directly when running (i.e., it's optional) so it doesn't matter if it's in the database directory or not (same with ribokmers.fa.gz). 
* `CheckV`, `GTDB-Tk`, and `geNomad` require specific versions so don't mess with those unless you know what you're doing.
* If you are low on disk space and already have these installed then just symlink them with the structure above. If so, them just comment out those sections of `download_databases.sh`.  

_______________________________________________________

If you have any issues, please create a GitHub issue with the prefix `[Database]` followed by the question.