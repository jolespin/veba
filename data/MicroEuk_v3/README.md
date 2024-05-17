# VEBA Microeukaryotic Protein Database (MicroEuk100/90/50)

Since we are all about [FAIR principles for scientific stewardship](https://www.go-fair.org/fair-principles/), here is a step-by-step guide to recreating VEBA's Microeukaryotic Protein Database (now referred to as MicroEuk100/90/50).

| Dataset     | Number of sequences | File size (Gzipped) | Compression of MicroEuk100 |
|-------------|---------------------|---------------------|----------------------------|
| MicroEuk100 | 79,920,431            | 19 GB               | 0%                         |
| MicroEuk90  | 51,767,730            | 13 GB               | 35.20%                     |
| MicroEuk50  | 29,898,853            | 6.5 GB              | 62.60%                     |

The following datasetes are used in the following priority:

1. [JGI MycoCosm](https://mycocosm.jgi.doe.gov/mycocosm/home)
2. [JGI PhycoCosm](https://phycocosm.jgi.doe.gov/phycocosm/home)
3. [EnsemblProtists](https://protists.ensembl.org/index.html)
4. [MMETSP](https://zenodo.org/records/257410)
5. [TARA_SAGv1](https://www.genoscope.cns.fr/tara/)
6. [EukProt](https://figshare.com/articles/dataset/EukProt_a_database_of_genome-scale_predicted_proteins_across_the_diversity_of_eukaryotic_life/12417881)
7. [EukZoo](https://zenodo.org/records/1476236)
8. [TARA_SMAGv1](https://www.genoscope.cns.fr/tara/)
9. [NR_Protists-Fungi](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/)

Duplicates in lower priority datasets are discarded.

For each dataset, build a source table that has the following columns: `source_taxonomy.tsv.gz`

```
['dataset', 'ncbi_taxid', 'lineage', 'class', 'order', 'family', 'genus', 'species', 'strain', 'notes', 'resolved_lineage']
```

Dataset is mandatory and having class, order, family, genus, and species is preferred.  Strain is not used.  Consensus eukaryotic classification can handle missing fields but, obviously, the more complete the better the classification.

**Note:** `TARA_SMAGv1` does not have species-level assignments so the MAG identifier was used for species.

**Example:** 

```
id_source	dataset	ncbi_taxid	lineage	class	order	family	genus	species	strain	notes	resolved_lineage
Aalte1	MycoCosm	5599	d__Eukaryota;p__Ascomycota;c__Dothideomycetes;o__Pleosporales;f__Pleosporaceae;g__Alternaria;s__Alternaria alternata	Dothideomycetes	Pleosporales	Pleosporaceae	Alternaria	Alternaria alternata			True
Aaoar1	MycoCosm	1450171	d__Eukaryota;p__Ascomycota;c__Dothideomycetes;o__Pleosporales;f__;g__Aaosphaeria;s__Aaosphaeria arxii	DothideomycetesPleosporales		Aaosphaeria	Aaosphaeria arxii			False
Abobi1	MycoCosm	137743	d__Eukaryota;p__Basidiomycota;c__Agaricomycetes;o__Polyporales;f__Podoscyphaceae;g__Abortiporus;s__Abortiporus biennis	Agaricomycetes	Polyporales	Podoscyphaceae	Abortiporus	Abortiporus biennis			True
Abobie1	MycoCosm	137743	d__Eukaryota;p__Basidiomycota;c__Agaricomycetes;o__Polyporales;f__Podoscyphaceae;g__Abortiporus;s__Abortiporus biennis	Agaricomycetes	Polyporales	Podoscyphaceae	Abortiporus	Abortiporus biennis			True
Abscae1	MycoCosm	90261	d__Eukaryota;p__Mucoromycota;c__Mucoromycetes;o__Mucorales;f__Cunninghamellaceae;g__Absidia;s__Absidia caerulea	Mucoromycetes	Mucorales	Cunninghamellaceae	Absidia	Absidia caerulea			True
Absrep1	MycoCosm	90262	d__Eukaryota;p__Mucoromycota;c__Mucoromycetes;o__Mucorales;f__Cunninghamellaceae;g__Absidia;s__Absidia repens	Mucoromycetes	Mucorales	Cunninghamellaceae	Absidia	Absidia repens			True
Acain1	MycoCosm	215250	d__Eukaryota;p__Basidiomycota;c__Exobasidiomycetes;o__Exobasidiales;f__Cryptobasidiaceae;g__Acaromyces;s__Acaromyces ingoldii	Exobasidiomycetes	Exobasidiales	Cryptobasidiaceae	Acaromyces	Acaromyces ingoldii			True
Acastr1	MycoCosm	1307806	d__Eukaryota;p__Ascomycota;c__Lecanoromycetes;o__Acarosporales;f__Acarosporaceae;g__Acarospora;s__Acarospora strigata	Lecanoromycetes	Acarosporales	Acarosporaceae	Acarospora	Acarospora strigata			True
Acema1	MycoCosm	886606	d__Eukaryota;p__Ascomycota;c__Leotiomycetes;o__Helotiales;f__Mollisiaceae;g__Acephala;s__Acephala macrosclerotiorum	Leotiomycetes	Helotiales	Mollisiaceae	Acephala	Acephala macrosclerotiorum			True
```


#### Calculating md5 hash for sequences (without stop codon)

There may an instance where a sequence has the same name in another database.  There is also the situation of more than one identifier referring to the same actual sequence.  To get around this, we want to use md5 hash identifiers that are unique to a sequence as the sequence identifier.  

After downloading all of the fasta files from the different protein databases, I gzipped them and made the following file structure:

```
[dataset]
[dataset]/proteins/[id_source].fasta.gz
[dataset]/id_to_hash.no_stop_codon.tsv
```

**Note:** `EukProt` has a missing linebreak in one or more of files.  To fix this, add an extra linebreak at the *beginning* of `EP00737_Agarophyton_chilense.fasta`.

```
cat MycoCosm/proteins/*/*.fasta.gz | seqkit replace --by-seq -p '\*$' | seqkit fx2tab -i -s -n > MycoCosm/id_to_hash.no_stop_codon.tsv

cat PhycoCosm/proteins/*/*.fasta.gz | seqkit replace --by-seq -p '\*$' | seqkit fx2tab -i -s -n > PhycoCosm/id_to_hash.no_stop_codon.tsv

cat EnsemblProtists/proteins/*/*.fa.gz | seqkit replace --by-seq -p '\*$' | seqkit fx2tab -i -s -n > EnsemblProtists/id_to_hash.no_stop_codon.tsv

cat MMETSP/proteins/*.pep.gz |  seqkit replace --by-seq -p '\*$' | seqkit fx2tab -i -s -n > MMETSP/id_to_hash.no_stop_codon.tsv

cat TARA_SAGv1/proteins/*.faa.gz | seqkit replace --by-seq -p '\*$' | seqkit fx2tab -i -s -n > TARA_SAGv1/id_to_hash.no_stop_codon.tsv

cat EukProt/proteins/*.fasta.gz | seqkit replace --by-seq -p '\*$' | seqkit fx2tab -i -s -n > EukProt/id_to_hash.no_stop_codon.tsv

cat EukZoo/proteins/*.fasta.gz | seqkit replace --by-seq -p '\*$' | seqkit fx2tab -i -s -n > EukZoo/id_to_hash.no_stop_codon.tsv

cat TARA_SMAGv1/proteins/*.fa.gz | seqkit replace --by-seq -p '\*$' | seqkit fx2tab -i -s -n > TARA_SMAGv1/id_to_hash.no_stop_codon.tsv

cat NR_Protists-Fungi/reference.no_metazoans.faa.gz |  seqkit replace --by-seq -p '\*$' | seqkit fx2tab -i -s -n > NR_Protists-Fungi/id_to_hash.no_stop_codon.tsv

```

#### Build an identifier mapping table for the sequences

Next we want to build an identifier mapping table that we can use for accounting of the protein sequences incorporated and where they came from with respect to the original datasets.

The `identifier_mapping.proteins.tsv.gz` must have the following fields fields (no header):

1. Dataset name
2. Source identifier for the organism
3. Original protein identifier
4. MD5 hash identifier for the protein without a stop codon


**Example:** 

```
MycoCosm	Lenfl1	jgi|Lenfl1|6812|CE6811_1407	370400ec9621f03f7a13a7f6f9741570
MycoCosm	Lenfl1	jgi|Lenfl1|473551|estExt_fgenesh1_pm.C_1_t10324	6bf5c7553a03c589a67aa512813fadb8
MycoCosm	Lenfl1	jgi|Lenfl1|285089|e_gw1.1.922.1	8be5aa97999c57aa9bc0d9d91be084d4
MycoCosm	Lenfl1	jgi|Lenfl1|381725|gm1.497_g	f0ffb14b7542da089c152dfc2b97b2e2
MycoCosm	Lenfl1	jgi|Lenfl1|353978|estExt_Genewise1Plus.C_1_t30174	5a26c5afd3d9045ed065a85c3e17198a
MycoCosm	Lenfl1	jgi|Lenfl1|381723|gm1.495_g	2888c78fce32572d2d8bc3da9980cc34
MycoCosm	Lenfl1	jgi|Lenfl1|6708|CE6707_463	37f2c8cf9696991b033312e334f8e2a0
MycoCosm	Lenfl1	jgi|Lenfl1|6709|CE6708_2421	fa03e514310bb31dff3733dfe22cf2c5
MycoCosm	Lenfl1	jgi|Lenfl1|398267|fgenesh1_pg.1_#_389	3877862f548aa27be71f25a89fd2cc5e
MycoCosm	Lenfl1	jgi|Lenfl1|381718|gm1.490_g	7af8e449bd547cc97fbc3a1fe1c209c8
```

#### Concatenate sequences to build the larger database


You can use this if fasta is correctly formatted (i.e., no linebreak issues.  See this [issue](https://github.com/shenwei356/seqkit/issues/422#issuecomment-1808731378).)

```
cat MycoCosm/proteins/*/*.fasta.gz PhycoCosm/proteins/*/*.fasta.gz EnsemblProtists/proteins/*/*.fa.gz MMETSP/proteins/*.pep.gz TARA_SAGv1/proteins/*.faa.gz EukProt/proteins/*.fasta.gz EukZoo/proteins/*.fasta.gz TARA_SMAGv1/proteins/*.pep.fa.gz NR_Protists-Fungi/reference.no_metazoans.faa.gz | seqkit seq -m 11 | seqkit replace --by-seq -p '\*$' -w 0 >  merged/references.gte11.no_hash.faa
```

Alternatively, you can use `clean_fasta.py` from VEBA to merge sequences, ignore duplicates (identifiers only), and make sure there's no linebreak issues (i.e., > characters in the sequence):

```
# clean_fasta.py is installed with VEBA
cat MycoCosm/proteins/*/*.fasta.gz PhycoCosm/proteins/*/*.fasta.gz EnsemblProtists/proteins/*/*.fa.gz MMETSP/proteins/*.pep.gz TARA_SAGv1/proteins/*.faa.gz EukProt/proteins/*.fasta.gz EukZoo/proteins/*.fasta.gz TARA_SMAGv1/proteins/*.pep.fa.gz NR_Protists-Fungi/reference.no_metazoans.faa.gz | gzip -d | clean_fasta.py -m 11 >  merged/references.gte11.no_hash.faa
```

#### Relabel proteins using the md5 hash (and log which proteins are excluded)

Now we can use Python 3.7+ for greater control of renaming the proteins using the md5 hash identifier while excluding duplicate sequences. 

##### Inputs:
* `identifier_mapping.proteins.tsv.gz` - Contains the dataset, source organism, protein id, and md5 hash id
* `references.gte11.no_hash.faa` - Concatenated protein sequences with protein ids in `identifier_mapping.proteins.tsv.gz`

##### Outputs:
* `MicroEuk100.faa` - Complete sequence set used for `MicroEuk100` with md5 hash as identifiers and each sequence is on a single line.
* `references.gte11.no_hash.discarded.list` - Protein ids that do not have a md5 hash
* `references.gte11.hash.duplicates.tsv` - Protein ids with a md5 hash but already exist in the database from a higher priority dataset.

```python
import gzip
from tqdm import tqdm
from Bio.SeqIO.FastaIO import SimpleFastaParser

id_to_hash = dict()
with gzip.open("identifier_mapping.proteins.tsv.gz", "rt") as f:
    for line in tqdm(f, total=79924653):
        line = line.strip()
        if line:
            dataset, id_source, id, id_hash = map(lambda x: x.strip(), line.split("\t"))
            assert id not in id_to_hash
            id_to_hash[id] = id_hash

f_out = open("MicroEuk100.faa", "w")
f_discarded = open("references.gte11.no_hash.discarded.list", "w")
f_duplicates = open("references.gte11.hash.duplicates.tsv", "w")

hash_ids = set()
with open("references.gte11.no_hash.faa", "r") as f_in:
    for id, seq in tqdm(SimpleFastaParser(f_in), total=127065066):
        id = id.split(" ")[0].strip()
        if id in id_to_hash:
            id_hash = id_to_hash[id]
            if id_hash not in hash_ids:
                hash_ids.add(id_hash)
                print(">{}\n{}".format(id_hash, seq), file=f_out)
            else:
                print(id, id_hash, sep="\t", file=f_duplicates)
        else:
            print(id, sep="\t", file=f_discarded)

f_out.close()
f_discarded.close()
f_duplicates.close()
```

#### Clustering MicroEuk100 for more economical databases

Performing the search space using ~80M proteins can require a fair amount of resources so we also provide clustered versions of the databases.  Here are the commands for running `MMSEQS2` for clustering using 80% query coverage. 


##### MicroEuk90

```
# mmseqs2_wrapper.py installed with VEBA
conda activate VEBA-cluster_env

FASTA=MicroEuk100.faa
mmseqs2_wrapper.py -i ${FASTA} -o mmseqs2_output/MicroEuk90/ --cluster_prefix MicroEuk90_ --no_sequences_and_header -f table -a easy-linclust -t 90 -p ${N_JOBS} --mmseqs2_options='--split-memory-limit 48G
```
##### MicroEuk50
```
conda activate VEBA-cluster_env

FASTA=mmseqs2_output/MicroEuk90/intermediate/mmseqs2_rep_seq.fasta.gz
mmseqs2_wrapper.py -i ${FASTA} -o mmseqs2_output/MicroEuk50/ --cluster_prefix MicroEuk50_ --no_sequences_and_header -f table -a easy-linclust -t 50 -p ${N_JOBS} --mmseqs2_options='--split-memory-limit 48G
```

#### Building the MMSEQS2 databases for use with MetaEuk

If you were to follow these instructions to build the database from scratch (not advised), you would do this: 

```
source activate VEBA-database_env

# MicroEuk100
mmseqs createdb --compressed 1 MicroEuk100.faa.gz metaeuk_databases/MicroEuk100

# MicroEuk90
mmseqs createdb --compressed 1 mmseqs2_output/MicroEuk90/intermediate/mmseqs2_rep_seq.fasta.gz metaeuk_databases/MicroEuk90

# MicroEuk50
mmseqs createdb --compressed 1 mmseqs2_output/MicroEuk50/intermediate/mmseqs2_rep_seq.fasta.gz metaeuk_databases/MicroEuk50
```

However, I provided a shortcut that obviates you needing to recluster using the precomputed cluster  files downloaded with the `.tar.gz` archive on Zenodo: 

```
# MicroEuk100
mmseqs createdb --compressed 1 MicroEuk100.faa.gz metaeuk_databases/MicroEuk100

# MicroEuk90
gzip -c -c MicroEuk90_clusters.tsv.gz | cut -f1 | sort -u > MicroEuk90.list
seqkit grep -f MicroEuk90.list MicroEuk100.faa.gz | gzip > MicroEuk90.faa
mmseqs createdb --compressed 1 MicroEuk90.faa metaeuk_databases/MicroEuk90

# MicroEuk50
gzip -c -c MicroEuk50_clusters.tsv.gz | cut -f1 | sort -u > MicroEuk50.list
seqkit grep -f MicroEuk50.list MicroEuk90.faa | gzip > MicroEuk50.faa
mmseqs createdb --compressed 1 MicroEuk50.faa metaeuk_databases/MicroEuk50
```

If you want to remove the fasta files: 

```
rm -f MicroEuk100.faa.gz MicroEuk90.faa MicroEuk50.faa
```

#### Building the pickle files used for VEBA's eukaryotic taxonomy classification

Last step is to build the Python pickle objects that are used for eukaryotic taxonomy classification in VEBA. 

```
# target_to_source.dict.pkl.gz {id_protein_hash:'id_source_organism'}
build_target_to_source_dictionary.py -i identifier_mapping.proteins.tsv.gz -o target_to_source.dict.pkl.gz

# source_to_lineage.dict.pkl.gz {id_source_organism:'c_class;o__[order];f__[family];g__[genus];s__[species]}'
build_source_to_lineage_dictionary.py -i source_lineage.tsv.gz -o source_to_lineage.dict.pkl.gz
```


#### Path to MicroEuk_v4:

* Remove redundancy in `source_taxonomy.tsv.gz` (MicroEuk_v3.1)
* Add missing lineage to source organisms when possible (e.g., `MicG_I_3`) (MicroEuk_v3.1)
* Add `TOPAZ` genomes (MicroEuk_v4)
