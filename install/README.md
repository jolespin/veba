<a name="readme-top"></a>

# Modules
[![Schematic](../images/Schematic.png)](../images/Schematic.pdf)

|   Status         |   Module                |   Environment                   |   Executable               |   Resources    |   Recommended Threads  |   Description                                                                                                      |
|------------------|-------------------------|---------------------------------|----------------------------|----------------|------------------------|--------------------------------------------------------------------------------------------------------------------|
|   Stable         |   [preprocess](#preprocesspy)            |   VEBA-preprocess_env           |   preprocess.py            |   4GB-16GB     |   4                    |   Fastq quality trimming, adapter removal, decontamination, and read statistics calculations (Short Reads)         |
|   Stable         |   preproces-long        |   VEBA-preprocess_env           |   preproces-long.py        |   4GB-16GB     |   4                    |   Fastq quality trimming, adapter removal, decontamination, and read statistics calculations (Long Reads)          |
|   Stable         |   assembly              |   VEBA-assembly_env             |   assembly.py              |   32GB-128GB+  |   16                   |   Assemble short reads, align reads to assembly, and count mapped reads                                            |
|   Stable         |   assembly-long         |   VEBA-assembly_env             |   assembly-long.py         |   32GB-128GB+  |   16                   |   Assemble long reads, align reads to assembly, and count mapped reads                                             |
|   Stable         |   coverage              |   VEBA-assembly_env             |   coverage.py              |   24GB         |   16                   |   Align short reads to (concatenated) reference and counts mapped reads                                            |
|   Stable         |   coverage-long         |   VEBA-assembly_env             |   coverage-long.py         |   24GB         |   16                   |   Align long reads to (concatenated) reference and counts mapped reads                                             |
|   Stable         |   binning-prokaryotic   |   VEBA-binning-prokaryotic_env  |   binning-prokaryotic.py   |   16GB         |   4                    |   Iterative consensus binning for recovering prokaryotic genomes with lineage-specific quality assessment          |
|   Stable         |   binning-eukaryotic    |   VEBA-binning-eukaryotic_env   |   binning-eukaryotic.py    |   128GB        |   4                    |   Binning for recovering eukaryotic genomes with exon-aware gene modeling and lineage-specific quality assessment  |
|   Stable         |   binning-viral         |   VEBA-binning-viral_env        |   binning-viral.py         |   16GB         |   4                    |   Detection of viral genomes and quality assessment                                                                |
|   Stable         |   classify-prokaryotic  |   VEBA-classify_env             |   classify-prokaryotic.py  |   64GB         |   32                   |   Taxonomic classification of prokaryotic genomes                                                                  |
|   Stable         |   classify-eukaryotic   |   VEBA-classify_env             |   classify-eukaryotic.py   |   32GB         |   1                    |   Taxonomic classification of eukaryotic genomes                                                                   |
|   Stable         |   classify-viral        |   VEBA-classify_env             |   classify-viral.py        |   16GB         |   4                    |   Taxonomic classification of viral genomes                                                                        |
|   Stable         |   cluster               |   VEBA-cluster_env              |   cluster.py               |   32GB+        |   32                   |   Species-level clustering of genomes and lineage-specific orthogroup detection                                    |
|   Stable         |   annotate              |   VEBA-annotate_env             |   annotate.py              |   64GB         |   32                   |   Annotates translated gene calls                                                      |
|   Stable         |   phylogeny             |   VEBA-phylogeny_env            |   phylogeny.py             |   16GB+        |   32                   |   Constructs phylogenetic trees given a marker set                                                                 |
|   Stable         |   index                 |   VEBA-mapping_env              |   index.py                 |   16GB         |   4                    |   Builds local or global index for alignment to genomes                                                            |
|   Stable         |   mapping               |   VEBA-mapping_env              |   mapping.py               |   16GB         |   4                    |   Aligns reads to local or global index of genomes                                                                 |
|   Stable         |   biosynthetic          |   VEBA-biosynthetic_env         |   biosynthetic.py          |   16GB         |   16                   |   Identify biosynthetic gene clusters in prokaryotes and fungi                                                     |
|   Stable         |   profile-taxonomy      |   VEBA-profile_env              |   profile-taxonomy.py      |   8-24GB       |   4                    |   Taxonomic profiling of *de novo* genomes                                                                           |
|   Stable         |   profile-pathway       |   VEBA-profile_env              |   profile-pathway.py       |   16-32GB      |   4                    |   Pathway profiling of *de novo* genomes                                                                             |
|   Deprecated     |   assembly-sequential   |   VEBA-assembly_env             |   assembly-sequential.py   |   32GB-128GB+  |   16                   |   Assemble metagenomes sequentially                                                                                |
|   Developmental  |   amplicon              |   VEBA-amplicon_env             |   amplicon.py              |   96GB         |   16                   |   Automated read trim position detection, DADA2 ASV detection, taxonomic classification, and file conversion       |

<p align="right"><a href="#readme-top">^__^</a></p>

______________________

### Stable

#### *preprocess.py*
**Fastq quality trimming, adapter removal, decontamination, and read statistics calculations**

The preprocess module is a wrapper around our [`fastq_preprocessor`](https://github.com/jolespin/fastq_preprocessor) which is a modernized reimplementation of [`KneadData`](https://github.com/biobakery/kneaddata) that relies on fastp for ultra-fast automated adapter removal and quality trimming. Pairing of the trimmed reads is assessed and corrected using `BBTools’ repair.sh`. If the user provides a contamination database (e.g., the human reference genome), then trimmed reads are aligned using `Bowtie2` and reads that do not map to the contamination database are stored. If the `--retain_contaminated_reads` flag is used then the contaminated reads are stored as well. Similarly, if a k-mer reference database is provided (e.g., ribosomal k-mers) then the trimmed or decontaminated reads are aligned against the reference database using `BBTools’ bbduk.sh` with an option for storing. By default, the none of the contaminated or k-mer analyzed reads are stored but regardless of the choice for retaining reads, the read sets are quantified using `seqkit` for accounting purposes (e.g., % contamination or % ribosomal). All sequences included were downloaded using `Kingfisher`, included in the preprocess environment, which a fast and flexible program for the procurement of sequencing files and their annotations from public data sources including ENA, NCBI SRA, Amazon AWS, and Google Cloud.

**Conda Environment:** `conda activate VEBA-preprocess_env`

```
usage: preprocess.py -1 <reads_1.fq> -2 <reads_2.fq> -n <name> -o <output_directory> |Optional| -x <reference_index> -k <kmer_database>

    Wrapper around github.com/jolespin/fastq_preprocessor
    Running: preprocess.py v2023.5.8 via Python v3.9.15 | /expanse/projects/jcl110/anaconda3/envs/VEBA-preprocess_env/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -1 FORWARD_READS, --forward_reads FORWARD_READS
                        path/to/reads_1.fastq
  -2 REVERSE_READS, --reverse_reads REVERSE_READS
                        path/to/reads_2.fastq
  -n NAME, --name NAME  Name of sample
  -o PROJECT_DIRECTORY, --project_directory PROJECT_DIRECTORY
                        path/to/project_directory [Default: veba_output/preprocess]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv. Must have at least 2 columns [name, executable] [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint
  -v, --version         show program's version number and exit

Fastp arguments:
  -m MINIMUM_READ_LENGTH, --minimum_read_length MINIMUM_READ_LENGTH
                        Fastp | Minimum read length [Default: 75]
  -a ADAPTERS, --adapters ADAPTERS
                        Fastp | path/to/adapters.fasta [Default: detect]
  --fastp_options FASTP_OPTIONS
                        Fastp | More options (e.g. --arg 1 ) [Default: '']

Bowtie2 arguments:
  -x CONTAMINATION_INDEX, --contamination_index CONTAMINATION_INDEX
                        Bowtie2 | path/to/contamination_index
                        (e.g., Human T2T CHM13 v2 in $VEBA_DATABASE/Contamination/chm13v2.0/chm13v2.0)
  --retain_trimmed_reads RETAIN_TRIMMED_READS
                        Retain fastp trimmed fastq after decontamination. 0=No, 1=yes [Default: 0]
  --retain_contaminated_reads RETAIN_CONTAMINATED_READS
                        Retain contaminated fastq after decontamination. 0=No, 1=yes [Default: 0]
  --bowtie2_options BOWTIE2_OPTIONS
                        Bowtie2 | More options (e.g. --arg 1 ) [Default: '']
                        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

BBDuk arguments:
  -k KMER_DATABASE, --kmer_database KMER_DATABASE
                        BBDuk | path/to/kmer_database
                        (e.g., Ribokmers in $VEBA_DATABASE/Contamination/kmers/ribokmers.fa.gz)
  --kmer_size KMER_SIZE
                        BBDuk | k-mer size [Default: 31]
  --retain_kmer_hits RETAIN_KMER_HITS
                        Retain reads that map to k-mer database. 0=No, 1=yes [Default: 0]
  --retain_non_kmer_hits RETAIN_NON_KMER_HITS
                        Retain reads that do not map to k-mer database. 0=No, 1=yes [Default: 0]
  --bbduk_options BBDUK_OPTIONS
                        BBDuk | More options (e.g., --arg 1) [Default: '']
                      
```

**Output:**

* cleaned_1.fastq.gz - Cleaned and trimmed fastq file (forward)
* cleaned_2.fastq.gz - Cleaned and trimmed fastq file (reverse)
* seqkit_stats.concatenated.tsv - Concatenated read statistics for all intermediate steps (e.g., fastp, bowtie2 removal of contaminated reads if provided, bbduk.sh removal of contaminated reads if provided)

<p align="right"><a href="#readme-top">^__^</a></p>

#### *assembly.py*
**Assemble reads, align reads to assembly, and count mapped reads**

The assembly module optimizes the output for typical metagenomics workflows. In particular, the module does the following: 1) assembles reads using either `metaSPAdes`, `SPAdes`, `rnaSPAdes`, any of the other task-specific assemblers installed with the `SPAdes` package or `MEGAHIT`; 2) builds a `Bowtie2` index for the scaffolds.fasta (or transcripts.fasta if `rnaSPAdes` is used); 3) aligns the reads using `Bowtie2` to the assembly; 4) pipes the alignment file into `Samtools` to produce a sorted BAM file (necessary for many coverage applications); 5) counts the reads mapping to each scaffold via `featureCounts`; and 6) `seqkit` for useful assembly statistics such as N50, number of scaffolds, and total assembly size. This module automates many critical yet overlooked workflows dealing with assemblies. 

**Conda Environment**: `conda activate VEBA-assembly_env`

```
usage: assembly.py -1 <forward_reads.fq> -2 <reverse_reads.fq> -n <name> -o <output_directory>

    Running: assembly.py v2023.5.15 via Python v3.9.15 | /expanse/projects/jcl110/anaconda3/envs/VEBA-preprocess_env/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -1 FORWARD_READS, --forward_reads FORWARD_READS
                        path/to/forward_reads.fq
  -2 REVERSE_READS, --reverse_reads REVERSE_READS
                        path/to/reverse_reads.fq
  -n NAME, --name NAME  Name of sample
  -o PROJECT_DIRECTORY, --project_directory PROJECT_DIRECTORY
                        path/to/project_directory [Default: veba_output/assembly]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit
  --tmpdir TMPDIR       Set temporary directory

Assembler arguments:
  -P PROGRAM, --program PROGRAM
                        Assembler |  {spades.py, metaspades.py, rnaspades.py, megahit, metaplasmidspades.py, plasmidspades.py, coronaspades.py}} [Default: 'metaspades.py']
  -s SCAFFOLD_PREFIX, --scaffold_prefix SCAFFOLD_PREFIX
                        Assembler |  Special options:  Use NAME to use --name.  Use NONE to not include a prefix. [Default: 'NAME__']
  -m MINIMUM_CONTIG_LENGTH, --minimum_contig_length MINIMUM_CONTIG_LENGTH
                        Minimum contig length.  Should be lenient here because longer thresholds can be used for binning downstream. Recommended for metagenomes to use 1000 here. [Default: 1]
  --assembler_options ASSEMBLER_OPTIONS
                        Assembler options for SPAdes-based programs and MEGAHIT (e.g. --arg 1 ) [Default: '']

SPAdes arguments:
  --spades_memory SPADES_MEMORY
                        SPAdes | RAM limit in Gb (terminates if exceeded). [Default: 250]

MEGAHIT arguments:
  --megahit_memory MEGAHIT_MEMORY
                        MEGAHIT | RAM limit in Gb (terminates if exceeded). [Default: 0.9]

Bowtie2 arguments:
  --bowtie2_index_options BOWTIE2_INDEX_OPTIONS
                        bowtie2-build | More options (e.g. --arg 1 ) [Default: '']
  --bowtie2_options BOWTIE2_OPTIONS
                        bowtie2 | More options (e.g. --arg 1 ) [Default: '']

featureCounts arguments:
  --featurecounts_options FEATURECOUNTS_OPTIONS
                        featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/
                        
```

**Output:**

* featurecounts.tsv.gz - featureCounts output for contig-level counts
* mapped.sorted.bam - Sorted BAM
* mapped.sorted.bam.bai - Sorted BAM index
* scaffolds.fasta - Assembly scaffolds (preferred over contigs by SPAdes documentation)
* scaffolds.fasta.\*.bt2 - Bowtie2 index of scaffolds
* scaffolds.fasta.saf - SAF formatted file for contig-level counts with featureCounts
* seqkit_stats.tsv.gz - Assembly statistics

<p align="right"><a href="#readme-top">^__^</a></p>

#### *coverage.py*
**Align reads to (concatenated) reference and counts mapped reads**

The coverage module further optimizes the output for typical metagenomics workflows. In particular, the module does the following: 1) filters contigs based on a size filter (default 1500 bp); 2) builds a `Bowtie2` index for the reference.fasta; 3) aligns the reads from all provided samples using `Bowtie2` to the assembly; 4) pipes the alignment file into `Samtools` to produce a sorted BAM file; 5) counts the reads mapping to each scaffold via `featureCounts`; and 6) `seqkit` for useful assembly statistics such as N50, number of scaffolds, and total assembly size. The preferred usage for this module is after prokaryotic, eukaryotic, and viral binning has been performed and the unbinned contigs are merged into a single coassembly used as input.  The outputs of this module are expected to be used as a final pass through prokaryotic and eukaryotic binning modules.  

**Conda Environment**: `conda activate VEBA-assembly_env`


```
usage: coverage.py -f <reference.fasta> -r <reads.tsv> -o <output_directory>

    Running: coverage.py v2023.5.16 via Python v3.9.15 | /expanse/projects/jcl110/anaconda3/envs/VEBA-preprocess_env/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -f FASTA, --fasta FASTA
                        path/to/reference.fasta. Recommended usage is for merging unbinned contigs. [Required]
  -r READS, --reads READS
                        path/to/reads_table.tsv with the following format: [id\_sample]<tab>[path/to/r1.fastq.gz]<tab>[path/to/r2.fastq.gz], No header
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        path/to/project_directory [Default: veba_output/assembly/multisample]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit
  --tmpdir TMPDIR       Set temporary directory

SeqKit seq arguments:
  -m MINIMUM_CONTIG_LENGTH, --minimum_contig_length MINIMUM_CONTIG_LENGTH
                        seqkit seq | Minimum contig length [Default: 1500]
  --seqkit_seq_options SEQKIT_SEQ_OPTIONS
                        seqkit seq | More options (e.g. --arg 1 ) [Default: '']

Bowtie2 arguments:
  --bowtie2_index_options BOWTIE2_INDEX_OPTIONS
                        bowtie2-build | More options (e.g. --arg 1 ) [Default: '']
  --one_task_per_cpu    Use GNU parallel to run GNU parallel with 1 task per CPU.  Useful if all samples are roughly the same size but inefficient if depth varies.
  --bowtie2_options BOWTIE2_OPTIONS
                        bowtie2 | More options (e.g. --arg 1 ) [Default: '']

featureCounts arguments:
  --featurecounts_options FEATURECOUNTS_OPTIONS
                        featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/
                        
```

**Output:**

* featurecounts.tsv.gz - featureCounts counts table of all samples
* [sample_id]/mapped.sorted.bam - Sorted BAM file under a subdirectory for each sample
* reference.fasta - Reference fasta (typically this would be the pseudo-coassembly of unbinned contigs)
* reference.fasta.saf - SAF formatted file for contig-level counts with featureCounts
* seqkit_stats.tsv - Assembly statistics

<p align="right"><a href="#readme-top">^__^</a></p>

#### *binning-prokaryotic.py*
**Iterative consensus binning for recovering prokaryotic genomes with lineage-specific quality assessment**

The prokaryotic binning module implements a novel iterative consensus binning procedure that uses `CoverM`  for fast coverage calculations, multiple binning algorithms (MaxBin2 (marker set = 107); `MaxBin2` (marker set = 40); `MetaBat2`; and `CONCOCT`, consensus dereplication and aggregate binning with `DAS Tool`, the consensus domain wrapper for `Tiara` for removing eukaryotes at the MAG level, and `CheckM2` for quality assessment where poor quality MAGs are removed (e.g., completeness < 50% and/or contamination ≥ 10) which has direct support for candidate phyla radiation (CPR). The novelty of this procedure is that the unbinned contigs are stored and fed back into the input of the binning procedure using a separate random seed state allowing for an exhaustive, yet effective, approach in extracting high quality and difficult to bin genomes; number of iterations specified by `--n_iter` option. Gene calls are performed using `Pyrodigal` and the gene models (GFF3 Format) are modified to include gene and contig identifiers for use with downstream feature counting software. `BARRNAP` and `tRNAscan-SE` are used for rRNA and tRNA detetion, respectively. Lastly, genome assembly a gene statistics such as GC, N50, number of scaffolds, and genome size are calculated using `seqkit`.  MAG naming scheme for prokaryotes follows `[SAMPLE]__[ALGORITHM]__P.[ITERATION]__[NAME]` (e.g., `SRR17458623__METABAT2__P.1__bin.1`)

**⚠️Notes:** 

* If you have a lot of samples and a lot of contigs then use the `--skip_maxbin2` flag because it takes MUCH longer to run.  For the *Plastisphere* it was going to take 40 hours per `MaxBin2` run (there are 2 `MaxBin2` runs) per iteration.  `Metabat2` and `CONCOCT` can do the heavy lifting much faster and often with better results so it's recommended to skip `MaxBin2` for larger datasets.


**Conda Environment**: `conda activate VEBA-binning-prokaryotic_env`


```
usage: binning-prokaryotic.py -f <scaffolds.fasta> -b <mapped.sorted.bam> -n <name> -o <output_directory>

    Running: binning-prokaryotic.py v2023.7.7 via Python v3.9.15 | /expanse/projects/jcl110/anaconda3/envs/VEBA-preprocess_env/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -f FASTA, --fasta FASTA
                        path/to/scaffolds.fasta
  -b BAM [BAM ...], --bam BAM [BAM ...]
                        path/to/mapped.sorted.bam files separated by spaces.
  -n NAME, --name NAME  Name of sample
  -o PROJECT_DIRECTORY, --project_directory PROJECT_DIRECTORY
                        path/to/project_directory [Default: veba_output/binning/prokaryotic/]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  -I N_ITER, --n_iter N_ITER
                        Number of iterations to run binning. Use -1 for brute force (-1 isn't ready yet) [Default: 3]
  --random_state RANDOM_STATE
                        Use -1 for completely random. Use 0 for consecutive random states.  Use any other positive integer for the same random state for all iterations [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit
  --tmpdir TMPDIR       Set temporary directory

Database arguments:
  --veba_database VEBA_DATABASE
                        VEBA database location.  [Default: $VEBA_DATABASE environment variable]

Binning arguments:
  -m MINIMUM_CONTIG_LENGTH, --minimum_contig_length MINIMUM_CONTIG_LENGTH
                        Minimum contig length.  Anything under 2500 will default to 2500 for MetaBat2 [Default: 1500]
  -s MINIMUM_GENOME_LENGTH, --minimum_genome_length MINIMUM_GENOME_LENGTH
                        Minimum genome length.  [Default: 150000]
  --concoct_fragment_length CONCOCT_FRAGMENT_LENGTH
                        CONCOCT | Fragment length [Default: 10000]
  --concoct_overlap_length CONCOCT_OVERLAP_LENGTH
                        CONCOCT | Fragment overlap length [Default: 0]
  --skip_maxbin2        MaxBin2 | Skip MaxBin2. Useful for large datasets
  --skip_concoct        CONCOCT | Skip CONCOCT. Useful when there's a lot of samples
  --maxbin2_options MAXBIN2_OPTIONS
                        MaxBin2 | More options (e.g. --arg 1 ) [Default: ''] | https://sourceforge.net/projects/maxbin/
  --metabat2_options METABAT2_OPTIONS
                        MetaBat2 | More options (e.g. --arg 1 ) [Default: ''] | https://bitbucket.org/berkeleylab/metabat/src/master/
  --concoct_options CONCOCT_OPTIONS
                        CONCOCT | More options (e.g. --arg 1 ) [Default: '']

Gene model arguments:
  --pyrodigal_minimum_gene_length PYRODIGAL_MINIMUM_GENE_LENGTH
                        Pyrodigal | Minimum gene length [Default: 90]
  --pyrodigal_minimum_edge_gene_length PYRODIGAL_MINIMUM_EDGE_GENE_LENGTH
                        Pyrodigal | Minimum edge gene length [Default: 60]
  --pyrodigal_maximum_gene_overlap_length PYRODIGAL_MAXIMUM_GENE_OVERLAP_LENGTH
                        Pyrodigal | Maximum gene overlap length [Default: 60]
  --pyrodigal_genetic_code PYRODIGAL_GENETIC_CODE
                        Pyrodigal -g translation table [Default: 11]

Evaluation arguments:
  --dastool_searchengine DASTOOL_SEARCHENGINE
                        DAS_Tool searchengine. [Default: diamond] | https://github.com/cmks/DAS_Tool
  --dastool_minimum_score DASTOOL_MINIMUM_SCORE
                        DAS_Tool score_threshold. Score threshold until selection algorithm will keep selecting bins. This is set to a relaxed setting because CheckM2 is run post hoc. [Default: 0.1] | https://github.com/cmks/DAS_Tool
  --dastool_options DASTOOL_OPTIONS
                        DAS_Tool | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/cmks/DAS_Tool
  --checkm2_completeness CHECKM2_COMPLETENESS
                        CheckM2 completeness threshold [Default: 50.0]
  --checkm2_contamination CHECKM2_CONTAMINATION
                        CheckM2 contamination threshold [Default: 10.0]
  --checkm2_options CHECKM2_OPTIONS
                        CheckM lineage_wf | More options (e.g. --arg 1 ) [Default: '']

barrnap arguments:
  --barrnap_length_cutoff BARRNAP_LENGTH_CUTOFF
                        barrnap | Proportional length threshold to label as partial [Default: 0.8]
  --barrnap_reject BARRNAP_REJECT
                        barrnap | Proportional length threshold to reject prediction [Default: 0.25]
  --barrnap_evalue BARRNAP_EVALUE
                        barrnap | Similarity e-value cut-off [Default: 1e-6]

tRNAscan-SE arguments:
  --trnascan_options TRNASCAN_OPTIONS
                        tRNAscan-SE | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/UCSC-LoweLab/tRNAscan-SE

featureCounts arguments:
  --long_reads          featureCounts | Use this if long reads are being used
  --featurecounts_options FEATURECOUNTS_OPTIONS
                        featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/

Domain classification arguments:
  --logit_transform LOGIT_TRANSFORM
                         Transformation for consensus_domain_classification: {softmax, tss} [Default: softmax]
  --tiara_minimum_length TIARA_MINIMUM_LENGTH
                        Tiara | Minimum contig length. Anything lower than 3000 is not recommended. [Default: 3000]
  --tiara_options TIARA_OPTIONS
                        Tiara | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/ibe-uw/tiara
                        
```

**Output:**

* binned.list - List of binned contigs
* bins.list - List of MAG identifiers
* checkm2_results.filtered.tsv - Filtered CheckM2 output
* featurecounts.orfs.tsv.gz - ORF-level counts table
* genome_statistics.tsv - Genome assembly statistics
* gene_statistics.cds.tsv - Gene sequence statistics (CDS)
* gene_statistics.rRNA.tsv - Gene sequence statistics (rRNA)
* gene_statistics.tRNA.tsv - Gene sequence statistics (tRNA)
* genomes/ - MAG subdirectory
* genomes/[id\_genome].fa - MAG assembly fasta
* genomes/[id\_genome].faa - MAG protein fasta
* genomes/[id\_genome].ffn - MAG CDS fasta
* genomes/[id\_genome].gff - MAG gene models for assembly, CDS, rRNA, and tRNA
* genomes/[id\_genome].rRNA - MAG rRNA fasta
* genomes/[id\_genome].tRNA - MAG tRNA fasta
* genomes/identifier\_mapping.tsv - Identifier mapping between [id\_orf, id\_contig, id_mag]
* scaffolds\_to\_bins.tsv - Identifier mapping between [id\_contig, id_mag]
* unbinned.fasta - Fasta of unbinned contigs that have passed length thresholding
* unbinned.list - List of unbinned contigs

<p align="right"><a href="#readme-top">^__^</a></p>

#### *binning-eukaryotic.py*
**Binning for recovering eukaryotic genomes with exon-aware gene modeling and lineage-specific quality assessment**
The eukaryotic binning module uses several checks and state-of-the-art software to ensure high quality genomes.  In particular, non-prokaryotic-biased binning algorithms `MetaBat2` [default] (coverage calculated with `CoverM`) or `CONCOCT` (coverage calculated using `CONCOCT` scripts) is used for binning out genomes followed by a genome size filter (2,000,000 bp is the default). The preliminary bins are run through the consensus domain wrapper for `Tiara` to predict eukaryotic MAGs. Contigs from the eukaryotic MAGs are input into `MetaEuk` easy-predict workflow using our custom consensus eukaryotic database. Although `MetaEuk` is a high-quality software suite, the identifiers from `MetaEuk` are very complex, long, and contain characters that are often problematic for downstream applications including parsing, file naming systems, and certain programs with simplified identifier requirements such as `Anvi’o`. In addition, the gene model GFF files are not intuitive, compatible with `P(y)rodigal(-gv)` GFF files or `featureCounts` without major modification. Therefore, we developed an essential wrapper for MetaEuk that simplifies identifiers (i.e., `[ContigID]_[GeneStart]:[GeneEnd]([strand])`), ensuring no duplicates are produced, creates a GFF file that can be concatenated with the `Pyrodigal` GFF file for use with `featureCounts`, and several identifier mapping tables to seamless convert between original and modified identifiers. `BARRNAP` and `tRNAscan-SE` are used for rRNA and tRNA detetion, respectively. `Tiara` is used to identfy mitochondrion and plastid contigs to run separate workflows for modeling genes (CDS, rRNA, and tRNA).  Lineage-specific genome quality estimation is performed using `BUSCO` where poor quality MAGs are removed (e.g., completeness < 50%, contamination > 10).  Gene counts are computed using `featureCounts` at the gene level.  Lastly, genome and gene statistics such as N50, number of scaffolds, and genome size are calculated using `seqkit`.  Iterative binning is not currently available as no consensus binning tool is available therefore iterative binning would result in diminishing returns. MAG naming scheme for eukaryotes follows `[SAMPLE]__[ALGORITHM]__E.[ITERATION]__[NAME]` (e.g., `ERR2002407__METABAT2__E.1__bin.2`).

**Conda Environment**: `conda activate VEBA-binning-eukaryotic_env`


```
usage: binning-eukaryotic.py -f <scaffolds.fasta> -b <mapped.sorted.bam> -n <name> -o <output_directory>

    Running: binning-eukaryotic.py v2023.7.6 via Python v3.9.15 | /expanse/projects/jcl110/anaconda3/envs/VEBA-preprocess_env/bin/python

optional arguments:
  -h, --help            show this help message and exit

I/O arguments:
  -f FASTA, --fasta FASTA
                        path/to/scaffolds.fasta
  -b BAM [BAM ...], --bam BAM [BAM ...]
                        path/to/mapped.sorted.bam files separated by spaces.
  -n NAME, --name NAME  Name of sample
  -l CONTIG_IDENTIFIERS, --contig_identifiers CONTIG_IDENTIFIERS
                        path/to/contigs.list [Optional]
  -o PROJECT_DIRECTORY, --project_directory PROJECT_DIRECTORY
                        path/to/project_directory [Default: veba_output/binning/eukaryotic]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit

Database arguments:
  --veba_database VEBA_DATABASE
                        VEBA database location.  [Default: $VEBA_DATABASE environment variable]

Binning arguments:
  -a ALGORITHM, --algorithm ALGORITHM
                        Binning algorithm: {concoct, metabat2}  [Default: metabat2]
  -m MINIMUM_CONTIG_LENGTH, --minimum_contig_length MINIMUM_CONTIG_LENGTH
                        Minimum contig length.  [Default: 1500]
  -s MINIMUM_GENOME_LENGTH, --minimum_genome_length MINIMUM_GENOME_LENGTH
                        Minimum genome length.  [Default: 2000000]
  --concoct_fragment_length CONCOCT_FRAGMENT_LENGTH
                        CONCOCT | Fragment length [Default: 10000]
  --concoct_overlap_length CONCOCT_OVERLAP_LENGTH
                        CONCOCT | Fragment overlap length [Default: 0]
  --concoct_options CONCOCT_OPTIONS
                        CONCOCT | More options (e.g. --arg 1 ) [Default: '']
  --metabat2_options METABAT2_OPTIONS
                        MetaBat2 | More options (e.g. --arg 1 ) [Default: ''] | https://bitbucket.org/berkeleylab/metabat/src/master/

Domain classification arguments:
  --logit_transform LOGIT_TRANSFORM
                         Transformation for consensus_domain_classification: {softmax, tss} [Default: softmax
  --tiara_minimum_length TIARA_MINIMUM_LENGTH
                        Tiara | Minimum contig length. Anything lower than 3000 is not recommended. [Default: 3000]
  --tiara_options TIARA_OPTIONS
                        Tiara | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/ibe-uw/tiara

MetaEuk arguments:
  --metaeuk_sensitivity METAEUK_SENSITIVITY
                        MetaEuk | Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive  [Default: 4.0]
  --metaeuk_evalue METAEUK_EVALUE
                        MetaEuk | List matches below this E-value (range 0.0-inf) [Default: 0.01]
  --metaeuk_options METAEUK_OPTIONS
                        MetaEuk | More options (e.g. --arg 1 ) [Default: ''] https://github.com/soedinglab/metaeuk

Pyrodigal arguments (Mitochondria):
  --pyrodigal_minimum_gene_length PYRODIGAL_MINIMUM_GENE_LENGTH
                        Pyrodigal | Minimum gene length [Default: 90]
  --pyrodigal_minimum_edge_gene_length PYRODIGAL_MINIMUM_EDGE_GENE_LENGTH
                        Pyrodigal | Minimum edge gene length [Default: 60]
  --pyrodigal_maximum_gene_overlap_length PYRODIGAL_MAXIMUM_GENE_OVERLAP_LENGTH
                        Pyrodigal | Maximum gene overlap length [Default: 60]
  --pyrodigal_mitochondrial_genetic_code PYRODIGAL_MITOCHONDRIAL_GENETIC_CODE
                        Pyrodigal -g translation table (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi/) [Default: 4] (The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code))
  --pyrodigal_plastid_genetic_code PYRODIGAL_PLASTID_GENETIC_CODE
                        Pyrodigal -g translation table (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi/) [Default: 11] (The Bacterial, Archaeal and Plant Plastid Code))

BUSCO arguments:
  --busco_completeness BUSCO_COMPLETENESS
                        BUSCO completeness [Default: 50.0]
  --busco_contamination BUSCO_CONTAMINATION
                        BUSCO contamination [Default: 10.0]
  --busco_evalue BUSCO_EVALUE
                        BUSCO | E-value cutoff for BLAST searches. Allowed formats, 0.001 or 1e-03 [Default: 1e-03]

barrnap arguments:
  --barrnap_length_cutoff BARRNAP_LENGTH_CUTOFF
                        barrnap | Proportional length threshold to label as partial [Default: 0.8]
  --barrnap_reject BARRNAP_REJECT
                        barrnap | Proportional length threshold to reject prediction [Default: 0.25]
  --barrnap_evalue BARRNAP_EVALUE
                        barrnap | Similarity e-value cut-off [Default: 1e-6]

tRNAscan-SE arguments:
  --trnascan_nuclear_options TRNASCAN_NUCLEAR_OPTIONS
                        tRNAscan-SE | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/UCSC-LoweLab/tRNAscan-SE
  --trnascan_mitochondrial_searchmode TRNASCAN_MITOCHONDRIAL_SEARCHMODE
                        tRNAscan-SE | Search mode [Default: '-O'] | Current best option according to developer: https://github.com/UCSC-LoweLab/tRNAscan-SE/issues/24
  --trnascan_mitochondrial_options TRNASCAN_MITOCHONDRIAL_OPTIONS
                        tRNAscan-SE | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/UCSC-LoweLab/tRNAscan-SE
  --trnascan_plastid_searchmode TRNASCAN_PLASTID_SEARCHMODE
                        tRNAscan-SE | Search mode [Default: '-O'] | https://github.com/UCSC-LoweLab/tRNAscan-SE
  --trnascan_plastid_options TRNASCAN_PLASTID_OPTIONS
                        tRNAscan-SE | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/UCSC-LoweLab/tRNAscan-SE

featureCounts arguments:
  --long_reads          featureCounts | Use this if long reads are being used
  --featurecounts_options FEATURECOUNTS_OPTIONS
                        featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/

```

**Output:**

* binned.list - List of binned contigs
* bins.list - List of MAG identifiers
* busco_results.filtered.tsv - Filtered BUSCO output
* featurecounts.orfs.tsv.gz - ORF-level counts table
* genome_statistics.tsv - Genome assembly statistics
* gene_statistics.cds.tsv - Gene sequence statistics (CDS)
* gene_statistics.rRNA.tsv - Gene sequence statistics (rRNA)
* gene_statistics.tRNA.tsv - Gene sequence statistics (tRNA)
* genomes/ - MAG subdirectory
* genomes/[id\_genome].fa - MAG assembly fasta
* genomes/[id\_genome].faa - MAG protein fasta
* genomes/[id\_genome].ffn - MAG CDS fasta
* genomes/[id\_genome].gff - MAG gene models for assembly, CDS, rRNA, and tRNA
* genomes/[id\_genome].rRNA - MAG rRNA fasta
* genomes/[id\_genome].tRNA - MAG tRNA fasta
* genomes/[id\_genome].seqtype.tsv - Identifier mapping between [id\_contig, sequence_type] {nuclear, mitochondrion, plastid}
* genomes/identifier\_mapping.tsv - Identifier mapping between [id\_orf, id\_contig, id_mag]
* identifier\_mapping.metaeuk.tsv - Identifier mapping between original MetaEuk identifiers and modified identifiers.  Includes fully parsed MetaEuk identifiers.
* scaffolds\_to\_bins.tsv - Identifier mapping between [id\_contig, id_mag]
* unbinned.fasta - Fasta of unbinned contigs that have passed length thresholding
* unbinned.list - List of unbinned contigs

<p align="right"><a href="#readme-top">^__^</a></p>

#### *binning-viral.py*
**Detection of viral genomes and quality assessment**

Viral binning is performed using either `geNomad` [default] or `VirFinder` to extract potential viral contigs (e.g., P < 0.05). The potential viral contigs are then input into `CheckV`  where quality assessment removes poor quality or low confidence viral predictions. The filtering scheme is based on author recommendations in which a viral contig is considered if it meets the following criteria: 1) number of viral genes ≥ 5 x number of host genes; 2) completeness ≥ 50%; 3) `CheckV` quality is either medium-quality, high-quality, or complete; and 4) `MIUViG` quality is either medium-quality, high-quality, or complete.  Proviruses can be included by using the `--include_proviruses` flag. After poor quality viral contigs are removed, `Prodigal-GV` is used for gene modeling and `seqkit` is used for useful genome statistics.  Iterative binning is not applicable for viral detection as algorithms are executed on a per-contig basis and all viral genomes will be identified on first pass.  MAG naming scheme for viruses follows `[SAMPLE]__[ALGORITHM]__[NAME]` (e.g., `SRR9668957__GENOMAD__Virus.1`).

**Conda Environment**: `conda activate VEBA-binning-viral_env`

```
usage: binning-viral.py -f <scaffolds.fasta> -l <contig_identifiers> -n <name> -o <output_directory> [Requires at least 20GB]

    Running: binning-viral.py v2023.7.7 via Python v3.9.15 | /expanse/projects/jcl110/anaconda3/envs/VEBA-preprocess_env/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -f FASTA, --fasta FASTA
                        path/to/scaffolds.fasta
  -n NAME, --name NAME  Name of sample
  -o PROJECT_DIRECTORY, --project_directory PROJECT_DIRECTORY
                        path/to/project_directory [Default: veba_output/binning/viral]
  -b BAM [BAM ...], --bam BAM [BAM ...]
                        path/to/mapped.sorted.bam files separated by spaces.

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit

Database arguments:
  --veba_database VEBA_DATABASE
                        VEBA database location.  [Default: $VEBA_DATABASE environment variable]

Binning arguments:
  -a ALGORITHM, --algorithm ALGORITHM
                        Binning algorithm to use: {genomad, virfinder}  [Default: genomad]
  -m MINIMUM_CONTIG_LENGTH, --minimum_contig_length MINIMUM_CONTIG_LENGTH
                        Minimum contig length.  [Default: 1500]
  --include_provirus_detection
                        Include provirus viral detection

Gene model arguments:
  --prodigal_genetic_code PRODIGAL_GENETIC_CODE
                        Prodigal-GV -g translation table (https://github.com/apcamargo/prodigal-gv) [Default: 11]

geNomad arguments
Using --relaxed mode by default.  Adjust settings according to the following table: https://portal.nersc.gov/genomad/post_classification_filtering.html#default-parameters-and-presets:
  --genomad_qvalue GENOMAD_QVALUE
                        Maximum accepted false discovery rate. [Default: 1.0; 0.0 < x ≤ 1.0]
  --sensitivity SENSITIVITY
                        MMseqs2 marker search sensitivity. Higher values will annotate more proteins, but the search will be slower and consume more memory. [Default: 4.0; x ≥ 0.0]
  --splits SPLITS       Split the data for the MMseqs2 search. Higher values will reduce memory usage, but will make the search slower. If the MMseqs2 search is failing, try to increase the number of splits. Also used for VirFinder. [Default: 0; x ≥ 0]
  --composition COMPOSITION
                        Method for estimating sample composition. (auto|metagenome|virome) [Default: auto]
  --minimum_score MINIMUM_SCORE
                        Minimum score to flag a sequence as virus or plasmid. By default, the sequence is classified as virus/plasmid if its virus/plasmid score is higher than its chromosome score, regardless of the value. [Default: 0; 0.0 ≤ x ≤ 1.0]
  --minimum_plasmid_marker_enrichment MINIMUM_PLASMID_MARKER_ENRICHMENT
                        Minimum allowed value for the plasmid marker enrichment score, which represents the total enrichment of plasmid markers in the sequence. Sequences with multiple plasmid markers will have higher values than the ones that encode few or no markers.[Default: -100]
  --minimum_virus_marker_enrichment MINIMUM_VIRUS_MARKER_ENRICHMENT
                        Minimum allowed value for the virus marker enrichment score, which represents the total enrichment of plasmid markers in the sequence. Sequences with multiple plasmid markers will have higher values than the ones that encode few or no markers. [Default: -100]
  --minimum_plasmid_hallmarks MINIMUM_PLASMID_HALLMARKS
                        Minimum number of plasmid hallmarks in the identified plasmids.  [Default: 0; x ≥ 0]
  --minimum_virus_hallmarks MINIMUM_VIRUS_HALLMARKS
                        Minimum number of virus hallmarks in the identified viruses.  [Default: 0; x ≥ 0]
  --maximum_universal_single_copy_genes MAXIMUM_UNIVERSAL_SINGLE_COPY_GENES
                        Maximum allowed number of universal single copy genes (USCGs) in a virus or a plasmid. Sequences with more than this number of USCGs will not be classified as viruses or plasmids, regardless of their score.  [Default: 100]
  --genomad_options GENOMAD_OPTIONS
                        geNomad | More options (e.g. --arg 1 ) [Default: '']

VirFinder arguments:
  --virfinder_pvalue VIRFINDER_PVALUE
                        VirFinder statistical test threshold [Default: 0.05]
  --mmseqs2_evalue MMSEQS2_EVALUE
                        Maximum accepted E-value in the MMseqs2 search. Used by genomad annotate when VirFinder is used as binning algorithm [Default: 1e-3]
  --use_qvalue          Use qvalue (FDR) instead of pvalue
  --use_minimal_database_for_taxonomy
                        Use a smaller marker database to annotate proteins. This will make execution faster but sensitivity will be reduced.
  --virfinder_options VIRFINDER_OPTIONS
                        VirFinder | More options (e.g. --arg 1 ) [Default: '']

CheckV arguments:
  --checkv_options CHECKV_OPTIONS
                        CheckV | More options (e.g. --arg 1 ) [Default: '']
  --multiplier_viral_to_host_genes MULTIPLIER_VIRAL_TO_HOST_GENES
                        Minimum number of viral genes [Default: 5]
  --checkv_completeness CHECKV_COMPLETENESS
                        Minimum completeness [Default: 50.0]
  --checkv_quality CHECKV_QUALITY
                        Comma-separated string of acceptable arguments between {High-quality,Medium-quality,Complete} [Default: High-quality,Medium-quality,Complete]
  --miuvig_quality MIUVIG_QUALITY
                        Comma-separated string of acceptable arguments between {High-quality,Medium-quality,Complete} [Default: High-quality,Medium-quality,Complete]

featureCounts arguments:
  --long_reads          featureCounts | Use this if long reads are being used
  --featurecounts_options FEATURECOUNTS_OPTIONS
                        featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/

```

**Output:**

* binned.list - List of binned contigs
* bins.list - List of MAG identifiers
* checkv_results.filtered.tsv - Filtered CheckV output
* featurecounts.orfs.tsv.gz - ORF-level counts table (If --bam file is provided)
* genome_statistics.tsv - Genome assembly statistics
* gene_statistics.cds.tsv - Gene sequence statistics (CDS)
* genomes/ - MAG subdirectory
* genomes/[id\_genome].fa - MAG assembly fasta
* genomes/[id\_genome].faa - MAG protein fasta
* genomes/[id\_genome].ffn - MAG CDS fasta
* genomes/[id\_genome].gff - MAG gene models
* genomes/identifier\_mapping.tsv - Identifier mapping between [id\_orf, id\_contig, id_mag]
* scaffolds\_to\_bins.tsv - Identifier mapping between [id\_contig, id_mag]
* unbinned.fasta - Fasta of unbinned contigs that have passed length thresholding
* unbinned.list - List of unbinned contigs

<p align="right"><a href="#readme-top">^__^</a></p>

#### *classify-prokaryotic.py*
**Taxonomic classification of prokaryotic genomes**

The prokaryotic classification module is a useful wrapper around `GTDB-Tk` which either combines the resulting archaea and bacteria summary tables or runs `GTDB-Tk lineage_wf` from the beginning.  If genome clusters are provided, then it performs consensus lineage classification.  `Krona` plots are generated showing taxonomic levels.

**Conda Environment**: `conda activate VEBA-classify_env`


```
usage: classify-prokaryotic.py -i <prokaryotic_binning_directory>|-g <genomes.list> -o <output_directory>

    Running: classify-prokaryotic.py v2023.6.16 via Python v3.9.15 | /expanse/projects/jcl110/anaconda3/envs/VEBA-preprocess_env/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -i PROKARYOTIC_BINNING_DIRECTORY, --prokaryotic_binning_directory PROKARYOTIC_BINNING_DIRECTORY
                        path/to/prokaryotic_binning_directory [Cannot be used with --genomes]
  -g GENOMES, --genomes GENOMES
                        path/to/genomes.list [Cannot be ued with --prokaryotic_binning_directory]
  -c CLUSTERS, --clusters CLUSTERS
                        path/to/clusters.tsv, Format: [id_mag]<tab>[id_cluster], No header.
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        path/to/output_directory [Default: veba_output/classify/prokaryotic]
  -x EXTENSION, --extension EXTENSION
                        Fasta file extension for genomes if a list is provided [Default: fa]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  --tmpdir TMPDIR       path/to/TMPDIR
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit

Database arguments:
  --veba_database VEBA_DATABASE
                        VEBA database location.  [Default: $VEBA_DATABASE environment variable]

GTDB-Tk arguments:
  --skip_ani_screen     Skip ANI screen [Default: Don't skip ANI screen]
  --gtdbtk_options GTDBTK_OPTIONS
                        GTDB-Tk | classify_wf options (e.g. --arg 1 ) [Default: '']

Consensus genome classification arguments:
  -l LENIENCY, --leniency LENIENCY
                        Leniency parameter. Lower value means more conservative weighting. A value of 1 indiciates no weight bias. A value greater than 1 puts higher weight on higher level taxonomic assignments. A value less than 1 puts lower weights on higher level taxonomic assignments.  [Default: 1.382]

```


**Output:**

* taxonomy.tsv - Prokaryotic genome classification based on GTDBTk
* taxonomy.clusters.tsv - Prokaryotic cluster classification (If --clusters are provided)
* krona.html - Krona plot for various taxonomic levels

<p align="right"><a href="#readme-top">^__^</a></p>
 
                
#### *classify-eukaryotic.py*
**Taxonomic classification of eukaryotic genomes**

The eukaryotic classification module can be performed on de novo genomes or utilize the target field of `MetaEuk` gene identifiers and the taxonomic lineage associated with each source genome.  The default marker set is `eukaryote_odb10` from `BUSCO` but custom marker sets are support along with the inclusion of all genes not just marker genes.  An option to include marker-specific noise cutoff scores is also available using the `--scores_cutoff` parameter which is default behavior with `BUSCO’s eukaryote_odb10` provided noise thresholds.  For each MAG, bitscores are accumulated for each taxonomic level and taxonomy is assigned with leniency specified by the leniency parameter with high leniency resulting higher order taxonomic assignments.  If genome clusters are provided, then it performs consensus lineage classification. `Krona` plots are generated showing taxonomic levels.

**Conda Environment**: `conda activate VEBA-classify_env` 

```
usage: classify-eukaryotic.py -i <eukaryotic_binning_directory>|-g <genomes.list> -o <output_directory>

    Running: classify-eukaryotic.py v2023.6.12 via Python v3.9.15 | /expanse/projects/jcl110/anaconda3/envs/VEBA-preprocess_env/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -i EUKARYOTIC_BINNING_DIRECTORY, --eukaryotic_binning_directory EUKARYOTIC_BINNING_DIRECTORY
                        path/to/eukaryotic_binning_directory [Cannot be used with --genomes]
  -g GENOMES, --genomes GENOMES
                        path/to/genomes.list where each line is a path to a genome.fasta [Cannot be ued with --eukaryotic_binning_directory]
  -c CLUSTERS, --clusters CLUSTERS
                        path/to/clusters.tsv, Format: [id_mag]<tab>[id_cluster], No header.
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        path/to/output_directory [Default: veba_output/classify/eukaryotic]
  -x EXTENSION, --extension EXTENSION
                        path/to/output_directory.  Does not support gzipped. [Default: fa]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit

HMMSearch arguments:
  --hmmsearch_threshold HMMSEARCH_THRESHOLD
                        HMMSearch | Threshold {cut_ga, cut_nc, gut_tc, e} [Default:  e]
  --hmmsearch_evalue HMMSEARCH_EVALUE
                        HMMSearch | E-Value [Default: 10.0]
  --hmmsearch_options HMMSEARCH_OPTIONS
                        HMMSearch | More options (e.g. --arg 1 ) [Default: '']

Database arguments:
  --veba_database VEBA_DATABASE
                        VEBA database location.  [Default: $VEBA_DATABASE environment variable]

Consensus genome arguments:
  -l LENIENCY, --leniency LENIENCY
                        Leniency parameter. Lower value means more conservative weighting. A value of 1 indiciates no weight bias. A value greater than 1 puts higher weight on higher level taxonomic assignments. A value less than 1 puts lower weights on higher level taxonomic assignments.  [Default: 1.382]
  --similarity_threshold SIMILARITY_THRESHOLD
                        Threshold for similarity analysis [Default: 0.8]
  --retain_unannotated RETAIN_UNANNOTATED
                        Consider unannotations (i.e., blank functions) in the scording system [Default: 1]
  --unannotated_weight UNANNOTATED_WEIGHT
                        Weight for unannotations (i.e., blank functions) in the scording system? [Default: 0.382]
  --representative_threshold REPRESENTATIVE_THRESHOLD
                        Score to consider as representative [Default: 0.618]

``` 

**Output:**

* taxonomy.tsv - Eukaryotic genome classification based on microeukaryotic protein database and BUSCO's eukaryota_odb10 marker set
* taxonomy.clusters.tsv - Eukaryotic cluster classification (If --clusters are provided)
* gene-source\_lineage.tsv - Gene source lineage and scores for classifying MAGs [id_gene, id_scaffold, id_mag, id_target, id_source, lineage, bitscore]
* krona.html - Krona plot for various taxonomic levels

<p align="right"><a href="#readme-top">^__^</a></p>

#### *classify-viral.py*
**Taxonomic classification for viral genomes**

Viral classification uses `geNomad's taxonomy` module. 

**Conda Environment**: `conda activate VEBA-classify_env`

```
usage: classify-viral.py -i <viral_binning_directory>|-g <genomes.list>  -o <output_directory>

    Running: classify-viral.py v2023.5.8 via Python v3.9.15 | /expanse/projects/jcl110/anaconda3/envs/VEBA-preprocess_env/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -i VIRAL_BINNING_DIRECTORY, --viral_binning_directory VIRAL_BINNING_DIRECTORY
                        Either: path/to/checkv/quality_summary.tsv or directory of veba_output/binning/viral
  -g GENOMES, --genomes GENOMES
                        path/to/genomes.list where each line is a path to a genome.fasta [Cannot be ued with --viral_binning_directory]
  -c CLUSTERS, --clusters CLUSTERS
                        path/to/clusters.tsv, Format: [id_mag]<tab>[id_cluster], No header.
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        path/to/output_directory [Default: veba_output/classify/viral]
  -x EXTENSION, --extension EXTENSION
                        path/to/output_directory.  Does not support gzipped. [Default: fa]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit

Database arguments:
  --veba_database VEBA_DATABASE
                        VEBA database location.  [Default: $VEBA_DATABASE environment variable]

Consensus genome classification arguments:
  -t THRESHOLD, --threshold THRESHOLD
                        Fraction of classifications for consensus [Default: 0.5]

```

**Output:**

* taxonomy.tsv - Viral genome classification based on geNomad classifications
* taxonomy.clusters.tsv - Viral cluster classification (If --clusters are provided)

<p align="right"><a href="#readme-top">^__^</a></p>

#### *cluster.py*
**Species-level clustering of genomes and lineage-specific orthogroup detection**

To leverage intra-sample genome analysis in an inter-sample analytical paradigm, genome clustering and lineage-specific orthogroup detection is necessary.  The cluster module utilizes 2 separate wrappers for global (inter-sample) and local (intra-sample) clustering.  `FastANI` is used to compute pairwise ANI and these are used to construct a `NetworkX graph` object where nodes are genomes and edges are ANI values.  This graph is converted into subgraphs of connected components whose edges are connected by a particular threshold such as 95% ANI [default] as recommended by the authors for species-level clustering.  These species-level clusters (SLC) are then partitioned and `MMSEQS2` is then run on each SLC panproteome to SLC-specific Protein Clusters (SSPC) previously referred to as Sample-specific Orthogroups (SSO). Pangenome tables are created for each SLC which includes the number of proteins for a protein-cluster are detected in each genome.  Each domain (i.e., prokaryotic, eukaryotic, viral) are clustered separately.

**Conda Environment**: `conda activate VEBA-cluster_env`

```
usage: cluster.py -m <mags> -a <proteins> -o <output_directory> -t 95

    Running: cluster.py v2023.6.14 via Python v3.9.15 | /expanse/projects/jcl110/anaconda3/envs/VEBA-preprocess_env/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -i INPUT, --input INPUT
                        path/to/input.tsv, Format: Must include the follow columns (No header) [organism_type]<tab>[id\_sample]<tab>[id_mag]<tab>[genome]<tab>[proteins] but can include additional columns to the right (e.g., [cds]<tab>[gene_models]).  Suggested input is from `compile_genomes_table.py` script.
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        path/to/project_directory [Default: veba_output/cluster]
  -e, --no_singletons   Exclude singletons
  --no_local_clustering
                        Only do global clustering

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit

FastANI arguments:
  -A ANI_THRESHOLD, --ani_threshold ANI_THRESHOLD
                        FastANI | Species-level cluster (SLC) ANI threshold (Range (0.0, 100.0]) [Default: 95.0]
  --genome_cluster_prefix GENOME_CLUSTER_PREFIX
                        Cluster prefix [Default: 'SLC-
  --genome_cluster_suffix GENOME_CLUSTER_SUFFIX
                        Cluster suffix [Default: '
  --genome_cluster_prefix_zfill GENOME_CLUSTER_PREFIX_ZFILL
                        Cluster prefix zfill. Use 7 to match identifiers from OrthoFinder.  Use 0 to add no zfill. [Default: 0]
  --fastani_options FASTANI_OPTIONS
                        FastANI | More options (e.g. --arg 1 ) [Default: '']

MMSEQS2 arguments:
  -a ALGORITHM, --algorithm ALGORITHM
                        MMSEQS2 | {easy-cluster, easy-linclust} [Default: easy-cluster]
  -t MINIMUM_IDENTITY_THRESHOLD, --minimum_identity_threshold MINIMUM_IDENTITY_THRESHOLD
                        MMSEQS2 | SLC-Specific Protein Cluster (SSPC, previously referred to as SSO) percent identity threshold (Range (0.0, 100.0]) [Default: 50.0]
  -c MINIMUM_COVERAGE_THRESHOLD, --minimum_coverage_threshold MINIMUM_COVERAGE_THRESHOLD
                        MMSEQS2 | SSPC coverage threshold (Range (0.0, 1.0]) [Default: 0.8]
  --protein_cluster_prefix PROTEIN_CLUSTER_PREFIX
                        Cluster prefix [Default: 'SSPC-
  --protein_cluster_suffix PROTEIN_CLUSTER_SUFFIX
                        Cluster suffix [Default: '
  --protein_cluster_prefix_zfill PROTEIN_CLUSTER_PREFIX_ZFILL
                        Cluster prefix zfill. Use 7 to match identifiers from OrthoFinder.  Use 0 to add no zfill. [Default: 0]
  --mmseqs2_options MMSEQS2_OPTIONS
                        MMSEQS2 | More options (e.g. --arg 1 ) [Default: '']

```
**Output:**

* global/feature\_compression\_ratios.tsv - Feature compression ratios for each domain.  **Also includes a summary of the number of genomes, genome clusters, proteins, and protein clusters.**
* global/genome\_clusters.tsv - Machine-readable table for genome clusters `[id_genome_cluster, number_of_components, number_of_samples_of_origin, components, samples_of_origin]`
* global/identifier\_mapping.genomes.tsv - Identifier mapping for genomes `[id_genome, organism_type, sample_of_origin, id_genome_cluster, number_of_proteins, number_of_singleton_protein_clusters, ratio_of_protein_cluster_are_singletons]`
* global/identifier\_mapping.proteins.tsv - Identifier mapping for proteins `[id_protein, organism_type, id_genome, sample_of_origin, id_genome_cluster, id_protein_cluster]`
* global/identifier\_mapping.scaffolds.tsv - Identifier mapping for contigs `[id_scaffold,	organism_type, id_genome, sample_of_origin, id_genome_cluster]`
* global/mags\_to\_slcs.tsv
* global/protein\_clusters.tsv - Machine-readable table for protein clusters `[id_protein_cluster, number_of_components, number_of_samples_of_origin, components, samples_of_origin]`
* global/proteins\_to\_orthogroups.tsv - Identifier mapping between proteins and protein clusters `[id_protein, id_protein-cluster]`
* global/representative\_sequences.faa - Protein sequences for cluster representatives.  Header follows the following format: `id_protein-cluster id_original_protein`
* global/scaffolds\_to\_mags.tsv - Identifier mapping between contigs and genomes `[id\_contig, id_genome]`
* global/scaffolds\_to\_slcs.tsv - Identifier mapping between contigs and genome clusters `[id\_contig, id_genome-cluster]`
* global/pangenome_tables/*.tsv.gz - Pangenome tables for each SLC with prevalence values
* global/serialization/*.dict.pkl - Python dictionaries for clusters
* global/serialization/*.networkx_graph.pkl - NetworkX graphs for clusters
* local/* - If `--no_local_clustering` is not selected then all of the files are generated for local clustering

<p align="right"><a href="#readme-top">^__^</a></p>


#### *annotate.py*
**Annotates translated gene calls against UniRef, Pfam, KOFAM, VFDB, MiBIG, AMRFinder, and AntiFam**

Annotation is performed using best hit annotations and profile HMMs. Proteins are aligned against `UniRef50/90`, `MiBIG`, and `VFDB` using `Diamond`. Protein domains are identified for `Pfam`, `NCBIfam-AMRFinder`,  and `AntiFam` using `HMMER3` and `KEGG` orthology using `KOFAMSCAN`.  Functionality for annotating all proteins or only protein cluster representatives and propogating annotations across cluster membership.


**Conda Environment**: `conda activate VEBA-annotate_env`

```
usage: annotate.py -a <proteins> -o <output_directory> |Optional: -i <identifier_mapping or -c <protein_clusters>]
Warnings:Proteins >100k will cause HMMSearch to crash.  Please prefilter using `seqkit seq -M 100000`

    Running: annotate.py v2023.6.20 via Python v3.9.15 | /expanse/projects/jcl110/anaconda3/envs/VEBA-preprocess_env/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -a PROTEINS, --proteins PROTEINS
                        Either path/to/proteins.faa or a directory of fasta files using [-x]
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        path/to/project_directory [Default: veba_output/annotation]
  -i IDENTIFIER_MAPPING, --identifier_mapping IDENTIFIER_MAPPING
                        Tab-seperated value table of [id_protein]<tab>[id\_contig]<tab>[id\_genome], No header.  Cannot be used with --protein_clusters
  -c PROTEIN_CLUSTERS, --protein_clusters PROTEIN_CLUSTERS
                        Tab-seperated value table of [id_protein]<tab>[id_protein_cluster].  Use this if the --proteins are representative sequences.  Cannot be used with --identifier_mapping
  -x EXTENSION, --extension EXTENSION
                        Fasta file extension for proteins if a directory is provided for --proteins [Default: faa]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  --keep_temporary_directory
                        Keep temporary directory [Default is to remove]
  -v, --version         show program's version number and exit

Database arguments:
  -u UNIREF, --uniref UNIREF
                        UniRef database to use {uniref90, uniref50}.  uniref90 receommended for well-characterized systems and uniref50 for less characterized systems [Default: uniref90]
  --veba_database VEBA_DATABASE
                        VEBA database location.  [Default: $VEBA_DATABASE environment variable]

Diamond arguments:
  --diamond_sensitivity DIAMOND_SENSITIVITY
                        Diamond | Sensitivity [Default:  '']
  --diamond_evalue DIAMOND_EVALUE
                        Diamond | E-Value [Default: 0.001]
  --diamond_options DIAMOND_OPTIONS
                        Diamond | More options (e.g. --arg 1 ) [Default: '']

HMMSearch arguments:
  --hmmsearch_options HMMSEARCH_OPTIONS
                        Diamond | More options (e.g. --arg 1 ) [Default: '']

KOFAMSCAN arguments:
  --kofamscan_options KOFAMSCAN_OPTIONS
                        Diamond | More options (e.g. --arg 1 ) [Default: '']

```

**Output:**

* annotations.tsv.gz - Concatenated annotations from Diamond (UniRef, MiBIG, VFDB), HMMSearch (Pfam, AntiFam), and KOFAMSCAN (KEGG)
* annotations.proteins.tsv.gz - Propogated annotations if clusters are provided
* module_completion_ratios.genomes.tsv - KEGG module completion ratios for each genomes [Only if --identifier_mapping is provided]
* module_completion_ratios.genome_clusters.tsv - KEGG module completion ratios for each genome clusters [Only if --identifier_mapping is provided]

<p align="right"><a href="#readme-top">^__^</a></p>

#### *phylogeny.py*
**Constructs phylogenetic trees given a marker set**

The phylogeny module is a tool used for phylogenetic inference and constructing phylogenetic trees for genomes given a reference marker set. This is performed by the following method: 1) identifying marker proteins using `HMMSearch` from the `HMMER3` suite; 2) creating protein alignments for each marker identified `MUSCLE`; 3) trimming the alignments using `ClipKIT`; 4) concatenating the alignments; 5) approximately-maximum-likelihood phylogenetic inference using `FastTree2` ; and 6) optional maximum likelihood phylogenetic inference using `IQ-TREE2`.  An option to include marker-specific noise cutoff scores is also available using the `--scores_cutoff` parameter.  Poor-quality genomes that do not meet a threshold in the proportion of markers in the reference are removed using the `--minimum_markers_aligned_ratio` parameter.  Similarly, non-informative markers that are not prevalent in the query genomes are removed using the `--minimum_genomes_aligned_ratio` parameter.

**Conda Environment**: `conda activate VEBA-phylogeny_env`


```
usage: phylogeny.py -d <database_hmms> -a <proteins> -o <output_directory>

    Running: phylogeny.py v2023.6.12 via Python v3.11.0 | /expanse/projects/jcl110/anaconda3/envs/VEBA-phylogeny_env/bin/python

options:
  -h, --help            show this help message and exit

Required I/O arguments:
  -d DATABASE_HMM, --database_hmm DATABASE_HMM
                        path/to/HMM database of markers
  -a PROTEINS, --proteins PROTEINS
                        Can be the following format: 1) Tab-seperated value table of [id_mag]<tab>[path/to/protein.fasta] (No header); 2) Files with list of filepaths [path/to/protein.fasta] (uses --extension); or 3) Directory of protein fasta  (uses --extension)
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        path/to/project_directory [Default: veba_output/phylogeny]
  -x EXTENSION, --extension EXTENSION
                        Fasta file extension for proteins if a list is provided [Default: faa]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit

HMMSearch arguments:
  --hmmsearch_threshold HMMSEARCH_THRESHOLD
                        HMMER | Threshold {cut_ga, cut_nc, gut_tc, e} [Default:  e]
  --hmmsearch_evalue HMMSEARCH_EVALUE
                        HMMER | E-Value [Default: 10.0]
  --hmmsearch_options HMMSEARCH_OPTIONS
                        HMMER | More options (e.g. --arg 1 ) [Default: '']
  -f HMM_MARKER_FIELD, --hmm_marker_field HMM_MARKER_FIELD
                        HMM reference type (accession, name) [Default: accession
  -s SCORES_CUTOFF, --scores_cutoff SCORES_CUTOFF
                        path/to/scores_cutoff.tsv. No header. [id_hmm]<tab>[score]

Alignment arguments:
  -A ALIGNMENT_ALGORITHM, --alignment_algorithm ALIGNMENT_ALGORITHM
                        Muscle alignment algorithm.  Align large input using Super5 algorithm if -align is too expensive. {align,super5} [Default: align]
  -g MINIMUM_GENOMES_ALIGNED_RATIO, --minimum_genomes_aligned_ratio MINIMUM_GENOMES_ALIGNED_RATIO
                        Minimum ratio of genomes include in alignment. This removes markers that are under represented. [Default: 0.95]
  -m MINIMUM_MARKERS_ALIGNED_RATIO, --minimum_markers_aligned_ratio MINIMUM_MARKERS_ALIGNED_RATIO
                        Minimum ratio of markers aligned. This removes genomes with few markers. Note, this is based on detected markers and NOT total markers in original HMM. [Default: 0.2]
  --muscle_options MUSCLE_OPTIONS
                        MUSCLE | More options (e.g. --arg 1 ) [Default: '']
  --clipkit_mode CLIPKIT_MODE
                        ClipKIT | Trimming mode [Default: smart-gap]
  --clipkit_options CLIPKIT_OPTIONS
                        ClipKIT | More options (e.g. --arg 1 ) [Default: '']

Tree arguments:
  --fasttree_options FASTTREE_OPTIONS
                        FastTree | More options (e.g. --arg 1 ) [Default: '']
  --no_iqtree           IQTree | Don't run IQTree
  --iqtree_model IQTREE_MODEL
                        IQTree | Model finder [Default: MFP]
  --iqtree_mset IQTREE_MSET
                        IQTree | Model set to choose from [Default: WAG,LG]
  --iqtree_bootstraps IQTREE_BOOTSTRAPS
                        IQTree | Bootstraps [Default: 1000]
  --iqtree_options IQTREE_OPTIONS
                        IQTree | More options (e.g. --arg 1 ) [Default: '']

```

**Output:**

* alignment_table.boolean.tsv.gz - Alignment table of (n = genomes, m = markers, ij=fasta pass qc)
* concatenated_alignment.fasta - Concatenated protein alignment of all marker hits
* concatenated_alignment.fasttree.nw - FastTree newick format based on concatenated alignment
* prefiltered_alignment_table.tsv.gz - Prefiltered alignment table of (n = genomes, m = markers, ij=fasta alignment)
* output.treefile - IQTREE2 newick format based on concatenated alignment (if --no_iqtree is not selected)

<p align="right"><a href="#readme-top">^__^</a></p>

#### *index.py*
**Builds local or global index for alignment to genomes**

The index module creates reference indices for alignments in both local or global paradigms. In the local paradigm, an index is created for all the assembled genomes concatenated together for each sample. This is useful in situations where perfectly paired metagenomics and metatranscriptomics are available where the metatranscriptomics can be mapped directly to the de novo reference generated from the metagenomics.  However, this is not applicable in all cases such as when there is not a perfect overlap between metagenomics and metatranscriptomics.  In this global paradigm, assembled genomes are concatenated across all samples and an alignment index is created for this concatenated reference.  Currently, `Bowtie2`  is the only alignment software packages supported. 

**Conda Environment**: `conda activate VEBA-mapping_env`

```
usage: index.py -i <mags> -o <output> --heatmap_output <pdf>

    Running: index.py v2023.5.8 via Python v3.11.0 | /expanse/projects/jcl110/anaconda3/envs/VEBA-phylogeny_env/bin/python

options:
  -h, --help            show this help message and exit

Required I/O arguments:
  -r REFERENCES, --references REFERENCES
                        local mode: [id\_sample]<tab>[path/to/reference.fa] and global mode: [path/to/reference.fa]
  -g GENE_MODELS, --gene_models GENE_MODELS
                        local mode: [id\_sample]<tab>[path/to/reference.gff] and global mode: [path/to/reference.gff]
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        path/to/project_directory [Default: veba_output/index]
  -m MINIMUM_CONTIG_LENGTH, --minimum_contig_length MINIMUM_CONTIG_LENGTH
                        Minimum contig length [Default: 1500]
  -M MODE, --mode MODE  Concatenate all references with global and build index or build index for each reference {global, local, infer}

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Use -1 for completely random. Use 0 for consecutive random states.  Use any other positive integer for the same random state for all iterations [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit

Bowtie2 Index arguments:
  --bowtie2_build_options BOWTIE2_BUILD_OPTIONS
                        bowtie2-build | More options (e.g. --arg 1 ) [Default: '']

```

**Output (global):**

* reference.fa.gz - Concatenated reference fasta
* reference.fa.gz.\*.bt2 - Bowtie2 index of reference fasta
* reference.gff - Concatenated gene models
* reference.saf - SAF format for reference

**Output (local):**

* [id\_sample]/reference.fa.gz - Concatenated reference fasta
* [id\_sample]/reference.fa.gz.\*.bt2 - Bowtie2 index of reference fasta
* [id\_sample]/reference.gff - Concatenated gene models
* [id\_sample]/reference.saf - SAF format for reference

<p align="right"><a href="#readme-top">^__^</a></p>

#### *mapping.py*
**Aligns reads to local or global index of genomes**

The mapping module uses local or global reference indices generated by the index module and aligns reads using `Bowtie2`.  The alignment files are sorted to produce sorted BAM files using `Samtools` which are then indexed. Coverage is calculated for contigs via `Samtools` and genome spatial coverage (i.e., ratio of bases covered in genome) is provided.  Reads from the sorted BAM files are then fed into featureCounts to produce gene-level counts, orthogroup-level counts, MAG-level counts, and SLC-level counts. 

**Conda Environment**: `conda activate VEBA-mapping_env`

```
usage: mapping.py -1 <reads_1.fq> -2 <reads_2.fq> -n <name> -o <output_directory> -x <reference_directory>

    Running: mapping.py v2023.5.15 via Python v3.11.0 | /expanse/projects/jcl110/anaconda3/envs/VEBA-phylogeny_env/bin/python

options:
  -h, --help            show this help message and exit

Required I/O arguments:
  -1 FORWARD_READS, --forward_reads FORWARD_READS
                        path/to/reads_1.fastq
  -2 REVERSE_READS, --reverse_reads REVERSE_READS
                        path/to/reads_2.fastq
  -n NAME, --name NAME  Name of sample
  -o PROJECT_DIRECTORY, --project_directory PROJECT_DIRECTORY
                        path/to/project_directory [Default: veba_output/mapping]

Reference arguments:
  -x REFERENCE_INDEX, --reference_index REFERENCE_INDEX
                        path/to/bowtie2_index. Either a file or directory. If directory, then it assumes the index is named `reference.fa.gz`
  -r REFERENCE_FASTA, --reference_fasta REFERENCE_FASTA
                        path/to/reference.fasta. If not provided then it is set to the --reference_index
  -a REFERENCE_GFF, --reference_gff REFERENCE_GFF
                        path/to/reference.gff. If not provided then --reference_index must be a directory that contains the file: 'reference.gff'
  -s REFERENCE_SAF, --reference_saf REFERENCE_SAF
                        path/to/reference.saf. If not provided then --reference_index must be a directory that contains the file: 'reference.saf'

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv. Must have at least 2 columns [name, executable] [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint
  -v, --version         show program's version number and exit

Bowtie2 arguments:
  --retain_unmapped_reads RETAIN_UNMAPPED_READS
                        Retain reads that do not map to reference. 0=No, 1=yes [Default: 1]
  --bowtie2_options BOWTIE2_OPTIONS
                        Bowtie2 | More options (e.g. --arg 1 ) [Default: '']
                        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

featureCounts arguments:
  -g ATTRIBUTE_TYPE, --attribute_type ATTRIBUTE_TYPE
                        Attribute type in GTF/GFF file. [Default: gene_id]
  -t FEATURE_TYPE, --feature_type FEATURE_TYPE
                        Feature type in GTF/GFF file. [Default: CDS]
  --retain_featurecounts RETAIN_FEATURECOUNTS
                        Retain feature counts output table (a slimmer version is output regardless). 0=No, 1=yes [Default: 0]
  --featurecounts_options FEATURECOUNTS_OPTIONS
                        featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/

Identifier arguments:
  --proteins_to_orthogroups PROTEINS_TO_ORTHOGROUPS
                        path/to/protein_to_orthogroup.tsv, [id\_orf]<tab>[id_orthogroup], No header
  --scaffolds_to_bins SCAFFOLDS_TO_BINS
                        path/to/scaffold_to_bins.tsv, [id_scaffold]<tab>[id_bin], No header
  --scaffolds_to_clusters SCAFFOLDS_TO_CLUSTERS
                        path/to/scaffold_to_cluster.tsv, [id_scaffold]<tab>[id_cluster], No header

```

**Output:**


* mapped.sorted.bam - Sorted BAM file
* mapped.sorted.bam.bai - Sorted BAM file index
* mapped.sorted.bam.coverage.tsv.gz - Samtools coverage table
* genome_spatial_coverage.tsv.gz - Spatial coverage for genome (i.e., ratio of bases covered) [Only if --scaffolds_to_bins is provided]
* counts.orfs.tsv.gz - ORF-level counts table
* counts.scaffolds.tsv.gz - Contig-level counts table
* counts.mags.tsv.gz - MAG-level counts table [Only if --scaffolds_to_bins is provided]
* counts.clusters.tsv.gz - SLC-level counts table [Only if --scaffolds_to_clusters is provided]
* counts.orthogroups.tsv.gz - Orthogroup-level counts table [Only if --orf_to_orthogroups is provided]
* unmapped_1.fastq.gz - Unmapped reads (forward)
* unmapped_2.fastq.gz - Unmapped reads (reverse)

<p align="right"><a href="#readme-top">^__^</a></p>

#### *biosynthetic.py*
**Identify biosynthetic gene clusters in prokaryotes and fungi**

The biosynthetic module is a wrapper around `antiSMASH`.  It produces a tabular output that is machine-readale and easier to parse than the GBK and JSON files produced by antiSMASH.  Novelty scores are calculated using the ratio of proteins that align to `MiBIG`. 

```

usage: biosynthetic.py -i <genomes_gene-models.tsv> -o <output_directory> -t bacteria | Suggested input is from `compile_genomes_table.py` script. Use cut -f3,4,7

    Running: biosynthetic.py v2023.7.10 via Python v3.9.9 | /Users/jespinoz/anaconda3/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -i INPUT, --input INPUT
                        path/to/input.tsv, Format: Must include the follow columns (No header) [id_mag]<tab>[genome]<tab>[gene_models].  Suggested input is from `compile_genomes_table.py` script. Use cut -f3,4,7
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        path/to/project_directory [Default: veba_output/biosynthetic]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit

Database arguments:
  --veba_database VEBA_DATABASE
                        VEBA database location.  [Default: $VEBA_DATABASE environment variable]

antiSMASH arguments:
  -t TAXON, --taxon TAXON
                        Taxonomic classification of input sequence {bacteria,fungi} [Default: bacteria]
  --minimum_contig_length MINIMUM_CONTIG_LENGTH
                        Minimum contig length.  [Default: 1500]
  -d ANTISMASH_DATABASE, --antismash_database ANTISMASH_DATABASE
                        antiSMASH | Database directory path [Default: /Users/jespinoz/anaconda3/lib/python3.9/site-packages/antismash/databases]
  -s HMMDETECTION_STRICTNESS, --hmmdetection_strictness HMMDETECTION_STRICTNESS
                        antiSMASH | Defines which level of strictness to use for HMM-based cluster detection {strict,relaxed,loose}  [Default: relaxed]
  --tta_threshold TTA_THRESHOLD
                        antiSMASH | Lowest GC content to annotate TTA codons at [Default: 0.65]
  --antismash_options ANTISMASH_OPTIONS
                        antiSMASH | More options (e.g. --arg 1 ) [Default: '']

Diamond arguments:
  --diamond_sensitivity DIAMOND_SENSITIVITY
                        Diamond | Sensitivity [Default:  '']
  --diamond_evalue DIAMOND_EVALUE
                        Diamond | E-Value [Default: 0.001]
  --diamond_options DIAMOND_OPTIONS
                        Diamond | More options (e.g. --arg 1 ) [Default: '']

Novelty score threshold arguments:
  --pident PIDENT       pident lower bound [float:0 ≤ x < 100] [Default: 0]
  --qcovhsp QCOVHSP     qcovhsp lower bound [float:0 ≤ x < 100] [Default: 0]
  --scovhsp SCOVHSP     scovhsp lower bound [float:0 ≤ x < 100] [Default: 0]
  --evalue EVALUE       e-value lower bound [float:0 < x < 1] [Default: 1e-3]

```

**Output:**

* bgc\_clusters.tsv - BGC to BGC nucleotide cluster
* bgc\_protocluster-types.tsv.gz - Summary of BGCs detected organized by type.  Also includes summary of BGCs that are NOT on contig edge.
* bgcs.representative\_sequences.fasta.gz - Full length BGC nucleotide cluster representatives
* component\_clusters.tsv - BGC protein to BGC protein cluster
* components.representative\_sequences.faa.gz - BGC protein cluster representatives
* fasta/[id/_genome].faa/fasta.gz - BGC sequences in protein and nucleotide space
* genbanks/[id\_genome]/*.gbk - Genbank formatted antiSMASH results
* homology.tsv.gz - Diamond results for MIBiG and VFDB
* identifier\_mapping.bgcs.tsv.gz - All of the BGCs in tabular format organized by genome, contig, region, and gene.
* identifier\_mapping.components.tsv.gz - All of the BGC components (i.e., genes in BGC) in tabular format organized by genome, contig, region, and gene.
* krona.html - HTML showing Krona plot for number of BGCs per protocluster-type.
* krona.tsv - Data to produce Krona plot
* prevalence\_tables/bgcs.tsv.gz - Genome vs. BGC nucleotide cluster prevalence table
* prevalence\_tables/components.tsv.gz - Genome vs. BGC protein cluster prevalence table


<p align="right"><a href="#readme-top">^__^</a></p>


#### *profile-taxonomy.py*
**Taxonomic profiling of *de novo* genomes**

The profile-taxonomy module does the following: 0) builds a Sylph sketch database (63) for non-viral and viral genomes using the compile_custom_sylph_sketch_database_from_genomes.py script prior to running the module;  1) converts paired reads to a query sketch database using Sylph; 2) profiles the genome sketch databases using the query sketch database generated from the reads; 3) reformats the Sylph output tables; and 4) aggregates abundances with respect to SLC if clusters are provided. 

```
usage: profile-taxonomy.py -1 <forward_reads.fq> -2 <reverse_reads.fq>|-s <sketch> -n <name> -o <output_directory> -d <db_1.syldb db_2.syldb ... db_n.syldb>

    Running: profile-taxonomy.py v2023.12.19 via Python v3.10.12 | /expanse/projects/jcl110/miniconda3/envs/VEBA-profile_env/bin/python

options:
  -h, --help            show this help message and exit

Required I/O arguments:
  -1 FORWARD_READS, --forward_reads FORWARD_READS
                        path/to/forward_reads.fq[.gz]
  -2 REVERSE_READS, --reverse_reads REVERSE_READS
                        path/to/reverse_reads.fq[.gz]]
  -s READS_SKETCH, --reads_sketch READS_SKETCH
                        path/to/reads_sketch.sylsp (e.g., sylph sketch output) (Cannot be used with --forward_reads and --reverse_reads)
  -n NAME, --name NAME  Name of sample
  -d SYLPH_DATABASES [SYLPH_DATABASES ...], --sylph_databases SYLPH_DATABASES [SYLPH_DATABASES ...]
                        Sylph database(s) with all genomes.  Can be multiple databases delimited by spaces.  Use compile_custom_sylph_sketch_database_from_genomes.py to build database.
  -o PROJECT_DIRECTORY, --project_directory PROJECT_DIRECTORY
                        path/to/project_directory [Default: veba_output/profiling/taxonomy]
  -c GENOME_CLUSTERS, --genome_clusters GENOME_CLUSTERS
                        path/to/mags_to_slcs.tsv. [id_genome]<tab>[id_genome-cluster], No header. Aggregates counts for genome clusters.
  -F {sketch,paired}, --input_reads_format {sketch,paired}
                        Input reads format {paired, sketch} [Default: auto]
  -x EXTENSION, --extension EXTENSION
                        Fasta file extension for bins. Assumes all genomes have the same file extension. [Default: fa]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit
  --tmpdir TMPDIR       Set temporary directory

Sylph sketch arguments (Fastq):
  --sylph_sketch_k {21,31}
                        Sylph sketch [Fastq] |  Value of k. Only k = 21, 31 are currently supported. [Default: 31]
  --sylph_sketch_minimum_spacing SYLPH_SKETCH_MINIMUM_SPACING
                        Sylph sketch [Fastq] |  Minimum spacing between selected k-mers on the genomes [Default: 30]
  --sylph_sketch_subsampling_rate SYLPH_SKETCH_SUBSAMPLING_RATE
                        Sylph sketch [Fastq] |  Subsampling rate.	 sylph runs without issues if the -c for all genomes is ≥ the -c for reads.  [Default: 100]
  --sylph_sketch_options SYLPH_SKETCH_OPTIONS
                        Sylph sketch [Fastq] | More options for `sylph sketch` (e.g. --arg 1 ) [Default: '']

Sylph profile arguments:
  --sylph_profile_minimum_ani SYLPH_PROFILE_MINIMUM_ANI
                        Sylph profile | Minimum adjusted ANI to consider (0-100). [Default: 95]
  --sylph_profile_minimum_number_kmers SYLPH_PROFILE_MINIMUM_NUMBER_KMERS
                        Sylph profile | Exclude genomes with less than this number of sampled k-mers.  Default is 50 in Sylph but lowering to 20 accounts for viruses and small CPR genomes. [Default: 20]
  --sylph_profile_minimum_count_correct SYLPH_PROFILE_MINIMUM_COUNT_CORRECT
                        Sylph profile | Minimum k-mer multiplicity needed for coverage correction. Higher values gives more precision but lower sensitivity [Default: 3]
  --sylph_profile_options SYLPH_PROFILE_OPTIONS
                        Sylph profile | More options for `sylph profile` (e.g. --arg 1 ) [Default: '']
  --header              Include header in taxonomic abundance tables

```

**Output:**

* reads.sylsp - Reads sketch if paired-end reads were provided
* sylph\_profile.tsv.gz - Output of `sylph profile`
* taxonomic_abundance.tsv.gz - Genome-level taxonomic abundance (No header)
* taxonomic_abundance.clusters.tsv.gz - SLC-level taxonomic abundance (No header, if --genome_clusters wer provided)

#### *profile-pathway.py*
**Pathway profiling of *de novo* genomes**
The profile-pathway module does the following: 0) builds a custom HUMAnN database based on protein annotations, identifier mapping tables, and taxonomy assignments using the compile_custom_humann_database_from_annotations.py script prior to running the module; 1) either accepts pre-joined reads, joins paired end reads using bbmerge.sh from BBSuite, or a BAM file of paired-end reads and joins them; 2) builds a Diamond database of proteins from the custom HUMAnN annotation table; 3) uses HUMAnN for pathway profiling of the joined reads using the custom HUMAnN database (16); and 4) reformats the output files.

```
usage: profile-pathway.py -1 <forward_reads.fq> -2 <reverse_reads.fq> -n <name> -o <output_directory>

    Running: profile-pathway.py v2023.11.30 via Python v3.10.12 | /expanse/projects/jcl110/miniconda3/envs/VEBA-profile_env/bin/python

options:
  -h, --help            show this help message and exit

Required reads arguments:
  -1 FORWARD_READS, --forward_reads FORWARD_READS
                        path/to/forward_reads.fq (Requires --reverse_reads, cannot be used with --joined_reads or --bam)
  -2 REVERSE_READS, --reverse_reads REVERSE_READS
                        path/to/reverse_reads.fq (Requires --forward_reads, cannot be used with --joined_reads or --bam)
  -j JOINED_READS, --joined_reads JOINED_READS
                        path/to/joined_reads.fq (e.g., bbmerge.sh output) (Cannot be used with --forward_reads, --reverse_reads, or --bam)
  -b BAM, --bam BAM     path/to/mapped.sorted.bam file aligned to genomes (Cannot be used with --forward_reads, --reverse_reads, or --joined_reads)
  -F INPUT_READS_FORMAT, --input_reads_format INPUT_READS_FORMAT
                        Input reads format {paired, joined, bam} [Default: auto]

Required database arguments:
  -i IDENTIFIER_MAPPING, --identifier_mapping IDENTIFIER_MAPPING
                        Identifier mapping which includes [id_protein]<tab>[id_uniref]<tab>[length]<tab>[lineage].  In VEBA, you can use `compile_custom_humann_database_from_annotations.py`.
                        https://github.com/biobakery/humann#custom-reference-database-annotations
  -f FASTA, --fasta FASTA
                        Protein fasta to build database
  -d DIAMOND_DATABASE, --diamond_database DIAMOND_DATABASE
                        Diamond database with all proteins from --identifier_mapping

Required I/O arguments:
  -n NAME, --name NAME  Name of sample
  -o PROJECT_DIRECTORY, --project_directory PROJECT_DIRECTORY
                        path/to/project_directory [Default: veba_output/profiling/pathways]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit
  --tmpdir TMPDIR       Set temporary directory

bbmerge.sh arguments:
  --minimum_merge_overlap MINIMUM_MERGE_OVERLAP
                        bbmerge.sh | Minimum number of overlapping bases to allow merging. [Default: 12]
  --bbmerge_options BBMERGE_OPTIONS
                        bbmerge.sh options (e.g. --arg 1 ) [Default: '']

HUMAnN arguments:
  --search_mode SEARCH_MODE
                        HUMAnN | Search for uniref50 or uniref90 gene families {uniref50, uniref90, auto} [Default: 'auto']
  --pathways PATHWAYS   HUMAnN | The database to use for pathway computations {metacyc, unipathway} [Default: 'metacyc']
  -e EVALUE, --evalue EVALUE
                        HUMAnN | The evalue threshold to use with the translated search [Default: 1.0]
  -m TRANSLATED_IDENTITY_THRESHOLD, --translated_identity_threshold TRANSLATED_IDENTITY_THRESHOLD
                        HUMAnN | Identity threshold for translated alignments [Default: Tuned automatically (based on uniref mode) unless a custom value is specified]
  -q TRANSLATED_QUERY_COVERAGE_THRESHOLD, --translated_query_coverage_threshold TRANSLATED_QUERY_COVERAGE_THRESHOLD
                        HUMAnN | Query coverage threshold for translated alignments [Default: 90.0]
  -s TRANSLATED_SUBJECT_COVERAGE_THRESHOLD, --translated_subject_coverage_threshold TRANSLATED_SUBJECT_COVERAGE_THRESHOLD
                        HUMAnN | Subject coverage threshold for translated alignments [Default: 50.0]
  --humann_memory HUMANN_MEMORY
                        HUMAnN | Memory use mode {minimum, maximum} [Default: 'minimum']
  --humann_options HUMANN_OPTIONS
                        HUMAnN options (e.g. --arg 1 ) [Default: '']
  ```

**Output:**

* reads.seqkit_stats.tsv - Sequence stats for input reads
* humann\_pathabundance.tsv - Stratified abundance of taxa-specific metabolic pathways
* humann\_pathcoverage.tsv - Stratified pathway completion ratio of taxa-specific metabolic pathways
* humann\_genefamilies.tsv - Stratified abundance of taxa-specific gene families
* humann\_diamond\_unaligned.fa.gz - Joined reads that did not align to database
* humann\_diamond\_aligned.tsv.gz - Aligned reads from translated blast search to database (blast6 format)
___________________________________________________________________
## Developmental


#### *amplicon.py*
**Automated read trim position detection, DADA2 ASV detection, taxonomic classification, and file conversion** 

The amplicon module is a wrapper around `QIIME2`'s implementation of the `DADA2` ASV pipeline which has been fairly standardized.  This works exclusively on paired-end short reads and is not designed for single-end reads nor long reads (the latter may be adapted later).  The experimental portion of this module is the automatic detection of forward and reverse trim.  This module first imports reads into a `QIIME2 Artifact` object, summarizes reads, and gets position-specific fastq statistics.  The amplicon module uses the position-specific fastq statistics to suggest forward and reverse trim positions (this part is experimental, please use `--inspect_trim_regions` to manually check the quality plots to ensure it is where you would cut). Next ASVs are detected via `DADA2` and denoising statistics are calculated.  After ASVs are detected, taxonomy is classified using classification modules provided by user (e.g., [silva-138-99-nb-classifier.qza](https://data.qiime2.org/2022.8/common/silva-138-99-nb-classifier.qza) followed by phylogenetic inference.  Finally, `QIIME2` and `BIOM` formatted files are converted into tab-separated value tables and fasta files.

```
usage: amplicon.py -i <reads_table.tsv> -c <classifier.qza> -o <output_directory>

    Running: amplicon.py v2022.10.24 via Python v3.9.7 | /Users/jespinoz/anaconda3/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -i READS_TABLE, --reads_table READS_TABLE
                        path/to/reads_table.tsv. 3 columns separated by tabs with the following header: [sample-id <tab> forward-absolute-filepath <tab> reverse-absolute-filepath]
                        A utility script is provided: compile_reads_table.py
  -c CLASSIFIER, --classifier CLASSIFIER
                        path/to/feature_classifier. Data Resources: https://docs.qiime2.org/2022.8/data-resources/
  -o PROJECT_DIRECTORY, --project_directory PROJECT_DIRECTORY
                        path/to/project_directory [Default: veba_output/amplicon]

Utility arguments:
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit
  --tmpdir TMPDIR       Set temporary directory

Trim detection arguments:
  --inspect_trim_regions
                        Manually inspect trim regions then rerun [PLEASE USE THIS TO CHECK TRIMMING SUGGESTIONS AS THEY ARE CURRENTLY EXPERIMENTAL]
  -f FORWARD_TRIM, --forward_trim FORWARD_TRIM
                        Specify forward trim position
  -r REVERSE_TRIM, --reverse_trim REVERSE_TRIM
                        Specify reverse trim position
  -q MINIMUM_QUALITY, --minimum_quality MINIMUM_QUALITY
                        Minimum quality value
  -m MINIMUM_LENGTH, --minimum_length MINIMUM_LENGTH
                        Minimum length.  If minimum quality value makes length shorter than this then an error will yield with which samples are responsible [Default: 100]
  -w WINDOW_SIZE, --window_size WINDOW_SIZE
                        Window size [Default: 4]
  -l MAXIMUM_AVERAGE_LOSS, --maximum_average_loss MAXIMUM_AVERAGE_LOSS
                        Maximum average loss for window size [Default: --window_size]

DADA2 arguments:
  --minimum_overlap MINIMUM_OVERLAP
                        DADA2 | The minimum length of the overlap required for merging the forward and reverse reads. [Default: 12]
  --dada2_options DADA2_OPTIONS
                        Additional DADA2 options. '--arg value'

Phylogeny arguments:
  --phylogeny_mode PHYLOGENY_MODE
                        QIIME2 phylogeny submodule [Default: align-to-tree-mafft-fasttree]
  --phylogeny_options PHYLOGENY_OPTIONS
                        Additional options. '--arg value'

Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)
```

**Output:**

* aligned-dna-sequences.fasta - Aligned ASV reference sequences
* dna-sequences.fasta - ASV reference sequences
* dna-sequences.with_taxonomy.fasta - ASV reference sequences with taxonomy/confidence info in description
* feature-table.biom - BIOM format ASV feature table
* feature-table.tsv - Tab-seperated values ASV table (rows=ASV, columns=samples, skiprows=1)
* stats.tsv - Read statistics
* taxonomy.tsv - Taxonomy table [ASV]<tab>[Lineage]<tab>[Confidence]
* tree.nwk - Newick formatted tree

<p align="right"><a href="#readme-top">^__^</a></p>

___________________________________________________________________
## Deprecated

#### *assembly-sequential.py*
**Assemble metagenomes sequentially**

This method first uses biosyntheticSPAdes followed by either 1) [default] metaSPAdes; or 2) metaplasmidSPAdes and metaSPAdes.  The reads that are not mapped to the scaffolds in step N are used for assembly in step N+1. The contigs/scaffolds are concatenated for the finally assembly. The purpose of this module is to properly handle biosynthetic gene clusters in addition to metagenomes.

```
usage: assembly-sequential.py -1 <forward_reads.fq> -2 <reverse_reads.fq> -n <name> -o <output_directory>

    Running: assembly-sequential.py v2022.11.14 via Python v3.9.7 | /Users/jespinoz/anaconda3/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -1 FORWARD_READS, --forward_reads FORWARD_READS
                        path/to/forward_reads.fq
  -2 REVERSE_READS, --reverse_reads REVERSE_READS
                        path/to/reverse_reads.fq
  -n NAME, --name NAME  Name of sample
  -o PROJECT_DIRECTORY, --project_directory PROJECT_DIRECTORY
                        path/to/project_directory [Default: veba_output/assembly_sequential]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit
  --tmpdir TMPDIR       Set temporary directory
  -S, --remove_intermediate_scaffolds
                        Remove intermediate scaffolds.fasta.*. If this option is chosen, output files are not validated [Default is to keep]
  -B, --remove_intermediate_bam
                        Remove intermediate mapped.sorted.bam.*. If this option is chosen, output files are not validated [Default is to keep]

SPAdes arguments:
  --run_metaplasmidspades
                        SPAdes | Run metaplasmidSPAdes.  This may sacrifice MAG completeness for plasmid completeness.
  -m MEMORY, --memory MEMORY
                        SPAdes | RAM limit in Gb (terminates if exceeded). [Default: 250]
  --spades_options SPADES_OPTIONS
                        SPAdes | More options (e.g. --arg 1 ) [Default: '']
                        http://cab.spbu.ru/files/release3.11.1/manual.html

Bowtie2 arguments:
  --bowtie2_index_options BOWTIE2_INDEX_OPTIONS
                        bowtie2-build | More options (e.g. --arg 1 ) [Default: '']
  --bowtie2_options BOWTIE2_OPTIONS
                        bowtie2 | More options (e.g. --arg 1 ) [Default: '']

featureCounts arguments:
  --featurecounts_options FEATURECOUNTS_OPTIONS
                        featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/

Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)
```

**Output:**

* featurecounts.tsv.gz - featureCounts output for contig-level counts
* mapped.sorted.bam - Sorted BAM
* mapped.sorted.bam.bai - Sorted BAM index
* scaffolds.fasta - Assembly scaffolds (preferred over contigs by SPAdes documentation)
* scaffolds.fasta.\*.bt2 - Bowtie2 index of scaffolds
* scaffolds.fasta.saf - SAF formatted file for contig-level counts with featureCounts
* seqkit_stats.tsv.gz - Assembly statistics

<p align="right"><a href="#readme-top">^__^</a></p>

