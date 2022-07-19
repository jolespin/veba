```
 _    _ _______ ______  _______
  \  /  |______ |_____] |_____|
   \/   |______ |_____] |     |
                               
                                                                                                                                                                                                                  
```
### Description
The *Viral Eukaryotic Bacterial Archaeal* (VEBA) is an open-source software suite developed with all domains of microorganisms as the primary objective (not *post hoc* adjustments) including prokaryotic, eukaryotic, and viral organisms.  VEBA is an end-to-end metagenomics analysis suite that can directly address eukaryotic and viral genomes in addition to prokaryotic genomes with support for CPR. VEBA implements a novel iterative binning procedure and hybrid sample-specific/multisample framework.  To optimize the microeukaryotic gene calling and taxonomic classifications, VEBA includes a consensus microeukaryotic database containing protists and fungi from several existing databases. VEBA also provides a unique clustering-based dereplication strategy allowing for sample-specific genomes and genes to be directly compared across non-overlapping biological samples.  In addition, VEBA automates the detection of candidate phyla radiation bacteria and implements the appropriate genome quality assessments for said organisms.  
___________________________________________________________________
### Citation
Espinoza et al. (In review)

___________________________________________________________________

### Development
VEBA is currently under active development. If you are interested in requesting features or wish to report a bug, please post a GitHub issue prefixed with the tag [Feature Request] and [Bug], respectively.
___________________________________________________________________

### Installation and databases
Please refer to the [*Installation and Database Configuration Guide*](install/README.md) for software installation and database configuration.

___________________________________________________________________
### Getting started with *VEBA*

Please refer to the *Walkthrough Guides*](walkthroughs/README.md) for tutorials and workflows on how to get started.

___________________________________________________________________


### Output structure
*VEBA*'s is built on the [GenoPype](https://github.com/jolespin/genopype) archituecture which creates a reproducible and easy-to-navigate directory structure.  *GenoPype*'s philosophy is to use the same names for all files but to have sample names as subdirectories.  This makes it easier to glob files for grepping, concatenating, etc.  *VEBA* 
e.g., 

Here is an example of *GenoPype*'s layout:

```
# Project directory
project_directory/

# Temporary directory
project_directory/tmp/

# Log directory
project_directory/logs/
project_directory/logs/[step]__[program-name].e
project_directory/logs/[step]__[program-name].o
project_directory/logs/[step]__[program-name].returncode

# Checkpoint directory
project_directory/checkpoints/
project_directory/checkpoints/

# Intermediate directories for each step
project_directory/intermediate/
project_directory/intermediate/[step]__[program-name]/

# Output directory
project_directory/output/

# Commands
project_directory/commands.sh
```

For *VEBA*, it has all the directories created by `GenoPype` above but is built for having multiple samples under the same project:

```
ID="sample_1"

# Main output directory
veba_output/

# Assembly directory
veba_output/assembly

# Assembly output for ${ID} sample
veba_output/assembly/${ID}/output/

# Prokaryotic binning for ${ID} sample
veba_output/binning/prokaryotic/${ID}/output/ 

# Eukaryotic binning
veba_output/binning/eukaryotic/${ID}/output/

# Viral binning
veba_output/binning/viral/${ID}/output/

```
* The above are default output locations but they can be customized.
___________________________________________________________________

### Modules
[![Schematic](images/Schematic.png)](images/Schematic.pdf)

* **preprocess** – Fastq quality trimming, adapter removal, decontamination, and read statistics calculations

* **assembly** – Assemble reads, align reads to assembly, and count mapped reads

* **coverage** – Align reads to (concatenated) reference and counts mapped reads

* **binning-prokaryotic** – Iterative consensus binning for recovering prokaryotic genomes with lineage-specific quality assessment

* **binning-eukaryotic** – Binning for recovering eukaryotic genomes with exon-aware gene modeling and lineage-specific quality assessment

* **binning-viral** – Detection of viral genomes and quality assessment

* **classify-prokaryotic** – Taxonomic classification and candidate phyla radiation adjusted quality 

* **classify-eukaryotic** – Taxonomic classification of eukaryotic genomes

* **classify-viral** – Taxonomic classification and isolation source of viral genomes

* **cluster** – Species-level clustering of genomes and lineage-specific orthogroup detection

* **annotate** – Annotates translated gene calls against NR, Pfam, and KOFAM

* **phylogeny** – Constructs phylogenetic trees given a marker set

* **index** – Builds local or global index for alignment to genomes
 
* **mapping** – Aligns reads to local or global index of genomes

#### preprocess – Fastq quality trimming, adapter removal, decontamination, and read statistics calculations
The preprocess module is a wrapper around our [fastq_preprocessor](https://github.com/jolespin/fastq_preprocessor) which is a modernized reimplementation of [KneadData](https://github.com/biobakery/kneaddata) that relies on fastp for ultra-fast automated adapter removal and quality trimming. Pairing of the trimmed reads is assessed and corrected using [BBTools’ reformat.sh](https://sourceforge.net/projects/bbmap). If the user provides a contamination database (e.g., the human reference genome), then trimmed reads are aligned using Bowtie2 (Langmead & Salzberg, 2012) and reads that do not map to the contamination database are stored. If the --retain_contaminated_reads flag is used then the contaminated reads are stored as well. Similarly, if a k-mer reference database is provided (e.g., ribosomal k-mers) then the trimmed or decontaminated reads are aligned against the reference database using BBTools’ bbduk.sh with an option for storing. By default, the none of the contaminated or k-mer analyzed reads are stored but regardless of the choice for retaining reads, the read sets are quantified using seqkit (Shen, Le, Li, & Hu, 2016) for accounting purposes (e.g., % contamination or % ribosomal). All sequences included were downloaded using Kingfisher (https://github.com/wwood/kingfisher-download), included in the preprocess environment, which a fast and flexible program for the procurement of sequencing files and their annotations from public data sources including ENA, NCBI SRA, Amazon AWS, and Google Cloud.

**Conda Environment:** `conda activate VEBA-preprocess_env`

```
usage: preprocess.py -1 <reads_1.fq> -2 <reads_2.fq> -n <name> -o <output_directory> |Optional| -x <reference_index> -k <kmer_database>

    Wrapper around github.com/jolespin/fastq_preprocessor
    Running: preprocess.py v2022.01.19 via Python v3.7.11 | /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-preprocess_env/bin/python

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
                        (e.g., Human GRCh38 from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids//GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz)
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
                        (e.g., ribokmers.fa.gz from https://drive.google.com/file/d/0B3llHR93L14wS2NqRXpXakhFaEk/view?usp=sharing)
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

#### assembly – Assemble reads, align reads to assembly, and count mapped reads
The assembly module optimizes the output for typical metagenomics workflows. In particular, the module does the following: 1) assembles reads using either metaSPAdes (Nurk, Meleshko, Korobeynikov, & Pevzner, 2017), SPAdes (A et al., 2012), rnaSPAdes (Bushmanova, Antipov, Lapidus, & Prjibelski, 2019), or any of the other task-specific assemblers installed with the SPAdes package (Antipov, Raiko, Lapidus, & Pevzner, 2020; Meleshko, Hajirasouliha, & Korobeynikov, 2021); 2) builds a Bowtie2 index for the scaffolds.fasta (or transcripts.fasta if rnaSPAdes is used); 3) aligns the reads using Bowtie2 to the assembly; 4) pipes the alignment file into Samtools (Li et al., 2009) to produce a sorted BAM file (necessary for many coverage applications); 5) counts the reads mapping to each scaffold via featureCounts (Liao, Smyth, & Shi, 2014); and 6) seqkit for useful assembly statistics such as N50, number of scaffolds, and total assembly size. This module automates many critical yet overlooked workflows dealing with assemblies. 

**Conda Environment**: `conda activate VEBA-assembly_env`

```
usage: assembly.py -1 <forward_reads.fq> -2 <reverse_reads.fq> -n <name> -o <output_directory>

    Running: assembly.py v2022.03.25 via Python v3.8.5 | /Users/jespinoz/anaconda3/bin/python

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

SPAdes arguments:
  -P SPADES_PROGRAM, --spades_program SPADES_PROGRAM
                        SPAdes | More options (e.g. --arg 1 ) [Default: 'metaspades.py']
                        http://cab.spbu.ru/files/release3.11.1/manual.html
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

```

**Output:**

* featurecounts.tsv.gz - featureCounts output for contig-level counts
* mapped.sorted.bam - Sorted BAM
* mapped.sorted.bam.bai - Sorted BAM index
* scaffolds.fasta - Assembly scaffolds (preferred over contigs by SPAdes documentation)
* scaffolds.fasta.\*.bt2 - Bowtie2 index of scaffolds
* scaffolds.fasta.saf - SAF formatted file for contig-level counts with featureCounts
* seqkit_stats.tsv.gz - Assembly statistics


#### coverage – Align reads to (concatenated) reference and counts mapped reads
The coverage module further optimizes the output for typical metagenomics workflows. In particular, the module does the following: 1) filters contigs based on a size filter (default 1500 bp); 2) builds a Bowtie2 index for the coassembly.fasta; 3) aligns the reads from all provided samples using Bowtie2 to the assembly; 4) pipes the alignment file into Samtools to produce a sorted BAM file; 5) counts the reads mapping to each scaffold via featureCounts; and 6) seqkit for useful assembly statistics such as N50, number of scaffolds, and total assembly size (Shen et al., 2016). The preferred usage for this module is after prokaryotic, eukaryotic, and viral binning has been performed and the unbinned contigs are merged into a single coassembly used as input.  The outputs of this module are expected to be used as a final pass through prokaryotic and eukaryotic binning modules.  

**Conda Environment**: `conda activate VEBA-assembly_env`


```
usage: coverage.py -f <reference.fasta> -r <reads.tsv> -o <output_directory>

    Running: coverage.py v2022.04.12 via Python v3.8.5 | /Users/jespinoz/anaconda3/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -f FASTA, --fasta FASTA
                        path/to/reference.fasta [Required]
  -r READS, --reads READS
                        path/to/reads_table.tsv with the following format: [id_sample]<tab>[path/to/r1.fastq.gz]<tab>[path/to/r2.fastq.gz], No header
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

#### binning-prokaryotic – Iterative consensus binning for recovering prokaryotic genomes with lineage-specific quality assessment
The prokaryotic binning module implements a novel iterative consensus binning procedure that uses CoverM (https://github.com/wwood/CoverM) for fast coverage calculations, multiple binning algorithms (MaxBin2 (marker set = 107); MaxBin2 (marker set = 40) (Wu et al., 2016); MetaBat2 (Kang et al., 2019); and CONCOCT (Alneberg et al., 2014), consensus dereplication and aggregate binning with DAS Tool (Sieber et al., 2018), the consensus domain wrapper for Tiara (Karlicki et al., 2022) for removing eukaryotes at the MAG level, and CheckM   for quality assessment where poor quality MAGs are removed (e.g., completeness < 50% and/or contamination ≥ 10). The novelty of this procedure is that the unbinned contigs are stored and fed back into the input of the binning procedure using a separate random seed state allowing for an exhaustive, yet effective, approach in extracting high quality and difficult to bin genomes; number of iterations specified by --n\_iter option. Gene calls are performed using Prodigal (Hyatt et al., 2010) and the gene models (GFF3 Format) are modified to include gene and contig identifiers for use with downstream feature counting software. Although CheckM can handle CPR it cannot do so with the typical, and recommended, lineage_wf directly in the current version but instead with a separate workflow. The prokaryotic binning module allows for basal bacteria to filter through intermediate genome quality checks, runs GTDB-Tk (Chaumeil et al., 2020) for genome classification, reruns CheckM CPR workflow for said genomes, and then updates the genome set with adjusted completeness and contamination scores.  The input alignment file is utilized using featureCounts to produce counts tables for the gene models and MAGs.  Lastly, genome statistics such as N50, number of scaffolds, and genome size are calculated using seqkit. Utility scripts, installed with VEBA, are run in the backend to modify prodigal gene models, consensus domain classification of MAGs using Tiara contig predictions, along with several fasta and pre/post-processing scripts. The input to this module is a fasta file (typically the scaffolds.fasta from metaSPAdes) and sorted BAM while the output includes the prokaryotic MAGs via Prodigal, gene models, identifier mappings, counts tables, CheckM output, GTDB-Tk output, and unbinned fasta.  MAG naming scheme for prokaryotes follows [SAMPLE]\_\_\[ALGORITHM\]\_\_P.\[ITERATION\]\_\_\[NAME] (e.g., SRR17458623\_\_METABAT2\_\_P.1\_\_bin.1)

**⚠️Note:** *If you have a lot of samples and a lot of contigs then use the --skip_maxbin2 flag because it takes MUCH longer to run.  For the Plastisphere it was going to take 40 hours per MaxBin2 run (there are 2 MaxBin2 runs) per iteration.  Metabat2 and CONCOCT can do the heavy lifting much faster and often with better results so it's recommended to skip MaxBin2 for larger datasets.*

**Conda Environment**: `conda activate VEBA-binning-prokaryotic_env`


```
usage: binning-prokaryotic.py -f <scaffolds.fasta> -b <mapped.sorted.bam> -n <name> -o <output_directory>

    Running: binning-prokaryotic.py v2022.7.8 via Python v3.8.5 | /Users/jespinoz/anaconda3/bin/python

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
  --maxbin2_options MAXBIN2_OPTIONS
                        MaxBin2 | More options (e.g. --arg 1 ) [Default: ''] | https://sourceforge.net/projects/maxbin/
  --metabat2_options METABAT2_OPTIONS
                        MetaBat2 | More options (e.g. --arg 1 ) [Default: ''] | https://bitbucket.org/berkeleylab/metabat/src/master/
  --concoct_options CONCOCT_OPTIONS
                        CONCOCT | More options (e.g. --arg 1 ) [Default: '']

Gene model arguments:
  --prodigal_genetic_code PRODIGAL_GENETIC_CODE
                        Prodigal -g translation table [Default: 11]

Evaluation arguments:
  --dastool_searchengine DASTOOL_SEARCHENGINE
                        DAS_Tool searchengine. [Default: diamond] | https://github.com/cmks/DAS_Tool
  --dastool_options DASTOOL_OPTIONS
                        DAS_Tool | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/cmks/DAS_Tool
  --pplacer_threads PPLACER_THREADS
                        Number of threads used for pplacer. Multithreaded uses a lot memory so don't use unless you have the resources [Default: 1]
  --checkm_tree CHECKM_TREE
                        CheckM tree type either 'reduced' or 'full' [Default: reduced]
  --checkm_completeness CHECKM_COMPLETENESS
                        CheckM completeness threshold [Default: 50]
  --checkm_contamination CHECKM_CONTAMINATION
                        CheckM contamination threshold [Default: 10]
  --checkm_strain_heterogeneity CHECKM_STRAIN_HETEROGENEITY
                        CheckM strain hetereogeneity threshold
  --checkm_options CHECKM_OPTIONS
                        CheckM lineage_wf | More options (e.g. --arg 1 ) [Default: '']
  --gtdbtk_options GTDBTK_OPTIONS
                        GTDB-Tk | classify_wf options (e.g. --arg 1 ) [Default: '']

featureCounts arguments:
  --featurecounts_options FEATURECOUNTS_OPTIONS
                        featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/

Domain classification arguments:
  --logit_transform LOGIT_TRANSFORM
                         Transformation for consensus_domain_classification: {softmax, tss} [Default: softmax
  --tiara_minimum_length TIARA_MINIMUM_LENGTH
                        Tiara | Minimum contig length. Anything lower than 3000 is not recommended. [Default: 3000]
  --tiara_options TIARA_OPTIONS
                        Tiara | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/ibe-uw/tiara

```

**Output:**

* binned.list - List of binned contigs
* bins.list - List of MAG identifiers
* checkm_output.filtered.tsv - Filtered CheckM output
* featurecounts.orfs.tsv.gz - ORF-level counts table
* genome_statistics.tsv - Genome assembly statistics
* genomes/ - MAG subdirectory
* genomes/\*.fa - MAG assembly fasta
* genomes/\*.faa - MAG protein fasta
* genomes/\*.ffn - MAG CDS fasta
* genomes/\*.gff - MAG gene models
* genomes/identifier_mapping.tsv - Identifier mapping between [id_orf, id_contig, id_mag]
* gtdbtk_output.filtered.tsv - Filtered GTDBTk output
* scaffolds_to_bins.tsv - Identifier mapping between [id_contig, id_mag]
* unbinned.fasta - Fasta of unbinned contigs that have passed length thresholding
* unbinned.list - List of unbinned contigs

#### binning-eukaryotic – Binning for recovering eukaryotic genomes with exon-aware gene modeling and lineage-specific quality assessment
The eukaryotic binning module uses several checks and state-of-the-art software to ensure high quality genomes.  In particular, non-prokaryotic-biased binning algorithms MetaBat2 [default] (coverage calculated with CoverM) or CONCOCT (coverage calculated using CONCOCT scripts) is used for binning out genomes followed by a genome size filter (2,000,000 bp is the default). The preliminary bins are run through the consensus domain wrapper for Tiara to predict eukaryotic MAGs. We found Tiara to be much more effective than EukRep (even when using an exhaustive parameter grid for EukRep) that we tested on known microeukaryotic contigs/transcripts from published assemblies and transcriptomes (Karlicki et al., 2022; West, Probst, Grigoriev, Thomas, & Banfield, 2018). Contigs from the eukaryotic MAGs are input into MetaEuk easy-predict workflow (Levy Karin et al., 2020) using our custom consensus eukaryotic database (see Database section in Methods). Although MetaEuk is a high-quality software suite, the identifiers from MetaEuk are very complex, long, and contain characters that are often problematic for downstream applications including parsing, file naming systems, and certain programs with simplified identifier requirements such as Anvi’o (Eren et al., 2015). In addition, the gene model GFF files are not intuitive, compatible with Prodigal GFF files or featureCounts without major modification. Therefore, we developed an essential wrapper for MetaEuk that simplifies identifiers (i.e., [ContigID]\_[GeneStart]:\[GeneEnd]([strand])), ensuring no duplicates are produced, creates a GFF file that can be concatenated with the Prodigal GFF file for use with featureCounts, and several identifier mapping tables to seamless convert between original and modified identifiers. Lineage-specific genome quality estimation is performed using BUSCO (Manni et al., 2021) where poor quality MAGs are removed (e.g., completeness < 50%).  Gene counts are computed using featureCounts at the gene level.  Lastly, genome statistics such as N50, number of scaffolds, and genome size are calculated using seqkit. The input to this module is a fasta file (typically the unbinned.fasta from the prokaryotic binning module) and sorted BAM while the output includes the eukaryotic MAGs, gene models via MetaEuk, identifier mappings, BUSCO output, counts tables, and unbinned fasta.  Iterative binning is not currently available as no consensus binning tool is available therefore iterative binning would result in diminishing returns. MAG naming scheme for eukaryotes follows \[SAMPLE]\_\_\[ALGORITHM]\_\_E.[ITERATION]\_\_\[NAME] (e.g., ERR2002407\_\_METABAT2\_\_E.1\_\_bin.2).

**Conda Environment**: `conda activate VEBA-binning-eukaryotic_env`


```
usage: binning-eukaryotic.py -f <scaffolds.fasta> -b <mapped.sorted.bam> -n <name> -o <output_directory>

    Running: binning-eukaryotic.py v2022.7.8 via Python v3.8.5 | /Users/jespinoz/anaconda3/bin/python

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

BUSCO arguments:
  --busco_completeness BUSCO_COMPLETENESS
                        BUSCO completeness [Default: 50.0]
  --busco_contamination BUSCO_CONTAMINATION
                        BUSCO contamination [Default: 10.0]
  --busco_evalue BUSCO_EVALUE
                        BUSCO | E-value cutoff for BLAST searches. Allowed formats, 0.001 or 1e-03 [Default: 1e-03]

featureCounts arguments:
  --featurecounts_options FEATURECOUNTS_OPTIONS
                        featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/

```

**Output:**

* binned.list - List of binned contigs
* bins.list - List of MAG identifiers
* busco_results.filtered.tsv - Filtered BUSCO output
* featurecounts.orfs.tsv.gz - ORF-level counts table
* genome_statistics.tsv - Genome assembly statistics
* genomes/ - MAG subdirectory
* genomes/\*.fa - MAG assembly fasta
* genomes/\*.faa - MAG protein fasta
* genomes/\*.ffn - MAG CDS fasta
* genomes/\*.gff - MAG gene models
* genomes/identifier_mapping.tsv - Identifier mapping between [id_orf, id_contig, id_mag]
* identifier_mapping.metaeuk.tsv - Identifier mapping between original MetaEuk identifiers and modified identifiers.  Includes fully parsed MetaEuk identifiers.
* scaffolds_to_bins.tsv - Identifier mapping between [id_contig, id_mag]
* unbinned.fasta - Fasta of unbinned contigs that have passed length thresholding
* unbinned.list - List of unbinned contigs

#### binning-viral – Detection of viral genomes and quality assessment
Viral binning is performed using VirFinder (Ren et al., 2017) to extract potential viral contigs (e.g., P < 0.05). The potential viral contigs are then input into CheckV (Nayfach et al., 2020) where quality assessment removes poor quality or low confidence viral predictions. The filtering scheme is based on author recommendations (Nayfach, 2021) in which a viral contig is considered if it meets the following criteria: 1) number of viral genes ≥ 5 x number of host genes; 2) completeness ≥ 50%; 3) CheckV quality is either medium-quality, high-quality, or complete; and 4) MIUViG quality is either medium-quality, high-quality, or complete (Roux et al., 2018).  Proviruses can be included by using the --include_proviruses flag. After poor quality viral contigs are removed, Prodigal is used for gene modeling and seqkit is used for useful genome statistics.  The input to this module is a fasta file (typically the unbinned.fasta from the eukaryotic binning module) while the output includes the viral MAGs, gene models via Prodigal, identifier mappings, and CheckV output.  Iterative binning is not applicable for viral detection as algorithms are executed on a per-contig basis and all viral genomes will be identified on first pass.  MAG naming scheme for viruses follows [SAMPLE]\_\_\[ALGORITHM]\_\_\[NAME] (e.g., SRR9668957\_\_VIRFINDER\_\_Virus.1).

**Conda Environment**: `conda activate VEBA-binning-viral_env`

```
usage: binning-viral.py -f <scaffolds.fasta> -l <contig_identifiers> -n <name> -o <output_directory>

    Running: binning-viral.py v2022.7.8 via Python v3.8.5 | /Users/jespinoz/anaconda3/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -f FASTA, --fasta FASTA
                        path/to/scaffolds.fasta
  -n NAME, --name NAME  Name of sample
  -l CONTIG_IDENTIFIERS, --contig_identifiers CONTIG_IDENTIFIERS
                        path/to/contigs.list [Optional]
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
  --remove_temporary_fasta
                        If contig identifiers were provided and a fasta is generated, remove this file

Database arguments:
  --veba_database VEBA_DATABASE
                        VEBA database location.  [Default: $VEBA_DATABASE environment variable]

Binning arguments:
  -m MINIMUM_CONTIG_LENGTH, --minimum_contig_length MINIMUM_CONTIG_LENGTH
                        Minimum contig length.  [Default: 1500]

Gene model arguments:
  --prodigal_genetic_code PRODIGAL_GENETIC_CODE
                        Prodigal -g translation table [Default: 11]

Virus arguments:
  --include_provirus    Include provirus viral detection
  --virfinder_pvalue VIRFINDER_PVALUE
                        VirFinder p-value threshold [Default: 0.05]
  --virfinder_options VIRFINDER_OPTIONS
                        VirFinder | More options (e.g. --arg 1 ) [Default: '']
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

```

**Output:**

* binned.list - List of binned contigs
* bins.list - List of MAG identifiers
* quality_summary.filtered.tsv - Filtered CheckV output
* featurecounts.orfs.tsv.gz - ORF-level counts table (If --bam file is provided)
* genome_statistics.tsv - Genome assembly statistics
* genomes/ - MAG subdirectory
* genomes/\*.fa - MAG assembly fasta
* genomes/\*.faa - MAG protein fasta
* genomes/\*.ffn - MAG CDS fasta
* genomes/\*.gff - MAG gene models
* genomes/identifier_mapping.tsv - Identifier mapping between [id_orf, id_contig, id_mag]
* scaffolds_to_bins.tsv - Identifier mapping between [id_contig, id_mag]
* unbinned.fasta - Fasta of unbinned contigs that have passed length thresholding
* unbinned.list - List of unbinned contigs


#### classify-prokaryotic – Taxonomic classification and candidate phyla radiation adjusted quality assessment of prokaryotic genomes
The prokaryotic classification module is a useful wrapper around GTDB-Tk which either combines the resulting archaea and bacteria summary tables or runs GTDB-Tk lineage_wf from the beginning.  If genome clusters are provided, then it performs consensus lineage classification.

```
usage: classify-prokaryotic.py -i <prokaryotic_binning_directory> -o <output_directory>

    Running: classify-prokaryotic.py v2022.06.07 via Python v3.8.5 | /Users/jespinoz/anaconda3/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -i PROKARYOTIC_BINNING_DIRECTORY, --prokaryotic_binning_directory PROKARYOTIC_BINNING_DIRECTORY
                        path/to/prokaryotic_binning_directory
  -c CLUSTERS, --clusters CLUSTERS
                        path/to/clusters.tsv, Format: [id_mag]<tab>[id_cluster], No header.
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        path/to/output_directory [Default: veba_output/classify/prokaryotic]
  -x EXTENSION, --extension EXTENSION
                        Fasta file extension for genomes if a list is provided [Default: fa]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  --tmpdir TMPDIR       path/to/TMPDIR [Default: /var/folders/fl/bh6r17mn52d33ycwvk70z2bc0003c0/T/]
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
  --gtdbtk_options GTDBTK_OPTIONS
                        GTDB-Tk | classify_wf options (e.g. --arg 1 ) [Default: '']

Consensus genome classification arguments:
  -l LENIENCY, --leniency LENIENCY
                        Leniency parameter. Lower value means more conservative weighting. A value of 1 indiciates no weight bias. A value greater than 1 puts higher weight on higher level taxonomic assignments. A value less than 1 puts lower weights on higher level taxonomic assignments.  [Default: 1.382]

**Conda Environment**: `conda activate VEBA-classify_env`

```


**Output:**

* prokaryotic_taxonomy.tsv - Prokaryotic genome classification based on GTDBTk
* prokaryotic_taxonomy.clusters.tsv - Prokaryotic cluster classification (If --clusters are provided)
                        
#### classify-eukaryotic – Taxonomic classification of eukaryotic genomes
The eukaryotic classification module utilizes the target field of MetaEuk gene identifiers and the taxonomic lineage associated with each source genome.  The default marker set is eukaryote_odb10 from BUSCO but custom marker sets are support along with the inclusion of all genes not just marker genes.  An option to include marker-specific noise cutoff scores is also available using the --scores_cutoff parameter which is default behavior with BUSCO’s eukaryote_odb10 provided noise thresholds.  For each MAG, bitscores are accumulated for each taxonomic level and taxonomy is assigned with leniency specified by the leniency parameter with high leniency resulting higher order taxonomic assignments.  If genome clusters are provided, then it performs consensus lineage classification. EUKulele (Krinos, Hu, Cohen, & Alexander, 2021) was attempted for this stage but a custom database could not be created with the software, likely due to the dependency of supergroup and division fields that were missing for certain taxa.  However, future implementations of VEBA may use this if such issues are resolved.  

**Conda Environment**: `conda activate VEBA-classify_env` 

```
usage: classify-eukaryotic.py -i <eukaryotic_binning_directory> -o <output_directory>

    Running: classify-eukaryotic.py v2022.7.8 via Python v3.8.5 | /Users/jespinoz/anaconda3/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -i EUKARYOTIC_BINNING_DIRECTORY, --eukaryotic_binning_directory EUKARYOTIC_BINNING_DIRECTORY
                        path/to/eukaryotic_binng_directory
  -c CLUSTERS, --clusters CLUSTERS
                        path/to/clusters.tsv, Format: [id_mag]<tab>[id_cluster], No header.
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        path/to/output_directory [Default: veba_output/classify/eukaryotic]

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
  --include_all_genes   Use if you want to include all genes for taxonomy classification instead of only core markers from BUSCO's eukaryota_odb10

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

* eukaryotic_taxonomy.tsv - Eukaryotic genome classification based on microeukaryotic protein database and BUSCO's eukaryota_odb10 marker set
* eukaryotic_taxonomy.clusters.tsv - Eukaryotic cluster classification (If --clusters are provided)
* gene-source_lineage.tsv - Gene source lineage and scores for classifying MAGs [id_gene, id_scaffold, id_mag, id_target, id_source, lineage, bitscore]


#### classify-viral – Taxonomic classification and isolation source of viral genomes
The viral classification module utilizes the CheckV database along with the best hit lineage and source habitat information from the CheckV output.  This includes a look up of CheckV identifiers based on direct terminal repeats and GenBank identifiers when applicable. If genome clusters are provided, then it performs consensus lineage classification and consensus habitat annotation. 

**Conda Environment**: `conda activate VEBA-classify_env`

```
usage: classify-viral.py -i <viral_binning_directory>  -o <output_directory>

    Running: classify-viral.py v2022.7.8 via Python v3.8.5 | /Users/jespinoz/anaconda3/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -i VIRAL_BINNING_DIRECTORY, --viral_binning_directory VIRAL_BINNING_DIRECTORY
                        Either: path/to/checkv/quality_summary.tsv or directory of veba_output/binning/viral
  -c CLUSTERS, --clusters CLUSTERS
                        path/to/clusters.tsv, Format: [id_mag]<tab>[id_cluster], No header.
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        path/to/output_directory [Default: veba_output/classify/viral]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit

Database arguments:
  --veba_database VEBA_DATABASE
                        VEBA database location.  [Default: $VEBA_DATABASE environment variable]

Consensus habitat/isolation source arguments:
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

* viral_taxonomy.tsv - Viral genome classification based on CheckV classifications
* viral_taxonomy.clusters.tsv - Viral cluster classification (If --clusters are provided)
* viral_isolation_source.clusters.tsv - Viral isolation source consensus (If --clusters are provided)

#### cluster – Species-level clustering of genomes and lineage-specific orthogroup detection
To leverage intra-sample genome analysis in an inter-sample analytical paradigm, genome clustering and lineage-specific orthogroup detection is necessary.  The cluster module first uses FastANI (Jain, Rodriguez-R, Phillippy, Konstantinidis, & Aluru, 2018) to compute pairwise ANI and these are used to construct a NetworkX graph object where nodes are genomes and edges are ANI values (Hagberg, Schult, & Swart, 2008).  This graph is converted into subgraphs of connected components whose edges are connected by a particular threshold such as 95% ANI [default] as recommended by the authors for species-level clustering.  These species-level clusters (SLC) are then partitioned and OrthoFinder (Emms & Kelly, 2019) is then run on each SLC panproteome. The input is a list of genome paths and list of protein fasta paths while the output includes identifier mappings between genomes, SLCs, scaffolds, proteins, and orthogroups.

**Conda Environment**: `conda activate VEBA-cluster_env`

```
usage: cluster.py -m <mags> -a <proteins> -o <output_directory> -t 95

    Running: cluster.py v2022.02.24 via Python v3.9.10 | /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-mapping_env/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -i SCAFFOLDS_TO_BINS, --scaffolds_to_bins SCAFFOLDS_TO_BINS
                        path/to/scaffolds_to_bins.tsv, Format: [id_scaffold]<tab>[id_bin], No header
  -m MAGS, --mags MAGS  Tab-seperated value table of [id_mag]<tab>[path/to/genome.fasta]
  -a PROTEINS, --proteins PROTEINS
                        Tab-seperated value table of [id_mag]<tab>[path/to/protein.fasta]
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        path/to/project_directory [Default: veba_output/cluster]
  --mags_extension MAGS_EXTENSION
                        Fasta file extension for --mags if a list is provided [Default: fa]
  --proteins_extension PROTEINS_EXTENSION
                        Fasta file extension for proteins if a list is provided [Default: faa]

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  -p N_JOBS, --n_jobs N_JOBS
                        Number of threads [Default: 1]
  --cluster_only        Only run FastANI
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit

FastANI arguments:
  -t ANI_THRESHOLD, --ani_threshold ANI_THRESHOLD
                        FastANI | Species-level clustering threshold [Default: 95.0]
  --fastani_options FASTANI_OPTIONS
                        FastANI | More options (e.g. --arg 1 ) [Default: '']
  --cluster_prefix CLUSTER_PREFIX
                        Cluster prefix [Default: 'SLC
  --cluster_suffix CLUSTER_SUFFIX
                        Cluster suffix [Default: '
  --copy_proteins       Copy instead of symlink
  --clone_label CLONE_LABEL
                        Singleton clone label [Default: '___clone

OrthoFinder arguments:
  --orthofinder_options ORTHOFINDER_OPTIONS
                        OrthoFinder | More options (e.g. --arg 1 ) [Default: '']

```
**Output:**

* clusters.tsv - Identifier mapping between MAGs and clusters [id_mag, id_cluster]
* identifier_mapping.orthogroups.tsv - Identifier mapping between ORFs, contigs, orthogroups, and clusters
* proteins_to_orthogroups.tsv - Identifier mapping between ORFs and orthogroups [id_orf, id_orthogroup]
* scaffolds_to_clusters.tsv - Identifier mapping between contigs and clusters [id_contig, id_cluster]

*Note: This should be run separately for prokaryotes, eukaryotes, and viruses. 

* veba_output/cluster/prokaryotes
* veba_output/cluster/prokaryotes
* veba_output/cluster/prokaryotes

#### annotate – Annotates translated gene calls against NR, Pfam, and KOFAM
Annotation is performed using best hit annotations and profile HMMs.  First proteins are aligned against NCBI non-redundant protein database (other databases are supported) using Diamond (Buchfink, Reuter, & Drost, 2021; Buchfink, Xie, & Huson, 2014). After annotation, protein domains are identified using the Pfam database (Mistry et al., 2021) via HMMER (Mistry, Finn, Eddy, Bateman, & Punta, 2013) and KEGG orthology is characterized via KOFAMSCAN (Aramaki et al., 2020).  Note, the `lineage_predictions.*.tsv` files generated here are experimental.

**Conda Environment**: `conda activate VEBA-annotate_env`

```
usage: annotate.py -i <identifier_mapping> -a <proteins> -o <output_directory>

    Running: annotate.py v2021.7.8 via Python v3.8.5 | /Users/jespinoz/anaconda3/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -a PROTEINS, --proteins PROTEINS
                        Either path/to/proteins.faa or a directory of fasta files using [-x]
  -i IDENTIFIER_MAPPING, --identifier_mapping IDENTIFIER_MAPPING
                        Tab-seperated value table of [id_orf]<tab>[id_contig]<tab>[id_mag]
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        path/to/project_directory [Default: veba_output/annotation]
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
  --hmmsearch_threshold HMMSEARCH_THRESHOLD
                        HMMSearch | Threshold {cut_ga, cut_nc, gut_tc} [Default:  cut_ga]
  --hmmsearch_options HMMSEARCH_OPTIONS
                        Diamond | More options (e.g. --arg 1 ) [Default: '']

KOFAMSCAN arguments:
  --kofamscan_options KOFAMSCAN_OPTIONS
                        Diamond | More options (e.g. --arg 1 ) [Default: '']

```

**Output:**

* annotations.orfs.tsv.gz - Concatenated annotations from Diamond (NR), HMMSearch (Pfam), and KOFAMSCAN (KOFAM)
* lineage_predictions.contigs.tsv.gz - [Experimental] Lineage predictions for contigs based on NR annotations.  This should be used only for experimentation as many lineages aren't defined properly. Contig-level classifications.
* lineage_predictions.mags.tsv.gz - [Experimental] Lineage predictions for contigs based on NR annotations.  This should be used only for experimentation as many lineages aren't defined properly. MAG-level classifications.



#### phylogeny – Constructs phylogenetic trees given a marker set
The phylogeny module is a tool used for phylogenetic inference and constructing phylogenetic trees for genomes given a reference marker set (see Databases section of Methods). This is performed by the following method: 1) identifying marker proteins using HMMSearch from the HMMER3 suite; 2) creating protein alignments for each marker identified MUSCLE (Edgar, 2004); 3) trimming the alignments using ClipKIT (Steenwyk, Buida, Li, Shen, & Rokas, 2020); 4) concatenating the alignments; 5) approximately-maximum-likelihood phylogenetic inference using FastTree2 (Price, Dehal, & Arkin, 2010); and 6) optional maximum likelihood phylogenetic inference using IQ-TREE2 (Minh et al., 2020).  An option to include marker-specific noise cutoff scores is also available using the --scores_cutoff parameter.  Poor-quality genomes that do not meet a threshold in the proportion of markers in the reference are removed using the --minimum_markers_aligned_ratio parameter.  Similarly, non-informative markers that are not prevalent in the query genomes are removed using the --minimum_genomes_aligned_ratio parameter.

**Conda Environment**: `conda activate VEBA-phylogeny_env`


```
usage: phylogeny.py -d <database_hmms> -a <proteins> -o <output_directory>

    Running: phylogeny.py v2022.06.22 via Python v3.8.5 | /Users/jespinoz/anaconda3/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -d DATABASE_HMM, --database_hmm DATABASE_HMM
                        path/to/HMM database of markers
  -a PROTEINS, --proteins PROTEINS
                        Can be the following format: 1) Tab-seperated value table of [id_mag]<tab>[path/to/protein.fasta]; 2) List of filepaths [path/to/protein.fasta]; or 3) Directory of protein fasta using --protein_extension
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
                        HMMSearch | Threshold {cut_ga, cut_nc, gut_tc, e} [Default:  e]
  --hmmsearch_evalue HMMSEARCH_EVALUE
                        Diamond | E-Value [Default: 10.0]
  --hmmsearch_options HMMSEARCH_OPTIONS
                        Diamond | More options (e.g. --arg 1 ) [Default: '']
  -f HMM_MARKER_FIELD, --hmm_marker_field HMM_MARKER_FIELD
                        HMM reference type (accession, name) [Default: accession
  -s SCORES_CUTOFF, --scores_cutoff SCORES_CUTOFF
                        path/to/scores_cutoff.tsv. No header. [id_hmm]<tab>[score]

Alignment arguments:
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



#### index – Builds local or global index for alignment to genomes
The index module creates reference indices for alignments in both local or global paradigms. In the local paradigm, an index is created for all the assembled genomes concatenated together for each sample. This is useful in situations where perfectly paired metagenomics and metatranscriptomics are available where the metatranscriptomics can be mapped directly to the de novo reference generated from the metagenomics.  However, this is not applicable in all cases such as when there is not a perfect overlap between metagenomics and metatranscriptomics.  In this global paradigm, assembled genomes are concatenated across all samples and an alignment index is created for this concatenated reference.  Currently, Bowtie2 (Langmead & Salzberg, 2012)  is the only alignment software packages supported. 

**Conda Environment**: `conda activate VEBA-index_env`

```
usage: index.py -i <mags> -o <output> --heatmap_output <pdf>

    Running: index.py v2022.02.17 via Python v3.9.10 | /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-mapping_env/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required I/O arguments:
  -r REFERENCES, --references REFERENCES
                        local mode: [id_sample]<tab>[path/to/reference.fa] and global mode: [path/to/reference.fa]
  -g GENE_MODELS, --gene_models GENE_MODELS
                        local mode: [id_sample]<tab>[path/to/reference.gff] and global mode: [path/to/reference.gff]
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
* [id_sample]/reference.fa.gz - Concatenated reference fasta
* [id_sample]/reference.fa.gz.\*.bt2 - Bowtie2 index of reference fasta
* [id_sample]/reference.gff - Concatenated gene models
* [id_sample]/reference.saf - SAF format for reference

#### mapping – Aligns reads to local or global index of genomes
The mapping module uses local or global reference indices generated by the index module and aligns reads using Bowtie2.  The alignment files are sorted to produce sorted BAM files using Samtools. Reads from the sorted BAM files are then fed into featureCounts to produce gene-level counts, orthogroup-level counts, and SLC-level counts. 

**Conda Environment**: `conda activate VEBA-mapping_env`

```
usage: mapping.py -1 <reads_1.fq> -2 <reads_2.fq> -n <name> -o <output_directory> -x <reference_directory>

    Wrapper around github.com/jolespin/fastq_preprocessor
    Running: mapping.py v2021.12.10 via Python v3.9.10 | /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/VEBA-mapping_env/bin/python

optional arguments:
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
  --orfs_to_orthogroups ORFS_TO_ORTHOGROUPS
                        path/to/orf_to_orthogroup.tsv, [id_orf]<tab>[id_orthogroup], No header
  --scaffolds_to_clusters SCAFFOLDS_TO_CLUSTERS
                        path/to/scaffold_to_cluster.tsv, [id_scaffold]<tab>[id_cluster], No header

```

**Output:**

* counts.orfs.tsv.gz - ORF-level counts table
* counts.scaffolds.tsv.gz - Contig-level counts table
* mapped.sorted.bam - Sorted BAM file
* unmapped_1.fastq.gz - Unmapped reads (forward)
* unmapped_2.fastq.gz - Unmapped reads (reverse)

___________________________________________________________________
#### References:


Alneberg, J., Bjarnason, B.S., De Bruijn, I., Schirmer, M., Quick, J., Ijaz, U.Z., Lahti, L., Loman, N.J., Andersson, A.F., and Quince, C. (2014). Binning metagenomic contigs by coverage and composition. Nat. Methods 11, 1144–1146.

Antipov, D., Raiko, M., Lapidus, A., and Pevzner, P.A. (2020). Metaviral SPAdes: assembly of viruses from metagenomic data. Bioinformatics 36, 4126–4129.

Aramaki, T., Blanc-Mathieu, R., Endo, H., Ohkubo, K., Kanehisa, M., Goto, S., and Ogata, H. (2020). KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold. Bioinformatics 36, 2251–2252.

Bankevich, A., Nurk, S., Antipov, D., Gurevich, A.A., Dvorkin, M., Kulikov, A.S., Lesin, V.M., Nikolenko, S.I., Pham, S., Prjibelski, A.D., et al. (2012). SPAdes: a new genome assembly algorithm and its applications to single-cell  sequencing. J. Comput. Biol.  a J. Comput. Mol. Cell  Biol. 19, 455–477.

Buchfink, B., Xie, C., and Huson, D.H. (2014). Fast and sensitive protein alignment using DIAMOND. Nat. Methods 2014 121 12, 59–60.

Buchfink, B., Reuter, K., and Drost, H.G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. Nat. Methods 2021 184 18, 366–368.

Bushmanova, E., Antipov, D., Lapidus, A., and Prjibelski, A.D. (2019). rnaSPAdes: a de novo transcriptome assembler and its application to RNA-Seq data. Gigascience 8, 1–13.

Chaumeil, P.-A., Mussig, A.J., Hugenholtz, P., and Parks, D.H. (2020). GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. Bioinformatics 36, 1925–1927.

Chen, S., Zhou, Y., Chen, Y., and Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, i884–i890.

Edgar, R.C. (2004). MUSCLE: A multiple sequence alignment method with reduced time and space complexity. BMC Bioinformatics 5, 1–19.

Emms, D.M., and Kelly, S. (2019). OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biol. 2019 201 20, 1–14.

Eren, A.M., Esen, O.C., Quince, C., Vineis, J.H., Morrison, H.G., Sogin, M.L., and Delmont, T.O. (2015). Anvi’o: An advanced analysis and visualization platformfor ’omics data. PeerJ 2015, e1319.

Hagberg, A.A., Schult, D.A., and Swart, P.J. (2008). Exploring Network Structure, Dynamics, and Function using NetworkX.

Hyatt, D., Chen, G.-L., Locascio, P.F., Land, M.L., Larimer, F.W., and Hauser, L.J. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119.

Jain, C., Rodriguez-R, L.M., Phillippy, A.M., Konstantinidis, K.T., and Aluru, S. (2018). High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. Nat. Commun. 2018 91 9, 1–8.

Kang, D.D., Li, F., Kirton, E., Thomas, A., Egan, R., An, H., and Wang, Z. (2019). MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. PeerJ 7.

Karlicki, M., Antonowicz, S., and Karnkowska, A. (2022). Tiara: deep learning-based classification system for eukaryotic sequences. Bioinformatics 38, 344–350.

Krinos, A., Hu, S., Cohen, N., and Alexander, H. (2021). EUKulele: Taxonomic annotation of the unsung eukaryotic microbes. J. Open Source Softw. 6, 2817.

Langmead, B., and Salzberg, S.L. (2012). Fast gapped-read alignment with Bowtie 2. Nat. Methods 9, 357.

Levy Karin, E., Mirdita, M., and Söding, J. (2020). MetaEuk-sensitive, high-throughput gene discovery, and annotation for large-scale eukaryotic metagenomics. Microbiome 8, 1–15.

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., and Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078–2079.
Liao, Y., Smyth, G.K., and Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics 30, 923–930.

Manni, M., Berkeley, M.R., Seppey, M., Simão, F.A., and Zdobnov, E.M. (2021). BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. Mol. Biol. Evol. 38, 4647–4654.

Meleshko, D., Hajirasouliha, I., and Korobeynikov, A. (2021). coronaSPAdes: from biosynthetic gene clusters to RNA viral assemblies. Bioinformatics.

Minh, B.Q., Schmidt, H.A., Chernomor, O., Schrempf, D., Woodhams, M.D., Von Haeseler, A., Lanfear, R., and Teeling, E. (2020). IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. Mol. Biol. Evol. 37, 1530–1534.

Mistry, J., Finn, R.D., Eddy, S.R., Bateman, A., and Punta, M. (2013). Challenges in homology search: HMMER3 and convergent evolution of coiled-coil regions. Nucleic Acids Res. 41, e121–e121.

Mistry, J., Chuguransky, S., Williams, L., Qureshi, M., Salazar, G.A., Sonnhammer, E.L.L., Tosatto, S.C.E., Paladin, L., Raj, S., Richardson, L.J., et al. (2021). Pfam: The protein families database in 2021. Nucleic Acids Res. 49, D412–D419.

Nayfach, S. (2021).  Recommended cutoffs for analyzing CheckV results?

Nayfach, S., Camargo, A.P., Schulz, F., Eloe-Fadrosh, E., Roux, S., and Kyrpides, N.C. (2020). CheckV assesses the quality and completeness of metagenome-assembled viral genomes. Nat. Biotechnol. 2020 395 39, 578–585.

Nurk, S., Meleshko, D., Korobeynikov, A., and Pevzner, P.A. (2017). metaSPAdes: a new versatile metagenomic assembler. Genome Res. 27, 824–834.

Price, M.N., Dehal, P.S., and Arkin, A.P. (2010). FastTree 2--approximately maximum-likelihood trees for large alignments. PLoS One 5, e9490.

Ren, J., Ahlgren, N.A., Lu, Y.Y., Fuhrman, J.A., and Sun, F. (2017). VirFinder: a novel k-mer based tool for identifying viral sequences from assembled  metagenomic data. Microbiome 5, 69.

Roux, S., Adriaenssens, E.M., Dutilh, B.E., Koonin, E. V., Kropinski, A.M., Krupovic, M., Kuhn, J.H., Lavigne, R., Brister, J.R., Varsani, A., et al. (2018). Minimum Information about an Uncultivated Virus Genome (MIUViG). Nat. Biotechnol. 2018 371 37, 29–37.

Shen, W., Le, S., Li, Y., and Hu, F. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS One 11, e0163962.

Sieber, C.M.K., Probst, A.J., Sharrar, A., Thomas, B.C., Hess, M., Tringe, S.G., and Banfield, J.F. (2018). Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. Nat. Microbiol. 2018 37 3, 836–843.

Steenwyk, J.L., Buida, T.J., Li, Y., Shen, X.X., and Rokas, A. (2020). ClipKIT: A multiple sequence alignment trimming software for accurate phylogenomic inference. PLOS Biol. 18, e3001007.

West, P.T., Probst, A.J., Grigoriev, I. V., Thomas, B.C., and Banfield, J.F. (2018). Genome-reconstruction for eukaryotes from complex natural microbial communities. Genome Res. 28, 569–580.

Wu, Y.-W., Simmons, B.A., and Singer, S.W. (2016). MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. Bioinformatics 32, 605–607.

