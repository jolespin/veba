### Coassembly setup for either metagenomics or metatranscriptomics
In some cases where you have low depth or sparse sequencing it may be preferable to go with a coassembly approach.  While this is not the main objective of *VEBA*, it is certainly supported. This walkthrough goes through setting up a coassembly including concatenating reads, the actual assembly, aligning reads to the coassembly to produce sorted BAM files, and counts tables.  The method will be the same for metagenomics and metatranscriptomics with the only difference being the spades program used.

What you'll end up with at the end of this is a coassembly, assembly statistics, and bam files for each sample  using the coassembly as reference.

Please refer to the [end-to-end metagenomics](end-to-end_metagenomics.md) or [recovering viruses from metatranscriptomics](recovering_viruses_from_metatranscriptomics.md) workflows for details on how to proceed with coassembly (e.g., binning, annotation, clustering, etc.).

_____________________________________________________

#### Steps:

1. Concatenate the forward reads
2. Concetenate the reverse reads
3. Coassembly using assembly.py
4. Align reads from each sample to the coassembly to create sorted BAM files that will be used for binning and counts tables.

#### 1. Concatenate forward and reverse reads separately

Refer to the [downloading and preprocessing reads workflow](download_and_preprocess_reads.md).  At this point, it's assumed you have the following: 

* A file with all of your identifiers on a separate line (e.g., `identifiers.list` but you can call it whatever you want)
* A directory to keep all your logs called `logs/`
* A directory of preprocessed reads: `veba_output/preprocess/${ID}/output/cleaned_1.fastq.gz` and `veba_output/preprocess/${ID}/output/cleaned_2.fastq.gz` where `${ID}` represents the identifiers in `identifiers.list`.

```bash
# Make a miscellaneous directory
mkdir -p veba_output/misc

# Compile reads table
# Note: Prior to 2022.10.24 `-a` and `--absolute` were used for absolute paths.  Now absolute paths are default but relative paths can be used via `-r` or `--relative` flags.  Older versions had a `-a` in the command below. 

compile_reads_table.py -i veba_output/preprocess/ > veba_output/misc/reads_table.tsv

# Forward reads
cat veba_output/preprocess/*/output/cleaned_1.fastq.gz > veba_output/misc/concatenated_1.fastq.gz

# Reverse reads
cat veba_output/preprocess/*/output/cleaned_2.fastq.gz > veba_output/misc/concatenated_2.fastq.gz

```

#### 2. Assemble reads, map reads to assembly, calculate assembly statistics, and index the assembly

Here we are going to coassemble all of the reads using `metaSPAdes` which is default but if you are using metatranscriptomics then use `--spades_program rnaSPAdes.py`.  

**Conda Environment:** `conda activate VEBA-assembly_env`

```
# Set the number of threads to use for each sample. Let's use 4
N_JOBS=4

# This is the default output directory 
OUT_DIR=veba_output/assembly

# Just making sure your log directory is available
mkdir -p logs 

# Run the assembly command
ID="coassembly"
	
# Set a useful name for log files
N="assembly__${ID}";
rm -f logs/${N}.*
	
# Get forward and reverse reads
R1=veba_output/misc/concatenated_1.fastq.gz
R2=veba_output/misc/concatenated_2.fastq.gz
	
# Set up command
CMD="source activate VEBA-assembly_env && assembly.py -1 ${R1} -2 ${R2} -n ${ID} -o ${OUT_DIR} -p ${N_JOBS}"

# Use this for metatranscriptomics
# CMD="source activate VEBA-assembly_env && assembly.py -1 ${R1} -2 ${R2} -n ${ID} -o ${OUT_DIR} -p ${N_JOBS} --spades_program rnaspades.py"
	
# Either run this command or use SunGridEnginge/SLURM
	

```

The following output files will produced for each sample: 

* featurecounts.tsv.gz - featureCounts output for contig-level counts
* mapped.sorted.bam - Sorted BAM
* mapped.sorted.bam.bai - Sorted BAM index
* scaffolds.fasta - Scaffold assemblies
* scaffolds.fasta.\*.bt2 - Bowtie2 index of scaffolds
* scaffolds.fasta.saf - SAF formatted file for scaffold-level counts with featureCounts
* seqkit_stats.tsv.gz - Assembly statistics

The main one we need is `scaffolds.fasta`  which we will use for binning.  Note that there is a `mapped.sorted.bam` file but we aren't going to use this because it loses all sample-specific information.  Instead, we are going to align sample-specific reads to the coassembly and use multiple sorted BAM files for binning instead.

#### 3. Align sample-specific reads to the coassembly



**Conda Environment:** `conda activate VEBA-assembly_env`

```
N_JOBS=4

ID="coassembly"
N="coverage__${ID}";
rm -f logs/${N}.*
FASTA=veba_output/assembly/${ID}/output/scaffolds.fasta
READS=veba_output/misc/reads_table.tsv
CMD="source activate VEBA-assembly_env && coverage.py -f ${FASTA} -r ${READS} -p ${N_JOBS} -m 1500 -o veba_output/coverage/${ID}"
	
# Either run this command or use SunGridEnginge/SLURM

```


The following output files will produced for each sample: 

* featurecounts.tsv.gz - featureCounts counts table of all samples
* [sample_id]/mapped.sorted.bam - Sorted BAM file under a subdirectory for each sample
* reference.fasta - Reference fasta (this is the filtered coassembly)
* reference.fasta.saf - SAF formatted file for contig-level counts with featureCounts
* seqkit_stats.tsv - Assembly statistics


_____________________________________________________

#### Next steps:

Now that you have a coassembly and multiple sorted BAM files, it's time for binning.  Start at step 3 of the [end-to-end metagenomics](end-to-end_metagenomics.md) or [recovering viruses from metatranscriptomics](recovering_viruses_from_metatranscriptomics.md) workflows depending on whether or not you have metagenomics or metatranscriptomics, respectively.  

**Please do not forget to adapt the BAM argument in the `binning-prokaryotic.py` command to include all the sample-specific sorted BAM files and not the concatenated sorted BAM.**  

More specifically, use `BAM="veba_output/coverage/coassembly/output/*/mapped.sorted.bam"` and not `BAM="veba_output/assembly/coassembly/output/mapped.sorted.bam"`.

Note that both will work but the former will yield much better results since coverage for different scaffolds/contigs are kept separate for each sample.
