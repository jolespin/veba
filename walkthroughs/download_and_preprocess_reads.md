### Download and preprocess reads
This walkthrough is for downloading reads from NCBI and running *VEBA's* `preprocess.py` module to decontaminate either metagenomic or metatranscriptomic reads.  Skip step 3 if you already have reads (e.g., a new project):

#### 1. Activate the preprocess environment 
Since *VEBA* has so many capabilities I have partitioned different walkthroughs into modules and environments for those modules.  In this case, we are going to use the `VEBA-preprocess_env`.  For assembly, we will use `VEBA-assembly_env`, etc.  This is to increase our flexibility in updating dependencies.
```
conda activate VEBA-preprocess_env
```

If you want to either remove human contamination or count ribosomal reads then make sure VEBA database in your path.  If it's not there then `export VEBA_DATABASE=/path/to/database/`.  If you just want to quality trim/remove adapters then you don't need the databases.

```
echo $VEBA_DATABASE
/expanse/projects/jcl110/db/veba/VDB_v6
# ^_^ Yours will be different obviously #
```

#### 2. Make a list of identifiers and create necessary directories
*VEBA* is built on structured data so if you're downloading from NCBI SRA then this should be really easy.  If not, you'll just have to explicitly get each forward read, reverse read, and name.  Let's just assume you're downloading from SRA for the sake of this walkthrough.

Here are the SRA identifiers for `PRJNA777294` which is the Plastisphere dataset from the paper:

```
mkdir -p logs/
mkdir -p Fastq/
cat identifiers.list
```


<details open>
<summary>Output:</summary>

```
SRR17458603
SRR17458604
SRR17458605
SRR17458606
SRR17458607
SRR17458608
SRR17458609
SRR17458610
SRR17458611
SRR17458612
SRR17458613
SRR17458614
SRR17458615
SRR17458616
SRR17458617
SRR17458618
SRR17458619
SRR17458620
SRR17458621
SRR17458622
SRR17458623
SRR17458624
SRR17458625
SRR17458626
SRR17458627
SRR17458628
SRR17458629
SRR17458630
SRR17458631
SRR17458632
SRR17458633
SRR17458634
SRR17458635
SRR17458636
SRR17458637
SRR17458638
SRR17458639
SRR17458640
SRR17458641
SRR17458642
SRR17458643
SRR17458644
SRR17458645
SRR17458646
```

</details>

#### 3. Create a directory for the raw Fastq reads and download the sequences from NCBI

**Note:** `kingfisher` is not officially part of the `VEBA` suite and is only provided for convenience.

```
# Activate environment
conda activate VEBA-preprocess_env

# Set the threads you want to use.  Here we are using 4.
N_JOBS=4

# Go into fastq directory temporarily
cd Fastq/

# Iterate through the identifier list and download each read set
for ID in $(cat ../identifiers.list);
	do kingfisher get -r $ID -m aws-http -t ${N_JOBS} -f fastq.gz
	done

# Get back out to main directory
cd ..
```

If you get an error with `kingfisher`, check out [#14 on the FAQ](https://github.com/jolespin/veba/blob/main/FAQ.md#14-why-am-i-getting-errors-for-kingfisher-saying-a-command-was-not-found).

#### 4. Perform quality/adapter trimming, remove human contamination, and count the ribosomal reads but don't remove them.

Here we are going to count the reads for the human contamination and ribosomal reads but not keep them.  If you wanted to keep the human reads to do some type of human-based study then do `--retain_contaminated_reads 1`.  If you wanted to just do read trimming and not remove any contamination or count anything then you would just leave out the `-x` and `-k` arguments.  For example, if you were using this to process some human reads.  


* ⚠️ If your host is not human then you will need to use a different contamination reference.  See item #22 in the [FAQ](https://github.com/jolespin/veba/blob/main/FAQ.md).


```
N_JOBS=4

# CHM13v2.0
HUMAN_INDEX=${VEBA_DATABASE}/Contamination/chm13v2.0/chm13v2.0

# Ribosomal k-mer fasta
RIBOSOMAL_KMERS=${VEBA_DATABASE}/Contamination/kmers/ribokmers.fa.gz

# Iterate through identifiers and run preprocessing
for ID in  $(cat identifiers.list); do

	# Get forward and reverse reads
	R1=Fastq/${ID}_1.fastq.gz
	R2=Fastq/${ID}_2.fastq.gz
	
	# Get a name to use for log files
	N=preprocessing__${ID}
	
	# Remove any preexisting log file
	rm -f logs/${N}.*
	
	# Set up the command (use source from base environment instead of conda because of the `init` issues)
	CMD="source activate VEBA && veba --module preprocess --params \"-n ${ID} -1 ${R1} -2 ${R2} -p ${N_JOBS} -x ${HUMAN_INDEX} -k ${RIBOSOMAL_KMERS} --retain_contaminated_reads 0 --retain_kmer_hits 0 --retain_non_kmer_hits 0 -o veba_output/preprocess\""
	
	# If you have SunGrid engine, do something like this:
	# qsub -o logs/${N}.o -e logs/${N}.e -cwd -N ${N} -j y -pe threaded ${N_JOBS} "${CMD}"
	
	# If you have SLURM engine, do something like this:
	# sbatch -J ${N} -N 1 -c ${N_JOBS} --ntasks-per-node=1 -o logs/${N}.o -e logs/${N}.e --export=ALL -t 12:00:00 --mem=20G --wrap="${CMD}"
	
	done
```
Note: `preprocess` is a wrapper around `fastq_preprocessor`.

It creates the following directory structure where each sample is it's own subdirectory.  Makes globbing much easier:


```
ls veba_output/preprocess/

SRR17458603  SRR17458606  SRR17458609  SRR17458612  SRR17458615  SRR17458618  SRR17458621  SRR17458624  SRR17458627  SRR17458630  SRR17458633  SRR17458636  SRR17458639  SRR17458642  SRR17458645
SRR17458604  SRR17458607  SRR17458610  SRR17458613  SRR17458616  SRR17458619  SRR17458622  SRR17458625  SRR17458628  SRR17458631  SRR17458634  SRR17458637  SRR17458640  SRR17458643  SRR17458646
SRR17458605  SRR17458608  SRR17458611  SRR17458614  SRR17458617  SRR17458620  SRR17458623  SRR17458626  SRR17458629  SRR17458632  SRR17458635  SRR17458638  SRR17458641  SRR17458644
```

The main files in each of the output directories are the following: 

* cleaned_1.fastq.gz - Cleaned and trimmed fastq file (forward)
* cleaned_2.fastq.gz - Cleaned and trimmed fastq file (reverse)
* seqkit_stats.concatenated.tsv - Concatenated read statistics for all intermediate steps (e.g., fastp, bowtie2 removal of contaminated reads if provided, bbduk.sh removal of contaminated reads if provided)

#### 5. If you want some read statistics then merge the dataframes.

```
concatenate_dataframes.py  veba_output/preprocess/*/output/seqkit_stats.concatenated.tsv
```

#### 6. All the reads trimmed and decontaminated reads will be here:
```
veba_output/preprocess/*/output/cleaned_1.fastq.gz
veba_output/preprocess/*/output/cleaned_2.fastq.gz
```

#### Next steps:

Now it's time to assemble these reads, recover genomes from the assemblies, and map reads to produce counts tables.  Please see the next walkthroughs.