#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, random, shutil
from collections import OrderedDict, defaultdict

import pandas as pd
import numpy as np

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2025.4.9"

def get_maxbin2_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Create dummy scaffolds_to_bins.tsv to overwrite later. 
    cmd = [
         "echo '' > {}".format(os.path.join(output_directory, "scaffolds_to_bins.tsv")),
         "&&",
         "mkdir -p {}".format(os.path.join(output_directory, "bins")),

    ]

    coverage_file = opts.coverage
    if opts.bam:
        coverage_file = os.path.join(output_directory,"intermediate", "coverage.tsv")
        
        # Coverage for MaxBin2
        cmd += [
        "&&",
        os.environ["coverm"],
        "contig",
        "--threads {}".format(opts.n_jobs),
         "--methods",
         "mean",
         "--bam-files",
        " ".join(opts.bam),
        "|",
        "tail",
        "-n",
        "+2",
        ">",
        coverage_file,
        ]

    cmd += [ 
        "&&",

        os.environ["run_MaxBin.pl"],
        "-contig {}".format(opts.fasta), # scaffolds.fasta
        "-out {}".format(os.path.join(output_directory, "intermediate", "bin")),
        "-abund {}".format(coverage_file ),
        "-min_contig_length {}".format(opts.minimum_contig_length),
        "-markerset {}".format(opts.maxbin2_markerset),
        "-thread {}".format(opts.n_jobs),
        "-verbose",
        opts.maxbin2_options,
    ]

    # Move bins and change extension
    cmd += [ 
"""

for FP in %s;
    do ID_GENOME=$(basename ${FP} .fasta);
    mv $FP %s/${ID_GENOME}.fa
    done

"""%( 
    os.path.join(output_directory, "intermediate", "*.fasta"),
    os.path.join(output_directory, "bins"),
    )
    ]

    # Remove small bins
    cmd += [ 
r"""
for FP in %s;
    do GENOME_SIZE=$(cat $FP | grep -v "^>" | tr -d "\n" | wc -m);
    echo "[GENOME SIZE] ${FP} = ${GENOME_SIZE}"

    if (( %d > ${GENOME_SIZE} )); then
        echo "[REMOVING] ${FP}"
        rm -f ${FP}
    fi

    done


"""%( 
    os.path.join(output_directory, "bins", "*.fa"),
    opts.minimum_genome_length,
    )
    
    ]

    return cmd

# Metabat2
def get_metabat2_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Create dummy scaffolds_to_bins.tsv to overwrite later. 
    cmd = [
         "echo '' > {}".format(os.path.join(output_directory, "scaffolds_to_bins.tsv")),
         "&&",
         "mkdir -p {}".format(os.path.join(output_directory, "bins")),

    ]

    # Coverage for MetaBat2 
    coverage_file = opts.coverage

    if opts.bam:
        coverage_file = os.path.join(output_directory,"intermediate", "coverage.tsv")
        cmd += [
        "&&",
        os.environ["coverm"],
        "contig",
        "--threads {}".format(opts.n_jobs),
         "--methods metabat",
         "--bam-files",
        " ".join(opts.bam),
        ">",
        coverage_file,
        ]
        
    cmd += [
        "&&",
       
        os.environ["metabat2"],
        "-i {}".format(opts.fasta), # scaffolds.fasta
        "-o {}".format(os.path.join(output_directory, "bins", "bin")),
        "-a {}".format(coverage_file),
        "-m {}".format(max(opts.minimum_contig_length, 1500)), # mininimum contig length must be >= 1500 nts for MetaBat2
        "-t {}".format(opts.n_jobs),
        "--minClsSize {}".format(opts.minimum_genome_length),
        # "--unbinned", # Don't forget about these...
        "--seed {}".format(opts.random_state),
        "--verbose",
        opts.metabat2_options,
        
        ]


    return cmd


def get_concoct_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [
        # Create dummy scaffolds_to_bins.tsv to overwrite later. 
        "echo '' > {}".format(os.path.join(output_directory, "scaffolds_to_bins.tsv")),
        "&&",
        "mkdir -p {}".format(os.path.join(output_directory, "bins")),
    ]
    
    cmd += [    
        "&&",
        # cut_up_fasta.py original_contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
        os.environ["cut_up_fasta.py"],
        opts.fasta, 
        "-c {}".format(opts.concoct_fragment_length),
        "-o {}".format(opts.concoct_overlap_length),
        "--merge_last",
        "-b {}".format(os.path.join(output_directory, "intermediate", "scaffolds_fragmented.bed")),
        ">",
        os.path.join(output_directory, "intermediate", "scaffolds_fragmented.fasta"),

        "&&",

        # concoct_coverage_table.py contigs_10K.bed mapping/Sample*.sorted.bam > coverage_table.tsv
        os.environ["concoct_coverage_table.py"],
        os.path.join(output_directory, "intermediate", "scaffolds_fragmented.bed"), 
        " ".join(opts.bam), 
        ">",
        os.path.join(output_directory, "intermediate", "coverage_table.tsv"), 

        "&&",

        # concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output/
        os.environ["concoct"],
        "--composition_file {}".format(os.path.join(output_directory, "intermediate", "scaffolds_fragmented.fasta")),
        "--coverage_file {}".format(os.path.join(output_directory, "intermediate", "coverage_table.tsv")),
        "--length_threshold {}".format(opts.minimum_contig_length),
        "--seed {}".format(opts.random_state),
        "--threads {}".format(opts.n_jobs),
        "-b {}".format(os.path.join(output_directory, "intermediate")),
        opts.concoct_options,

        "&&",
        # merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv

        os.environ["merge_cutup_clustering.py"],
        os.path.join(output_directory, "intermediate", "clustering_gt{}.csv".format(opts.minimum_contig_length)),
        ">",
        os.path.join(output_directory, "intermediate", "clustering_gt{}.merged.csv".format(opts.minimum_contig_length)),


        # extract_fasta_bins.py original_contigs.fa concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins

        "&&",

        os.environ["extract_fasta_bins.py"],
        opts.fasta,
        os.path.join(output_directory, "intermediate", "clustering_gt{}.merged.csv".format(opts.minimum_contig_length)),
        "--output_path {}".format(os.path.join(output_directory, "bins/")),

    ]

    # Remove small bins
    cmd += [ 
r"""

for FP in %s;
    do GENOME_SIZE=$(cat $FP | grep -v "^>" | tr -d "\n" | wc -m);
    echo "[GENOME SIZE] ${FP} = ${GENOME_SIZE}"

    if (( %d > ${GENOME_SIZE} )); then
        echo "[REMOVING] ${FP}"
        rm -f ${FP}
    fi

    done

"""%( 
    os.path.join(output_directory, "bins", "*.fa"),
    opts.minimum_genome_length,
    )
    ]

    return cmd

# SemiBin2
def get_semibin2_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Create dummy scaffolds_to_bins.tsv to overwrite later. 
    cmd = [
         "echo '' > {}".format(os.path.join(output_directory, "scaffolds_to_bins.tsv")),
         "&&",
         "mkdir -p {}".format(os.path.join(output_directory, "bins")),
    ]

    cmd += [
        "&&",
        
        os.environ["SemiBin2"],
        "single_easy_bin",
        "-i {}".format(opts.fasta),
        "-b {}".format(" ".join(opts.bam)) if opts.bam else "",
        "--depth-metabat2 {}".format(opts.coverage) if opts.coverage else "",
        "-o {}".format(os.path.join(output_directory, "intermediate")),
        "-m {}".format(max(opts.minimum_contig_length, 1000)),
        "-p {}".format(opts.n_jobs),
        "--environment {}".format(opts.semibin2_biome) if opts.semibin2_biome else "",
        "--orf-finder {}".format(opts.semibin2_orf_finder),
        # "--prodigal-output-faa {}".format(opts.proteins) if opts.proteins else "",
        "--minfasta-kbs {}".format(opts.minimum_genome_length//1000),
        "--random-seed {}".format(opts.random_state),
        "--verbose",
        "--tmpdir {}".format(directories["tmp"]),
        "--engine {}".format(opts.semibin2_engine),
        "--compression none",
        "--tag-output bin",
        "--sequencing-type {}".format(opts.semibin2_sequencing_type),
        opts.semibin2_options,
        ]
    
    cmd += [
        "&&",
        "mv",
        os.path.join(output_directory, "intermediate", "output_bins", "*.fa"),
        os.path.join(output_directory, "bins"),
    ]


    return cmd

# MetaDecoder
def get_metadecoder_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Create dummy scaffolds_to_bins.tsv to overwrite later. 
    cmd = [
         "echo '' > {}".format(os.path.join(output_directory, "scaffolds_to_bins.tsv")),
         "&&",
         "mkdir -p {}".format(os.path.join(output_directory, "bins")),
    ]

    coverage_file = opts.coverage
    if not opts.coverage:
        coverage_file = os.path.join(output_directory, "intermediate", "coverage.tsv")
        cmd += [

          "&&",
        os.environ["metadecoder"],
        "coverage",
        "-b",
        " ".join(opts.bam),
        "-o {}".format(coverage_file),
        "--threads {}".format(opts.n_jobs),
        opts.metadecoder_coverage_options,
    ]
    
    cmd += [

            "&&",
    
        os.environ["metadecoder"],
        "seed",
        "-f {}".format(opts.fasta),
        "-o {}".format(os.path.join(output_directory, "intermediate", "seed.tsv")),
        "--proteins {}".format(opts.proteins) if opts.proteins else "",
        "--proteins_to_contigs {}".format(opts.proteins_to_contigs) if opts.proteins_to_contigs else "",
        "--hmm_database {}".format(opts.hmm_database) if opts.hmm_database else "",
        "--hmm_marker_field {}".format(opts.hmm_marker_field),
        "--score_type {}".format(opts.score_type),
        "--threshold_method {}".format(opts.threshold_method),
        "--evalue {}".format(opts.evalue),

            "&&",

        os.environ["metadecoder"],
        "cluster",
        "-f {}".format(opts.fasta),
        "-c {}".format(coverage_file),
        "-s {}".format(os.path.join(output_directory, "intermediate", "seed.tsv")),
        "-o {}".format(os.path.join(output_directory, "bins", "bin")),
        "--min_sequence_length {}".format(max(opts.minimum_contig_length, 2000)),
        "--min_cluster_size {}".format(opts.minimum_genome_length),
        "--random_number {}".format(opts.random_state),
        opts.metadecoder_cluster_options,
        
        # Move bins and change extension
"""

for FP in %s;
    do ID_GENOME=$(basename ${FP} .fasta);
    mv $FP %s/${ID_GENOME}.fa
    done

"""%( 
    os.path.join(output_directory, "bins", "*.fasta"),
    os.path.join(output_directory, "bins"),
    )
    ]
    

    return cmd


# MetaCoAG
def get_metacoag_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Create dummy scaffolds_to_bins.tsv to overwrite later. 
    cmd = [
         "echo '' > {}".format(os.path.join(output_directory, "scaffolds_to_bins.tsv")),
         "&&",
        "mkdir -p {}".format(os.path.join(output_directory, "intermediate")),


    ]

    # Coverage for MetaCoAG 
    coverage_file = opts.coverage

    if opts.bam:
        coverage_file = os.path.join(output_directory,"intermediate", "coverage.tsv")
        cmd += [
        "&&",
        os.environ["coverm"],
        "contig",
        "--threads {}".format(opts.n_jobs),
         "--methods",
         "mean",
         "--bam-files",
        " ".join(opts.bam),
        "|",
        "tail",
        "-n",
        "+2",
        ">",
        coverage_file,
        ]
        
    cmd += [
        "&&",
       
        os.environ["metacoag"],
        "--assembler {}".format(opts.metacoag_assembler),
        "--contigs {}".format(opts.fasta), # scaffolds.fasta
        "--abundance {}".format(coverage_file),
        "--graph {}".format(opts.metacoag_graph), 
        "--paths {}".format(opts.metacoag_paths) if opts.metacoag_paths else "", 
        "--proteins {}".format(opts.proteins) if opts.proteins else "", 
        "--proteins_to_contigs {}".format(opts.proteins_to_contigs) if opts.proteins_to_contigs else "", 
        "--output {}".format(os.path.join(output_directory, "intermediate")),
        "--min_bin_size {}".format(opts.minimum_genome_length),
        "--min_length {}".format(opts.minimum_contig_length),
        "--nthreads {}".format(opts.n_jobs),
        "--hmm {}".format(opts.hmm_database) if opts.hmm_database else "",
        "--hmm_marker_field {}".format(opts.hmm_marker_field),
        "--score_type {}".format(opts.score_type),
        "--threshold_method {}".format(opts.threshold_method),
        "--evalue {}".format(opts.evalue),
        opts.metacoag_options,
        
        "&&",
        
        "mv",
        os.path.join(output_directory, "intermediate", "bins"),
        os.path.join(output_directory),
        
        
        ]


    return cmd

# VAMB
# The coverage table needs to match the contigs exactly
def get_vamb_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    # Create dummy scaffolds_to_bins.tsv to overwrite later. 
    cmd = [
         "echo '' > {}".format(os.path.join(output_directory, "scaffolds_to_bins.tsv")),
         "&&",
         "mkdir -p {}".format(os.path.join(output_directory, "bins")),
    ]
    
    # Get contigs for coverage
    cmd += [ 
            "&&",
            # Get contig list
            "cat",
            opts.fasta,
            "|",
            "grep",
            '"^>"',
            "|",
            "cut",
            "-c2-",
            "|",
            "cut",
            "-d",
            '" "',
            "-f1",
            ">",
            os.path.join(output_directory, "intermediate", "contigs.list"),
    ]

    # Coverage for VAMB 
    coverage_file = opts.coverage

    if opts.bam:
        coverage_file = os.path.join(output_directory,"intermediate", "coverage.tsv")
        cmd += [


            "&&",
            
            os.environ["coverm"],
            "contig",
            "--threads {}".format(opts.n_jobs),
            "--methods metabat",
            "--bam-files",
            " ".join(opts.bam),
            "|",
            os.environ["convert_metabat2_coverage.py"],
            "--index_name"
            "contigname",
            "--identifiers",
            os.path.join(output_directory, "intermediate", "contigs.list"),
            ">",
            coverage_file,
        ]
        
    else:
        coverage_file = os.path.join(output_directory,"intermediate", "coverage.tsv")
        cmd += [
            "&&",
            os.environ["subset_table.py"],
            "-i",
            os.path.join(output_directory, "intermediate", "contigs.list"),
            "-t",
            opts.coverage,
            "--index_name",
            "contigname",
            ">",
            coverage_file,
        ]
        
    cmd += [

            "&&",

        os.environ["vamb"],
        "bin",
        "default",
        "--fasta",
        opts.fasta,
        "--abundance_tsv",
        coverage_file,
        "--outdir",
        os.path.join(output_directory, "intermediate", "results"),
        "-p",
        opts.n_jobs,
        "-m",
        opts.minimum_contig_length,
        "--minfasta",
        opts.minimum_genome_length,
        "--seed",
        opts.random_state,
        {"cpu":"","gpu":"--cuda"}[opts.vamb_engine],
        "-o",
        opts.vamb_options,
        
            
        # Move bins and change extension
"""

for FP in %s;
    do ID_GENOME=$(basename ${FP} .fna);
    mv $FP %s/${ID_GENOME}.fa
    done

"""%( 
    os.path.join(output_directory, "intermediate", "results","bins", "*.fna"),
    os.path.join(output_directory, "bins"),
    )

    ]
    return cmd

def get_compile_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # scaffolds_to_bins
    cmd = [
        os.environ["scaffolds_to_bins.py"],
        "-x fa",
        "-i {}".format(os.path.join(output_directory, "bins")),
    ]
    if opts.bin_prefix:
        cmd += ["--bin_prefix {}".format(opts.bin_prefix)]

    cmd += [
        ">",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),

        # "binned.list"
       "&&",

        "cut",
        "-f1",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),
        ">",
        os.path.join(output_directory, "binned.list"),

        # bins.list
        "&&",

        "cut -f2",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),
        "|",
        "sort",
        "-u",
        ">",
        os.path.join(output_directory, "bins.list"),

        # unbinned.list
        "&&",
        "cat",
        opts.fasta,
        "|",
        "grep",
        '"^>"',
        "|",
        "cut -c2-",
        ">",
        os.path.join(output_directory, "unbinned.list"),
    ]

    # unbinned.fasta
    if opts.retain_unbinned_fasta:
        cmd += [ 
        "&&",
        "cat",
        opts.fasta,
        "|",
        os.environ["seqkit"],
        "grep",
        "--pattern-file {}".format(os.path.join(output_directory, "binned.list")),
        "-v",
        ">",
        os.path.join(output_directory, "unbinned.fasta"),
        ]

    # relabel with prefix
    if opts.bin_prefix:
        cmd += [ 
    """

    for FP in %s;
        do ID_GENOME=$(basename ${FP} .fa);
        mv $FP %s/%s${ID_GENOME}.fa
        done

    """%( 
        os.path.join(output_directory, "bins", "*.fa"),
        os.path.join(output_directory, "bins"),
        opts.bin_prefix,
        )
        ]
    else:
        cmd += ["&&"]

    cmd += [ 
        # "genome_statistics.tsv"
        # "&&",
        os.environ["seqkit"],
        "stats",
        "--basename",
        "--all",
        "-T",
        "-j {}".format(opts.n_jobs),
        os.path.join(output_directory, "bins", "*.fa"),

            "|",

        # Remove extension
        """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: x[:-3]); df.to_csv(sys.stdout, sep="\t")'"""

            ">",

        os.path.join(output_directory,  "genome_statistics.tsv"),
    ]


    if opts.remove_bins:
        cmd += [ 
            "&&",
            "rm -rf {}".format(os.path.join(output_directory,"bins")),
        ]


    return cmd

def get_multisplit_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [

        os.environ["partition_multisplit_bins.py"],
        "--scaffolds_to_samples {}".format(opts.multisplit),
        "--binning_directory {}".format(input_filepaths[0]),
        "--output_directory {}".format(output_directory),
        "-x fa",
        "-d {}".format(opts.delimiter),

    ]

    if not opts.remove_bins:
        cmd += [

            "&&",

        os.environ["seqkit"],
        "stats",
        "--basename",
        "--all",
        "-T",
        "-j {}".format(opts.n_jobs),
        os.path.join(output_directory, "*", "bins", "*.fa"),

            "|",

        # Remove extension
        """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: x[:-3]); df.to_csv(sys.stdout, sep="\t")'"""

            ">",

        os.path.join(output_directory,  "genome_statistics.tsv"),

            "&&",

        "for FP in %s; do DIR=$(dirname $FP); %s -i $FP -t %s -o ${DIR}/genome_statistics.tsv; done"%(
            os.path.join(output_directory, "*", "bins.list"),
            os.environ["subset_table.py"],
            os.path.join(output_directory,  "genome_statistics.tsv"),
            ),
    ]

    return cmd

# ============
# Run Pipeline
# ============
def create_pipeline(opts, directories, f_cmds):

    # .................................................................
    # Primordial
    # .................................................................
    # Commands file
    pipeline = ExecutablePipeline(name=__program__, description=opts.algorithm, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])


    if opts.random_state == -1:
        opts.random_state = random.choice(range(1,10000))


    # ==========
    # Binning
    # ==========
    step = 1

    program = "binning_{}".format(opts.algorithm)

    program_label = "{}__{}".format(step, program)

    # Add to directories


    # Info
    description = "Binning via {}".format(opts.algorithm)
    
    # i/o
    output_directory = directories["output"]
    output_filepaths = [
        os.path.join(output_directory, "bins", "*.fa"),
    ]

    # Metabat2
    if opts.algorithm == "metabat2":

        if opts.bam:
            input_filepaths = [
                opts.fasta, 
                *opts.bam,
                ]
        else:
            input_filepaths = [
                opts.fasta, 
                opts.coverage,
                ]            


        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_metabat2_cmd(**params)
        
    # SemiBin2
    elif opts.algorithm == "semibin2":
        if opts.bam:
            input_filepaths = [
                opts.fasta, 
                *opts.bam,
                ]
        else:
            input_filepaths = [
                opts.fasta, 
                opts.coverage,
                ]            

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_semibin2_cmd(**params)
        
    # MetaDecoder
    elif opts.algorithm == "metadecoder":
        if opts.bam:
            input_filepaths = [
                opts.fasta, 
                *opts.bam,
                ]
        else:
            input_filepaths = [
                opts.fasta, 
                opts.coverage,
                ]            


        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_metadecoder_cmd(**params)
        
    # VAMB
    elif opts.algorithm == "vamb":
        if opts.bam:
            input_filepaths = [
                opts.fasta, 
                *opts.bam,
                ]
        else:
            input_filepaths = [
                opts.fasta, 
                opts.coverage,
                ]            


        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_vamb_cmd(**params)
        
    # MetaCoAG
    elif opts.algorithm == "metacoag":
        if opts.bam:
            input_filepaths = [
                opts.fasta, 
                *opts.bam,
                ]
        else:
            input_filepaths = [
                opts.fasta, 
                opts.coverage,
                ]            


        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_metacoag_cmd(**params)

    elif opts.algorithm == "maxbin2":
        # Output directory
        output_filepaths = [
            os.path.join(output_directory, "bins", "*.fa"),
        ]
        if opts.bam:
            input_filepaths = [
                opts.fasta, 
                *opts.bam,
                ]
        else:
            input_filepaths = [
                opts.fasta, 
                opts.coverage,
                ]            


        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_maxbin2_cmd(**params)

    elif opts.algorithm == "concoct":
        # Output directory


        input_filepaths = [
            opts.fasta, 
            *opts.bam,
            ]

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_concoct_cmd(**params)

    # Add step to pipeline
    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
                log_prefix=program_label,
    )


   # =============
    # Compile
    # =============
    step = 2
    program = "compile"

    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories["output"]

    # Info

    description = "Compiling results for output"

    # i/o
    input_filepaths = [
        os.path.join(directories["output"],  "bins","*.fa"),
    ]


    output_filenames =  [
        "scaffolds_to_bins.tsv", 
        "binned.list",
        "unbinned.list",
        "bins.list",
        "genome_statistics.tsv",
    ]
    if not opts.remove_bins:
        output_filenames += ["bins/*.fa"]
    if opts.retain_unbinned_fasta:
        output_filenames += ["unbinned.fasta"]

    output_filepaths = list(map(lambda fn:os.path.join(directories["output"], fn), output_filenames))

    
    params = {
    "input_filepaths":input_filepaths,
    "output_filepaths":output_filepaths,
    "output_directory":output_directory,
    "opts":opts,
    "directories":directories,
    }

    cmd = get_compile_cmd(**params)
    pipeline.add_step(
            id=program_label,
            description = description,
            step=step,
            cmd=cmd,
            input_filepaths = input_filepaths,
            output_filepaths = output_filepaths,
            validate_inputs=True,
            validate_outputs=True,
            log_prefix=program_label,

    )

    if opts.multisplit:
        step = 3
        program = "multisplit"

        program_label = "{}__{}".format(step, program)

        # Add to directories
        output_directory = os.path.join(directories["output"], "multisplit")

        # Info

        description = "Partition bins by sample"

        # i/o
        input_filepaths = [
            directories["output"], 
        ]


        output_filenames =  [
            "*/scaffolds_to_bins.tsv", 
            "*/binned.list",
            "*/bins.list",
        ]
        if not opts.remove_bins:
            output_filenames += [ 
            "*/genome_statistics.tsv",
            "*/bins/*.fa",
            ]


        output_filepaths = list(map(lambda fn:os.path.join(output_directory, fn), output_filenames))

        
        params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        }

        cmd = get_multisplit_cmd(**params)
        pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
                log_prefix=program_label,

        )

    return pipeline



# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = {
        "scaffolds_to_bins.py",
        "partition_multisplit_bins.py",
        "subset_table.py",
        "convert_metabat2_coverage.py",
        "cut_table_by_column_index.py",
        "partition_bins.py",
    }

    required_executables={
                "seqkit",
     } 
               
    # Metabat2
    if opts.algorithm == "metabat2":
        required_executables |= {
            "coverm",
            "metabat2",
            }
    
    # SemiBin2
    if opts.algorithm == "semibin2":
            required_executables |= {
                "SemiBin2",
                }

    # MetaDecoder
    if opts.algorithm == "metadecoder":
        required_executables |= {
            "metadecoder",
            "pyhmmsearch",
            "pyrodigal",
            }
        
    # MetaCoAG
    if opts.algorithm == "metacoag":
        required_executables |= {
            "coverm",
            "metacoag",
            "pyhmmsearch",
            "pyrodigal",
            }

    # MaxBin2
    if opts.algorithm == "maxbin2":
        required_executables |= {
            "coverm",
            "run_MaxBin.pl",
            }
        
    # CONCOCT
    if opts.algorithm == "concoct":
         required_executables  |= {
                "concoct",
                "cut_up_fasta.py", 
                "concoct_coverage_table.py",
                "merge_cutup_clustering.py", 
                "extract_fasta_bins.py", 
            } 
            
    if opts.algorithm == "vamb":
        required_executables |= {
            "coverm",
            "vamb",
            }

    required_executables |= accessory_scripts

    if opts.path_config == "CONDA_PREFIX":
        executables = dict()
        for name in sorted(required_executables):
            if name not in accessory_scripts:
                executables[name] = os.path.join(os.environ["CONDA_PREFIX"], "bin", name)
    else:
        if opts.path_config is None:
            opts.path_config = os.path.join(opts.script_directory, "veba_config.tsv")
        opts.path_config = format_path(opts.path_config)
        assert os.path.exists(opts.path_config), "config file does not exist.  Have you created one in the following directory?\n{}\nIf not, either create one, check this filepath:{}, or give the path to a proper config file using --path_config".format(opts.script_directory, opts.path_config)
        assert os.stat(opts.path_config).st_size > 1, "config file seems to be empty.  Please add 'name' and 'executable' columns for the following program names: {}".format(required_executables)
        df_config = pd.read_csv(opts.path_config, sep="\t")
        assert {"name", "executable"} <= set(df_config.columns), "config must have `name` and `executable` columns.  Please adjust file: {}".format(opts.path_config)
        df_config = df_config.loc[:,["name", "executable"]].dropna(how="any", axis=0).applymap(str)
        # Get executable paths
        executables = OrderedDict(zip(df_config["name"], df_config["executable"]))
        assert required_executables <= set(list(executables.keys())), "config must have the required executables for this run.  Please adjust file: {}\nIn particular, add info for the following: {}".format(opts.path_config, required_executables - set(list(executables.keys())))

    # Display

    for name in sorted(accessory_scripts):
        executables[name] = "'{}'".format(os.path.join(opts.script_directory, name)) # Can handle spaces in path


    print(format_header( "Adding executables to path from the following source: {}".format(opts.path_config), "-"), file=sys.stdout)
    for name, executable in executables.items():
        if name in required_executables:
            print(name, executable, sep = " --> ", file=sys.stdout)
            os.environ[name] = executable.strip()
    print("", file=sys.stdout)


# Configure parameters
def configure_parameters(opts, directories):
    # assert opts.reference_assembly is not None, "Must include --reference_assembly"
    assert not all([bool(opts.bam), bool(opts.coverage)]), "Cannot have both --bam and --coverage"
    if opts.algorithm in {"metabat2", "maxbin2", "semibin2", "metacoag", "metadecoder"}:
        
        assert  any([bool(opts.bam), bool(opts.coverage)]), f"Must have either --bam or --coverage for --algorithm  {opts.algorithm}"
    if opts.algorithm in {"concoct"}:
        assert opts.bam is not None, f"Must provide --bam for --algorithm  {opts.algorithm}"
        
    if opts.semibin2_biome == "NONE":
        opts.semibin2_biome = None

    if opts.bin_prefix == "DEFAULT":
        if opts.algorithm == "maxbin2":
            opts.bin_prefix = "{}-{}__".format(opts.algorithm.upper(), opts.maxbin2_markerset)

        else:
            opts.bin_prefix = "{}__".format(opts.algorithm.upper())

    if opts.bin_prefix == "NONE":
        opts.bin_prefix = None
        
    if opts.algorithm == "metacoag":
        # assert opts.metacoag_assembler is not None, "Must provide --assembler if --algorithm = metacoag"
        assert opts.metacoag_graph is not None, "Must provide --graph if --algorithm = metacoag"

        if opts.metacoag_assembler in {"spades", "flye"}:
            assert opts.metacoag_paths is not None, "Must provide --paths if --algorithm = metacoag and --assembler = spades or flye"
    

    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(argv=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -f <scaffolds.fasta> -b <mapped.sorted.bam>|-c <coverage.tsv> -a <algorithm> -o <output_directory>".format(__program__)

    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-f","--fasta", type=str, required=True, help = "path/to/scaffolds.fasta")
    parser_io.add_argument("-b","--bam", type=str, required=False, nargs="+", help = "path/to/mapped.sorted.bam files separated by spaces. Can be used with any binning algorithm but cannot be used with --coverage")
    parser_io.add_argument("-c","--coverage", type=str, required=False, help = "path/to/coverage.tsv. Can be used with --algorithm metabat2 maxbin2 but not available for --algorithm concoct and cannot be used with --bam")
    parser_io.add_argument("-o","--output_directory", type=str,  help = "path/to/binning_output [Default: binning_output/[algorithm_name]]")
    parser_io.add_argument("--proteins", type=str, help = "path/to/proteins.fasta used to bypass gene prediction.  Must be Prodigal format where [id_contig]_[gene_number]")
    parser_io.add_argument("--proteins_to_contigs", type=str, help = "path/to/proteins_to_contigs.tsv Tab-delimited file mapping proteins to contigs [id_protein]<tab>[id_contig].  If --proteins and provided without --proteins_to_contigs then id_protein formatting is assumed to be [id_contig]_[gene_number]")

    parser_multisplit = parser.add_argument_group('Multisplit arguments')
    parser_multisplit.add_argument("-M","--multisplit", type=str,  help = "path/to/scaffolds_to_samples.tsv. Use this to perform multisplit binning. Expected table input is the following format: [id_scaffold]<tab>[id_sample], No header. [Optional]")
    parser_multisplit.add_argument("-d", "--delimiter", type=str,  default="__", help="Delimiter between [id_sample]<delimiter>[id_bin].  Only used with multisplit. [Default: __ which is [id_sample]__[id_bin]")

    # Multi-split will make the output/sample_x/ and then all of the files that are normally in output but per sample
    
    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')

    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Binning
    parser_binning = parser.add_argument_group('Binning arguments')
    parser_binning.add_argument("-a", "--algorithm", type=str, default="metabat2", choices={"concoct", "metabat2", "maxbin2", "semibin2", "metadecoder", "metacoag", "vamb"}, help="Binning algorithm [Default: metabat2] ")
    parser_binning.add_argument("-m", "--minimum_contig_length", type=int, default=1500, help="Minimum contig length.  [Default: 1500] ")
    parser_binning.add_argument("-s", "--minimum_genome_length", type=int, default=200000, help="Minimum genome length.  [Default: 200000] ")
    parser_binning.add_argument("-P","--bin_prefix", type=str,  default="DEFAULT", help = "Prefix for bin names. Special strings include: 1) --bin_prefix NONE which does not include a bin prefix; and 2) --bin_prefix DEFAULT then prefix is [ALGORITHM_UPPERCASE]__")
    parser_binning.add_argument("-B", "--remove_bins", action="store_true", help="Remove fasta files for bins")
    parser_binning.add_argument("-U", "--retain_unbinned_fasta", action="store_true", help="Retain unbinned fasta")
    parser_binning.add_argument("-R", "--remove_intermediate_files", action="store_true", help="Remove intermediate files. Warning: Cannot use checkpoints after this.")

    parser_hmms = parser.add_argument_group('HMM arguments for Metadecoder and MetaCoag')
    parser_hmms.add_argument("--hmm_database", type=str, help="Path to the HMM reference database.  Only include if you want to use a custom HMM and not recommended unless you know what you're doing.")
    parser_hmms.add_argument("--hmm_marker_field", default="accession", type=str, choices={"accession", "name"}, help="HMM reference type (accession, name) [Default: accession]")
    parser_hmms.add_argument("--score_type",  default="full", type=str, choices = {"full", "domain"}, help="{full, domain} [Default: full]")
    parser_hmms.add_argument("--threshold_method", type=str, default="trusted",choices={"gathering", "noise", "trusted", "e"},  help="Cutoff threshold method [Default:  trusted]")
    parser_hmms.add_argument("--evalue", type=float, default=10.0,  help = "E-value threshold [Default: 10.0]")
    
    # Metabat2
    parser_metabat2 = parser.add_argument_group('Metabat2 arguments')
    parser_metabat2.add_argument("--metabat2_options", type=str, default="", help="MetaBat2 | More options (e.g. --arg 1 ) [Default: ''] | https://bitbucket.org/berkeleylab/metabat/src/master/")

    # SemiBin2
    parser_semibin2 = parser.add_argument_group('SemiBin2 arguments')
    parser_semibin2.add_argument("--semibin2_orf_finder", type=str, choices={'fast-naive', 'prodigal', 'fraggenescan'}, default="fast-naive", help="SemiBin2 | ORF finder used to estimate the number of bins  [Default: fast-naive]")
    parser_semibin2.add_argument("--semibin2_biome", type=str, choices={'ocean', 'wastewater', 'global', 'pig_gut', 'human_oral', 'cat_gut', 'soil', 'chicken_caecum', 'human_gut', 'built_environment', 'dog_gut', 'mouse_gut', 'NONE'}, default="global", help="SemiBin2 | Biome/environment for the built-in model.  Use 'NONE' to implement Semi-Supervised training (takes longer with more compute) [Default: global]")
    parser_semibin2.add_argument("--semibin2_engine", type=str, choices={'auto', 'cpu', 'gpu'}, default="auto", help="SemiBin2 | Device used to train the model [Default: auto]")
    parser_semibin2.add_argument("--semibin2_sequencing_type", type=str, choices={'short_read', 'long_read'}, default="short_read", help="SemiBin2 | Sequencing type [Default: short_read]")
    parser_semibin2.add_argument("--semibin2_options", type=str, default="", help="SemiBin2 | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/BigDataBiology/SemiBin")

    # MetaDecoder
    parser_metadecoder = parser.add_argument_group('MetaDecoder arguments')
    parser_metadecoder.add_argument("--metadecoder_coverage_options", type=str, default="", help="MetaDecoder | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/jolespin/metadecoder-nal")
    parser_metadecoder.add_argument("--metadecoder_cluster_options", type=str, default="", help="MetaDecoder | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/jolespin/metadecoder-nal")

    # MetaCoAG
    parser_metacoag = parser.add_argument_group('MetaCoAG arguments')
    parser_metacoag.add_argument("--metacoag_assembler", type=str, choices={"auto", "spades", "megahit", "flye"}, default="auto", help="MeteaCoAG | Assembler used during assembly [Required if MetaCoAG is used]")
    parser_metacoag.add_argument("--metacoag_graph", default="auto", type=str, help="MetaCoAG | de Bruijn graph from SPAdes, MEGAHIT, or metaFlye [Required if MetaCoAG is used, if `auto` then assembly graphs will be looked]")
    parser_metacoag.add_argument("--metacoag_paths", default="auto", type=str, help="MetaCoAG | de Bruijn graph paths from SPAdes or metaFlye [Required if MetaCoAG is used with SPAdes or metaFlye]")
    parser_metacoag.add_argument("--metacoag_options", type=str, default="", help="MetaCoAG | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/jolespin/metacoag-nal")

    # VAMB
    parser_vamb = parser.add_argument_group('VAMB arguments')
    parser_vamb.add_argument("--vamb_engine", type=str, choices={'cpu', 'gpu'}, default="cpu", help="VAMB | Device used to train & cluster [Default: cpu]")
    parser_vamb.add_argument("--vamb_options", type=str, default="", help="VAMB | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/RasmussenLab/vamb")

    # MaxBin2
    parser_maxbin2 = parser.add_argument_group('MaxBin2 arguments')
    parser_maxbin2.add_argument("--maxbin2_markerset", type=int, default=107, help="MaxBin2 marker set: {107, 40} [Default: 107] ")
    parser_maxbin2.add_argument("--maxbin2_options", type=str, default="", help="MaxBin2 | More options (e.g. --arg 1 ) [Default: ''] | https://sourceforge.net/projects/maxbin2/")

    # CONCOCT
    parser_concoct = parser.add_argument_group('CONCOCT arguments')
    parser_concoct.add_argument("--concoct_fragment_length", type=int, default=10000, help="CONCOCT | Fragment length [Default: 10000] ")
    parser_concoct.add_argument("--concoct_overlap_length", type=int, default=0, help="CONCOCT | Fragment overlap length [Default: 0] ")
    parser_concoct.add_argument("--concoct_options", type=str, default="", help="CONCOCT | More options (e.g. --arg 1 ) [Default: '']")

    # Options
    opts = parser.parse_args(argv)

    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    # Threads
    if opts.n_jobs == -1:
        from multiprocessing import cpu_count 
        opts.n_jobs = cpu_count()
    assert opts.n_jobs >= 1, "--n_jobs must be ≥ 1.  To select all available threads, use -1."

    # Directories
    if opts.multisplit:
        markerbin_based_algorithms = {"maxbin2", "semibin2", "metadecoder", "metacoag"}
        assert opts.algorithm not in markerbin_based_algorithms, f"{opts.algorithm} is marker-based and cannot be used for multisplit binning"

    if not opts.output_directory:
        if opts.algorithm == "maxbin2":
            opts.output_directory = "binning_output/{}-{}".format(opts.algorithm, opts.maxbin2_markerset)
        else:
            opts.output_directory = "binning_output/{}".format(opts.algorithm)
            
    

    directories = dict()
    directories["output"] = create_directory(opts.output_directory)
    directories["intermediate"] = create_directory(os.path.join(directories["output"], "intermediate"))
    directories["log"] = create_directory(os.path.join(directories["output"],"intermediate", "log"))
    directories["tmp"] = create_directory(os.path.join(directories["output"],"intermediate", "tmp"))
    directories["checkpoints"] = create_directory(os.path.join(directories["output"],"intermediate", "checkpoints"))
    # os.environ["TMPDIR"] = directories["tmp"]

    # Info
    print(format_header(__program__, "="), file=sys.stdout)
    print(format_header("Configuration:", "-"), file=sys.stdout)
    print("Python version:", sys.version.replace("\n"," "), file=sys.stdout)
    print("Python path:", sys.executable, file=sys.stdout) #sys.path[2]
    print("Script version:", __version__, file=sys.stdout)
    print("Moment:", get_timestamp(), file=sys.stdout)
    print("Directory:", os.getcwd(), file=sys.stdout)
    print("Commands:", list(filter(bool,sys.argv)),  sep="\n", file=sys.stdout)
    configure_parameters(opts, directories)
    sys.stdout.flush()


    # Run pipeline
    with open(os.path.join(directories["output"], "commands.sh"), "w") as f_cmds:
        pipeline = create_pipeline(
                    opts=opts,
                    directories=directories,
                    f_cmds=f_cmds,
        )
        pipeline.compile()
        pipeline.execute(restart_from_checkpoint=opts.restart_from_checkpoint)

    if opts.remove_intermediate_files:
        shutil.rmtree(directories["intermediate"], ignore_errors=True)

if __name__ == "__main__":
    main(sys.argv[1:])





# Extra:
