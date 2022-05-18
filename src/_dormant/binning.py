#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
from collections import OrderedDict, defaultdict

import pandas as pd
import numpy as np

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.08.12"

# .............................................................................
# Primordial
# .............................................................................


# Assembly
def get_coverage_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]

    # Command
    cmd = [
    "(",
    # Coverage for MetaBat2 & VAMB
    os.environ["jgi_summarize_bam_contig_depths"],
    "--outputDepth {}".format(output_filepaths[0]),
    input_filepaths[0],
    ")",
        "&&",
    "(",
    # Coverage for MaxBin2
    "cut -f1,4",
    output_filepaths[0],
    "|",
    "tail -n +2",
    ">",
    output_filepaths[1],
    ")",
        "&&",
    "(",
    "cat {}".format(output_filepaths[0]),
    "|",
    "python -c 'import sys, pandas as pd; pd.read_csv(sys.stdin, sep=\"\t\", index_col=0).query(\"contigLen >= {}\").to_csv(sys.stdout, sep=\"\t\")'".format(opts.minimum_contig_length),
    ">",
    output_filepaths[2],
    ")",
    ]
    return cmd

        # Metabat2
    	# jgi_summarize_bam_contig_depths --outputDepth ${out}/work_files/metabat_depth.txt ${out}/work_files/*.bam
        # Maxbin2
        # jgi_summarize_bam_contig_depths --outputDepth ${out}/work_files/mb2_master_depth.txt --noIntraDepthVariance ${out}/work_files/*.bam
        #BinSanity
        # Binsanity-profile -i fasta_file -s {sam,bam}_file -c output_file
	# metabat2 -i ${out}/work_files/assembly.fa -a ${out}/work_files/metabat_depth.txt\
	#  -o ${out}/metabat2_bins/bin -m $metabat_len -t $threads --unbinned



# Metabat2
def get_metabat2_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, prefix):
    cmd = [
        "(",
        os.environ["metabat2"],
        "-i {}".format(input_filepaths[0]), # scaffolds.fasta
        "-o {}".format(os.path.join(output_directory, "bin")),
        "-a {}".format(input_filepaths[1]),
        "-m {}".format(max(opts.minimum_contig_length, 1500)), # mininimum contig length must be >= 1500 nts for MetaBat2
        "-t {}".format(opts.n_jobs),
        # "--unbinned", # Don't forget about these...
        "--seed {}".format(opts.random_state),
        "--verbose",
        opts.metabat2_options,
        "&&",
        os.environ["scaffolds_to_bins.py"],
        "-x fa",
        "-i {}".format(output_directory),
        # "--sep '\t'",
        "--bin_prefix {}__METABAT2__{}-".format(opts.name, prefix),
        ">",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),
        ")",
        "&&",
        "cut",
        "-f1",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),
        ">",
        os.path.join(output_directory, "binned.list"),
        "&&",
        "(",
        "rm -f {}".format(os.path.join(output_directory, "*.fa")), # Remove because the info is in the scaffolds to bins file
        ")",
        # "&&",
        # "(mkdir -p {})".format(os.path.join(output_directory, "unbinned")),
        # "&&",
        # "(mv {} {})".format(os.path.join(output_directory, "*unbinned*"),os.path.join(output_directory, "unbinned")),
    ]
    return cmd

# Maxbin2
def get_maxbin2_107_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [
        "(",
        os.environ["run_MaxBin.pl"],
        "-contig {}".format(input_filepaths[0]), # scaffolds.fasta
        "-out {}".format(os.path.join(output_directory, "bin")),
        "-abund {}".format(input_filepaths[1]),
        "-min_contig_length {}".format(opts.minimum_contig_length),
        "-markerset {}".format(107),
        "-thread {}".format(opts.n_jobs),
        "-verbose",
        opts.maxbin2_options,
        "&&",
        os.environ["scaffolds_to_bins.py"],
        "-x fasta",
        "-i {}".format(output_directory),
        # "--sep '\t'",
        "--bin_prefix {}__MAXBIN2-107__P-".format(opts.name),
        ">",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),
        ")",
        "&&",
        "cut",
        "-f1",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),
        ">",
        os.path.join(output_directory, "binned.list"),
        "&&",
        "(",
        "rm -f {}".format(os.path.join(output_directory, "*.fasta")), # Remove because the info is in the scaffolds to bins file
        "&&",
        "rm -f {}".format(os.path.join(output_directory, "bin.noclass")), 
        "&&",
        "rm -f {}".format(os.path.join(output_directory, "bin.tooshort")), 
        ")",
    ]
    return cmd

def get_maxbin2_40_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [
        "(",
        os.environ["run_MaxBin.pl"],
        "-contig {}".format(input_filepaths[0]), # scaffolds.fasta
        "-out {}".format(os.path.join(output_directory, "bin")),
        "-abund {}".format(input_filepaths[1]),
        "-min_contig_length {}".format(opts.minimum_contig_length),
        "-markerset {}".format(40),
        "-thread {}".format(opts.n_jobs),
        "-verbose",
        opts.maxbin2_options,
        "&&",
        os.environ["scaffolds_to_bins.py"],
        "-x fasta",
        "-i {}".format(output_directory),
        # "--sep '\t'",
        "--bin_prefix {}__MAXBIN2-40__P-".format(opts.name),
        ">",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),
        ")",
        "&&",
        "cut",
        "-f1",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),
        ">",
        os.path.join(output_directory, "binned.list"),
        "&&",
        "(",
        "rm -f {}".format(os.path.join(output_directory, "*.fasta")), # Remove because the info is in the scaffolds to bins file
        "&&",
        "rm -f {}".format(os.path.join(output_directory, "bin.noclass")), 
        "&&",
        "rm -f {}".format(os.path.join(output_directory, "bin.tooshort")), 
        ")",
    ]
    return cmd

def get_concoct_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [
        # cut_up_fasta.py original_contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
        "(",
        os.environ["cut_up_fasta.py"],
        input_filepaths[0], 
        "-c {}".format(opts.concoct_fragment_length),
        "-o {}".format(opts.concoct_overlap_length),
        "--merge_last",
        "-b {}".format(os.path.join(output_directory, "scaffolds_fragmented.bed")),
        ">",
        os.path.join(output_directory, "scaffolds_fragmented.fasta"),
        ")",
        "&&",
        "(",
        # concoct_coverage_table.py contigs_10K.bed mapping/Sample*.sorted.bam > coverage_table.tsv
        os.environ["concoct_coverage_table.py"],
        os.path.join(output_directory, "scaffolds_fragmented.bed"), 
        input_filepaths[1], 
        ">",
        os.path.join(output_directory, "coverage_table.tsv"), 
        ")",
        "&&",
        # concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output/
        "(",
        os.environ["concoct"],
        "--composition_file {}".format(os.path.join(output_directory, "scaffolds_fragmented.fasta")),
        "--coverage_file {}".format(os.path.join(output_directory, "coverage_table.tsv")),
        "--length_threshold {}".format(opts.minimum_contig_length),
        "--seed {}".format(opts.random_state),
        "--threads {}".format(opts.n_jobs),
        "-b {}".format(output_directory),
        opts.concoct_options,
        ")",
        "&&",
        # merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv
        "(",
        os.environ["merge_cutup_clustering.py"],
        os.path.join(output_directory, "clustering_gt{}.csv".format(opts.minimum_contig_length)),
        ">",
        os.path.join(output_directory, "clustering_gt{}.merged.csv".format(opts.minimum_contig_length)),
        ")",
        "&&",
        # extract_fasta_bins.py original_contigs.fa concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins
        "(",
        "mkdir -p {}".format(os.path.join(output_directory, "bins/")),
        "&&",
        os.environ["extract_fasta_bins.py"],
        input_filepaths[0],
        os.path.join(output_directory, "clustering_gt{}.merged.csv".format(opts.minimum_contig_length)),
        "--output_path {}".format(os.path.join(output_directory, "bins/")),
        ")",
        "&&",
        "(",
        os.environ["scaffolds_to_bins.py"],
        "-x fa",
        "-i {}".format(os.path.join(output_directory, "bins")),
        # "--sep '\t'",
        "--bin_prefix {}__CONCOCT__P-".format(opts.name),
        ">",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),
        ")",
        "&&",
        "cut",
        "-f1",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),
        ">",
        os.path.join(output_directory, "binned.list"),
        "&&",
        "(",
        "rm -rf {}".format(os.path.join(output_directory, "bins")), # Remove because the info is in the scaffolds to bins file
        "&&",
        "rm -f {}".format(os.path.join(output_directory, "scaffolds_fragmented.fasta")), 
        "&&",
        "rm -f {}".format(os.path.join(output_directory, "*.csv")), 
        "&&",
        "rm -f {}".format(os.path.join(output_directory, "*.bed")), 
        "&&",
        "rm -f {}".format(os.path.join(output_directory, "coverage_table.tsv")), 
        ")",
    ]
    return cmd

# # VAMB
# def get_vamb_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
#     os.environ["TMPDIR"] = directories["tmp"]
#     # Command
#     cmd = [
#     "(echo $CONDA_PREFIX && echo $PATH)",
#     "&&",
#     # VAMB
#     "(",
#     "rm -r {}".format(output_directory), # There can't be an existing directory for some reason
#     "&&",
#     os.environ["vamb"],
#     "--fasta {}".format(input_filepaths[0]),
#     "--jgi {}".format(input_filepaths[1]),
#     "-m {}".format(opts.minimum_contig_length),
#     # "-p {}".format(opts.n_jobs),
#     "--outdir {}".format(output_directory),
#     opts.vamb_options,
#         "&&",
#     os.environ["scaffolds_to_bins.py"],
#     "-x fna",
#     "-i {}".format(output_directory),
#     # "--sep '\t'",
#     "--bin_prefix VAMB__",
#     ">",
#     os.path.join(output_directory, "scaffolds_to_bins.tsv"),
#     ")",
#     ]
#     return cmd

# DAS_Tool
def get_dastool_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [
        # Get non-empty scaffolds to bins
        "S2B=$({} -i {},{},{},{} -n metabat2,maxbin2-107,maxbin2-40,concoct)".format(
            os.environ["check_scaffolds_to_bins.py"],
            input_filepaths[0],
            input_filepaths[1],
            input_filepaths[2],
            input_filepaths[3],
           ),
           "&&",
           'IFS=" " read -r -a S2B_ARRAY <<< "$S2B"',
           "&&",

        # "echo ${S2B_ARRAY[0]} ${S2B_ARRAY[1]}",
        "(",
        os.environ["DAS_Tool"],
        "--bins ${S2B_ARRAY[0]}", # scaffolds.fasta
        "--contigs {}".format(input_filepaths[4]),
        "--outputbasename {}".format(os.path.join(output_directory, "_")),
        "--labels ${S2B_ARRAY[1]}",
        "--search_engine {}".format(opts.dastool_searchengine),
        "--write_bins 1",
        "--threads {}".format(opts.n_jobs),
        opts.dastool_options,
        ")",
        "&&",
        "cut",
        "-f1",
        os.path.join(output_directory, "__DASTool_scaffolds2bin.txt"),
        ">",
        os.path.join(output_directory, "binned.list"),
    ]
    return cmd

# Prodigal
def get_prodigal_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [
        "(",
        "cat",
        os.path.join(input_filepaths[0], "*.fa"),
        "|",
        os.environ["prodigal"],
        "-p meta",
        "-g {}".format(opts.prodigal_genetic_code),
        "-f gff",
        "-d {}".format(os.path.join(output_directory, "gene_models.ffn")),
        "-a {}".format(os.path.join(output_directory, "gene_models.faa")),
        "|",
        os.environ["append_geneid_to_prodigal_gff.py"],
        "-a gene_id",
        ">",
        os.path.join(output_directory, "gene_models.gff"),
        ")",
        "&&",
        "(",
        os.environ["partition_gene_models.py"],
        "-i {}".format(input_filepaths[1]),
        "-g {}".format(os.path.join(output_directory, "gene_models.gff")),
        "-d {}".format(os.path.join(output_directory, "gene_models.ffn")),
        "-a {}".format(os.path.join(output_directory, "gene_models.faa")),
        "-o {}".format(output_directory),
        ")",
        "&&",
        "(",
        "rm -f {}".format(os.path.join(output_directory, "gene_models.gff")),
        "&&",
        "rm -f {}".format(os.path.join(output_directory, "gene_models.ffn")),
        "&&",
        "rm -f {}".format(os.path.join(output_directory, "gene_models.faa")),
        ")",

    ]
    return cmd

# CheckM
def get_checkm_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # checkm lineage_wf --tab_table -f checkm_output/${ID}/output.tab --pplacer_threads ${N_JOBS} -t ${N_JOBS} -x fa -r ${BINS} ${OUT_DIR}
    cmd = [
        "(",
        os.environ["checkm"],
        "lineage_wf",
        {"full":"", "reduced":"-r"}[opts.checkm_tree],
        "-g", 
        "--tab_table",
        "-f {}".format(output_filepaths[0]),
        "--pplacer_threads {}".format(opts.n_jobs),
        "-t {}".format(opts.n_jobs),
        "-x faa",
        "--tmpdir {}".format(directories["tmp"]),
        input_filepaths[0],
        output_directory,
        ")",
        "&&",
        "(",
        os.environ["filter_checkm_results.py"],
        "-i {}".format(output_filepaths[0]),
        "-b {}".format(input_filepaths[1]),
        "-o {}".format(os.path.join(output_directory, "filtered")),
        "-f {}".format(opts.fasta),
        "-m {}".format(opts.minimum_contig_length),
        "--completeness {}".format(opts.checkm_completeness),
        "--contamination {}".format(opts.checkm_contamination),
        {True:"--strain_heterogeneity {}".format(opts.checkm_strain_heterogeneity),False:""}[bool(opts.checkm_strain_heterogeneity)],
        "-x fa",
        ")",
    ]
    return cmd

# # VirSorter2
# def get_virsorter2_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
#     # cat [] | seqkit seq -M [] | seqkit grep --pattern-file .list > unbinned.fasta
#     cmd = [
#         "(",
#         "cat",
#         opts.fasta,
#         "|",
#         os.environ["seqkit"],
#         "seq",
#         "-m {}".format(opts.minimum_contig_length),
#         "|",
#         os.environ["seqkit"],
#         "grep",
#         "--pattern-file {}".format(input_filepaths[0]),
#         ">",
#         os.path.join(output_directory, "unbinned.fasta"),
#         ")",
#         "&&"
#         "(",
#         os.environ["virsorter"],
#         "run",
#         "--tmpdir {}".format(directories["tmp"]),
#         "--rm-tmpdir",
#         # "-d {}".format(opts.virsorter2_database),
#         {True:"", False:"--provirus-off"}[bool(opts.include_provirus)],
#         "--use-conda-off",
#         "-j {}".format(opts.n_jobs),
#         "-w {}".format(output_directory),
#         "-i {}".format(os.path.join(output_directory, "unbinned.fasta")),
#         "--min-length {}".format(opts.minimum_contig_length),
#         ")",
#         "&&",
#         "rm {}".format(os.path.join(output_directory, "unbinned.fasta")),
#     ]
#     return cmd

# def get_viralverify_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
#     cmd = [
#         "(",
#         "cat",
#         opts.fasta,
#         "|",
#         os.environ["seqkit"],
#         "seq",
#         "-m {}".format(opts.minimum_contig_length),
#         "|",
#         os.environ["seqkit"],
#         "grep",
#         "--pattern-file {}".format(input_filepaths[0]),
#         ">",
#         os.path.join(directories["tmp"], "unbinned.fasta"),
#         ")",
#         "&&",
#         "rm -rf {}".format(output_directory),
#         "&&",
#         "(",
#         os.environ["viral_verify"],
#         "-i {}".format(os.path.join(directories["tmp"], "unbinned.fasta")),
#         "-o {}".format(output_directory),
#         "-H {}".format(opts.pfam_database),
#         "-t {}".format(opts.n_jobs),
#         opts.viralverify_options, 
#         ")",
#         "&&",
#         "rm {}".format(os.path.join(directories["tmp"], "unbinned.fasta")),
#     ]
#     return cmd

def get_checkv_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # checkv end_to_end ${FASTA} ${OUT_DIR} -t ${N_JOBS} --restart

    cmd = [
        "(",
        "cat",
        opts.fasta,
        "|",
        os.environ["seqkit"],
        "seq",
        "-m {}".format(opts.minimum_contig_length),
        "|",
        os.environ["seqkit"],
        "grep",
        "--pattern-file {}".format(input_filepaths[0]),
        ">",
        os.path.join(directories["tmp"], "unbinned_contig_for_checkv.fasta"),
        ")",
        "&&",
        "rm -rf {}".format(output_directory),
        "&&",
        "(",
        os.environ["checkv"],
        "end_to_end",
        os.path.join(directories["tmp"], "unbinned_contig_for_checkv.fasta"),
        output_directory,
        "-t {}".format(opts.n_jobs),
        "-d {}".format(opts.checkv_database),
        "--restart",
        opts.checkv_options, 
        ")",
        "&&",
        "(",
        os.environ["filter_checkv_results.py"],
        "-i {}".format(os.path.join(output_directory, "quality_summary.tsv")),
        "-f {}".format(os.path.join(directories["tmp"], "unbinned_contig_for_checkv.fasta")),
        "-o {}".format(os.path.join(output_directory, "filtered")),
        "-m {}".format(opts.minimum_contig_length),
        "--unbinned",
        "-p {}__Virus.".format(opts.name),
        "--multiplier_viral_to_host_genes {}".format(opts.multiplier_viral_to_host_genes),
        "--completeness {}".format(opts.checkv_completeness),
        "--checkv_quality {}".format(opts.checkv_quality),
        "--miuvig_quality {}".format(opts.miuvig_quality),
        ")",
        "&&",
        "rm {}".format(os.path.join(directories["tmp"], "unbinned_contig_for_checkv.fasta")),
    ]
    return cmd


# def get_metaeuk_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
#     # checkv end_to_end ${FASTA} ${OUT_DIR} -t ${N_JOBS} --restart

#     cmd = [
#         "(",
#         "cat",
#         opts.fasta,
#         "|",
#         os.environ["seqkit"],
#         "seq",
#         "-m {}".format(opts.minimum_contig_length),
#         "|",
#         os.environ["seqkit"],
#         "grep",
#         "--pattern-file {}".format(input_filepaths[0]),
#         ">",
#         os.path.join(directories["tmp"], "unbinned_contig_for_checkv.fasta"),
#         ")",
#         "&&",
#         "rm -rf {}".format(output_directory),
#         "&&",
#         "(",
#         os.environ["checkv"],
#         "end_to_end",
#         os.path.join(directories["tmp"], "unbinned_contig_for_checkv.fasta"),
#         output_directory,
#         "-t {}".format(opts.n_jobs),
#         "-d {}".format(opts.checkv_database),
#         "--restart",
#         opts.checkv_options, 
#         ")",
#         "&&",
#         "(",
#         os.environ["filter_checkv_results.py"],
#         "-i {}".format(os.path.join(output_directory, "quality_summary.tsv")),
#         "-f {}".format(os.path.join(directories["tmp"], "unbinned_contig_for_checkv.fasta")),
#         "-o {}".format(os.path.join(output_directory, "filtered")),
#         "-m {}".format(opts.minimum_contig_length),
#         "--unbinned",
#         "-p {}__Virus.".format(opts.name),
#         "--multiplier_viral_to_host_genes {}".format(opts.multiplier_viral_to_host_genes),
#         "--completeness {}".format(opts.checkv_completeness),
#         "--checkv_quality {}".format(opts.checkv_quality),
#         "--miuvig_quality {}".format(opts.miuvig_quality),
#         ")",
#         "&&",
#         "rm {}".format(os.path.join(directories["tmp"], "unbinned_contig_for_checkv.fasta")),
#     ]
#     return cmd



# # Symlink
# def get_symlink_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
#     # Command
#     cmd = ["("]
#     for filepath in input_filepaths:
#         cmd.append("ln -f -s {} {}".format(os.path.realpath(filepath), os.path.realpath(output_directory)))
#         cmd.append("&&")
#     cmd[-1] = ")"
#     return cmd

# ============
# Run Pipeline
# ============



def create_pipeline(opts, directories, f_cmds):

    # .................................................................
    # Primordial
    # .................................................................
    # Commands file
    pipeline = ExecutablePipeline(name=__program__, description=opts.name, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])

    # ==========
    # Coverage
    # ==========
    step = 1

    program = "coverage"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Calculating coverage for assembly"

    # i/o
    input_filepaths = [
        opts.bam,
    ]
    output_filenames = ["coverage_metabat.tsv", "coverage_noheader.tsv", "coverage_vamb.tsv"]
    # output_filenames = ["coverage_metabat.tsv", "coverage_noheader.tsv"]

    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_coverage_cmd(**params)
    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
                errors_ok=False,
    )

    # ==========
    # MetaBat2
    # ==========
    step  = 2

    program = "binning_metabat2"
    program_label = "{}__{}".format(step, program)
    
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Binning via MetaBat2 [All contigs]"

    # i/o
    input_filepaths = [
        opts.fasta,
        os.path.join(directories[("intermediate",  "1__coverage")], "coverage_metabat.tsv"),
    ]

    output_filenames = [
        # "bin*.fa", 
        "scaffolds_to_bins.tsv",
        ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        "prefix":"P",
    }

    cmd = get_metabat2_cmd(**params)
    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=False,
                errors_ok=False,
    )
    

    # =============
    # MaxBin2 (107)
    # =============
    step = 3

    program = "binning_maxbin2-107"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Binning via MaxBin2 [Marker Set=107]"

    # i/o
    input_filepaths = [
        opts.fasta,
        os.path.join(directories[("intermediate",  "1__coverage")], "coverage_noheader.tsv"),
    ]

    output_filenames = [
        # "bin*.fasta", 
        "scaffolds_to_bins.tsv",
        ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_maxbin2_107_cmd(**params)
    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=False,
                errors_ok=False,
    )

    # ============
    # MaxBin2 (40)
    # ============
    step = 4

    program = "binning_maxbin2-40"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Binning via MaxBin2 [Marker Set=40]"

    # i/o
    input_filepaths = [
        opts.fasta,
        os.path.join(directories[("intermediate",  "1__coverage")], "coverage_noheader.tsv"),
    ]

    output_filenames = [
        # "bin*.fasta", 
        "scaffolds_to_bins.tsv",
        ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_maxbin2_40_cmd(**params)
    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=False,
                errors_ok=False,
    )

    # ============
    # CONCOCT
    # ============
    step = 5

    program = "binning_concoct"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Binning via CONCOCT"

    # i/o
    input_filepaths = [
        opts.fasta,
        opts.bam,
    ]

    output_filepaths = [
        # os.path.join(output_directory, "bins/*"),
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),

    ]
    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_concoct_cmd(**params)
    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=False,
                errors_ok=False,
    )


    # PyTorch gets hung up when using subprocess
    # # ==========
    # # VAMB
    # # ==========
    # program = "binning_vamb"
    # # Add to directories
    # output_directory = directories[("intermediate",  program)] = create_directory(os.path.join(directories["intermediate"], "{}_output".format(program)))

    # # Info
    # step = 6
    # description = "Binnig via VAMB"

    # # i/o
    # input_filepaths = [
    #     opts.fasta, 
    #     os.path.join(directories[("intermediate",  "coverage")], "coverage_vamb.tsv"),
    # ]

    # output_filenames = ["bins/*.fna", "scaffolds_to_bins.tsv"]
    # output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    # params = {
    #     "input_filepaths":input_filepaths,
    #     "output_filepaths":output_filepaths,
    #     "output_directory":output_directory,
    #     "opts":opts,
    #     "directories":directories,
    # }

    # cmd = get_vamb_cmd(**params)
    # pipeline.add_step(
    #             id=program,
    #             description = description,
    #             step=step,
    #             cmd=cmd,
    #             input_filepaths = input_filepaths,
    #             output_filepaths = output_filepaths,
    #             validate_inputs=True,
    #             validate_outputs=False,
    #             errors_ok=False,
    # )

    # ==========
    # DAS_Tool
    # ==========
    step = 6

    program = "dastool"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Evaluation via DAS_Tool"
    # i/o
    input_filepaths = [
        os.path.join(directories[("intermediate",  "2__binning_metabat2")], "scaffolds_to_bins.tsv"),
        os.path.join(directories[("intermediate",  "3__binning_maxbin2-107")], "scaffolds_to_bins.tsv"),
        os.path.join(directories[("intermediate",  "4__binning_maxbin2-40")], "scaffolds_to_bins.tsv"),
        os.path.join(directories[("intermediate",  "5__binning_concoct")], "scaffolds_to_bins.tsv"),
        opts.fasta,
    ]

    output_filenames = [
        "__DASTool_bins",
        "__DASTool_summary.txt",
        "__DASTool_scaffolds2bin.txt",
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    } 

    cmd = get_dastool_cmd(**params)
    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=False,
                validate_outputs=True,
                errors_ok=False,
    )

    # ==========
    # Prodigal
    # ==========
    step = 7

    program = "gene-models_prodigal-prokaryote"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Gene calls via Prodigal [Prokaryote]"
    # i/o
    input_filepaths = [
        os.path.join(directories[("intermediate",  "6__dastool")], "__DASTool_bins"),
        os.path.join(directories[("intermediate",  "6__dastool")], "__DASTool_scaffolds2bin.txt"),

    ]

    output_filenames = [
        os.path.join(output_directory, "*.gff"),
        os.path.join(output_directory, "*.faa"),
        os.path.join(output_directory, "*.ffn"),
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_prodigal_cmd(**params)
    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
                errors_ok=False,
    )

   # ==========
    # CheckM
    # ==========
    step = 8

    program = "checkm"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info 
    description = "Evaluation via CheckM [Prokaryote]"
    # i/o
    input_filepaths = [
        directories[("intermediate",  "7__gene-models_prodigal-prokaryote")],
        os.path.join(directories[("intermediate",  "6__dastool")], "__DASTool_bins"),
    ]

    output_filenames = [
         "output.tsv",
         "filtered/*.list",
        #  "filtered/genomes/*.fa",
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_checkm_cmd(**params)
    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
                errors_ok=False,
    )

    # if not opts.no_virus:
    # ==========
        # # VirSorter2
        # # ==========
        # program = "viral_virsorter2"

        # # Add to directories
        # output_directory = directories[("intermediate",  program)] = create_directory(os.path.join(directories["intermediate"], "{}_output".format(program)))

        # # Info
        # step += 1
        # description = "Viral detection via VirSorter2 [VirSorter2]"
        # # i/o
        # input_filepaths = [
        #     os.path.join(directories[("intermediate",  "evaluation_checkm")],"filtered", "unbinned.list")
        # ]

        # output_filenames = [
        # ]
        # output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        # params = {
        #     "input_filepaths":input_filepaths,
        #     "output_filepaths":output_filepaths,
        #     "output_directory":output_directory,
        #     "opts":opts,
        #     "directories":directories,
        # }

        # cmd = get_virsorter2_cmd(**params)
        # pipeline.add_step(
        #             id=program,
        #             description = description,
        #             step=step,
        #             cmd=cmd,
        #             input_filepaths = input_filepaths,
        #             output_filepaths = output_filepaths,
        #             validate_inputs=True,
        #             validate_outputs=True,
        #             errors_ok=False,
        # )

        # # ViralVerify
        # # ==========
        # program = "viral_viralverify"

        # # Add to directories
        # output_directory = directories[("intermediate",  program)] = create_directory(os.path.join(directories["intermediate"], "{}_output".format(program)))

        # # Info
        # step += 1
        # description = "Viral detection via ViralVerify [Virus]"
        # # i/o
        # input_filepaths = [
        #     os.path.join(directories[("intermediate",  "evaluation_checkm")],"filtered", "unbinned.list")
        # ]

        # output_filenames = [
        # ]
        # output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        # params = {
        #     "input_filepaths":input_filepaths,
        #     "output_filepaths":output_filepaths,
        #     "output_directory":output_directory,
        #     "opts":opts,
        #     "directories":directories,
        # }

        # cmd = get_viralverify_cmd(**params)
        # pipeline.add_step(
        #             id=program,
        #             description = description,
        #             step=step,
        #             cmd=cmd,
        #             input_filepaths = input_filepaths,
        #             output_filepaths = output_filepaths,
        #             validate_inputs=True,
        #             validate_outputs=True,
        #             errors_ok=False,
        # )

    # CheckV
    # ==========
    step = 9 

    program = "checkv"

    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Viral verification with CheckV [Virus]"
    # i/o
    input_filepaths = [
        os.path.join(directories[("intermediate",  "8__checkm")],"filtered", "unbinned.list")
    ]

    output_filenames = [
        "quality_summary.tsv",
        "filtered/unbinned.list",
        "filtered/unbinned.fasta",
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_checkv_cmd(**params)
    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
                errors_ok=False,
    )

    # ==========
    # Prodigal
    # ==========
    step = 10

    program = "gene-models_prodigal-viral"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Gene calls via Prodigal [Viral]"
    # i/o
    input_filepaths = [
        os.path.join(directories[("intermediate",  "9__checkv")],"filtered", "genomes"),
        os.path.join(directories[("intermediate",  "9__checkv")],"filtered", "scaffolds_to_bins.tsv"),

    ]

    output_filenames = [
        os.path.join(output_directory, "*.gff"),
        os.path.join(output_directory, "*.faa"),
        os.path.join(output_directory, "*.ffn"),
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_prodigal_cmd(**params)
    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=False,
                validate_outputs=False,
                errors_ok=True,
    )
    
    # # ==========
    # # MetaBat2
    # # ==========
    # step  = 11

    # program = "binning_metabat2-eukaryotic"
    # program_label = "{}__{}".format(step, program)
    
    # # Add to directories
    # output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # # Info
    # description = "Binning via MetaBat2 [Eukaryotic]"

    # # i/o
    # input_filepaths = [
    #     os.path.join(directories[("intermediate",  "9__checkv")],"filtered","unbinned.fasta"),
    #     os.path.join(directories[("intermediate",  "1__coverage")], "coverage_metabat.tsv"),
    # ]

    # output_filenames = [
    #     # "bin*.fa", 
    #     "scaffolds_to_bins.tsv",
    #     ]
    # output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    # params = {
    #     "input_filepaths":input_filepaths,
    #     "output_filepaths":output_filepaths,
    #     "output_directory":output_directory,
    #     "opts":opts,
    #     "directories":directories,
    #     "prefix":"E",

    # }

    # cmd = get_metabat2_cmd(**params)
    # pipeline.add_step(
    #             id=program,
    #             description = description,
    #             step=step,
    #             cmd=cmd,
    #             input_filepaths = input_filepaths,
    #             output_filepaths = output_filepaths,
    #             validate_inputs=True,
    #             validate_outputs=False,
    #             errors_ok=False,
    # )

    
    

    # # =============
    # # Symlink
    # # =============
    # program = "symlink"
    # # Add to directories
    # output_directory = directories["output"]

    # # Info
    # step = 3
    # description = "Symlinking relevant output files"

    # # i/o
    # input_filepaths = [
    #     os.path.join(directories[("intermediate", "bowtie2")], "mapped.sorted.bam"),
    #     os.path.join(directories[("intermediate", "coverage")], "coverage.tsv"),

    # ]

    # output_filenames =  map(lambda fp: fp.split("/")[-1], input_filepaths)
    # output_filepaths = list(map(lambda fn:os.path.join(directories["output"], fn), output_filenames))
    #     # Prodigal
    #     # os.path.join(directories["output"], "*"),
    
    # params = {
    # "input_filepaths":input_filepaths,
    # "output_filepaths":output_filepaths,
    # "output_directory":output_directory,
    # "opts":opts,
    # "directories":directories,
    # }

    # cmd = get_symlink_cmd(**params)
    # pipeline.add_step(
    #         id=program,
    #         description = description,
    #         step=step,
    #         cmd=cmd,
    #         input_filepaths = input_filepaths,
    #         output_filepaths = output_filepaths,
    #         validate_inputs=True,
    #         validate_outputs=False,
    # )

    return pipeline


def create_viral_pipeline(opts, directories, f_cmds):

    # .................................................................
    # Primordial
    # .................................................................
    # Commands file
    pipeline = ExecutablePipeline(name="Viral", description=opts.name, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])
    
    # ==========
    # VirFinder
    # ==========
    step = 1

    program = "virfinder"

    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Viral identification with VirFinder [Viral]"
    # i/o
    input_filepaths = [
        os.path.join(directories[("intermediate",  "8__checkm")],"filtered", "unbinned.list")
    ]

    output_filenames = [
        "quality_summary.tsv",
        "filtered/unbinned.list",
        "filtered/unbinned.fasta",
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_checkv_cmd(**params)
    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
                errors_ok=False,
    )

    # CheckV
    # ==========
    step = 9 

    program = "checkv"

    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Viral verification with CheckV [Virus]"
    # i/o
    input_filepaths = [
        os.path.join(directories[("intermediate",  "8__checkm")],"filtered", "unbinned.list")
    ]

    output_filenames = [
        "quality_summary.tsv",
        "filtered/unbinned.list",
        "filtered/unbinned.fasta",
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_checkv_cmd(**params)
    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
                errors_ok=False,
    )

    # ==========
    # Prodigal
    # ==========
    step = 10

    program = "gene-models_prodigal-viral"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Gene calls via Prodigal [Viral]"
    # i/o
    input_filepaths = [
        os.path.join(directories[("intermediate",  "9__checkv")],"filtered", "genomes"),
        os.path.join(directories[("intermediate",  "9__checkv")],"filtered", "scaffolds_to_bins.tsv"),

    ]

    output_filenames = [
        os.path.join(output_directory, "*.gff"),
        os.path.join(output_directory, "*.faa"),
        os.path.join(output_directory, "*.ffn"),
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_prodigal_cmd(**params)
    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=False,
                validate_outputs=False,
                errors_ok=True,
    )
    

    return pipeline
# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """

    required_executables={
                # 1
                "jgi_summarize_bam_contig_depths", # coverm instead? 
                # 2
                "metabat2",
                # 3-4
                "run_MaxBin.pl",
                # 5
                "cut_up_fasta.py",
                "concoct_coverage_table.py",
                "concoct",
                "merge_cutup_clustering.py",
                "extract_fasta_bins.py",
                # 6
                # "vamb",
                # 6
                "DAS_Tool",
                # 7
                "prodigal",
                # 8
                "checkm",
                # 9 
                # "virsorter",
                # "viral_verify",
                "checkv",
                "seqkit",

                # "gtdbtk",
                # -1
                "scaffolds_to_bins.py",
                "check_scaffolds_to_bins.py",
                "partition_gene_models.py",
                "append_geneid_to_prodigal_gff.py",
                "filter_checkm_results.py",
                "filter_checkv_results.py",
                "VirFinder_wrapper.R",

     }

    if opts.path_config == "CONDA_PREFIX":
        executables = dict()
        for name in required_executables:
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
    accessory_scripts = [
        "scaffolds_to_bins.py",
        "check_scaffolds_to_bins.py",
        "partition_gene_models.py",
        "append_geneid_to_prodigal_gff.py",
        "filter_checkm_results.py",
        "filter_checkv_results.py",
        "VirFinder_wrapper.R",

    ]
    for name in accessory_scripts:
        executables[name] = "python " + os.path.join(opts.script_directory, "scripts", name)

    print(format_header( "Adding executables to path from the following source: {}".format(opts.path_config), "-"), file=sys.stdout)
    for name, executable in executables.items():
        if name in required_executables:
            print(name, executable, sep = " --> ", file=sys.stdout)
            os.environ[name] = executable.strip()
    print("", file=sys.stdout)


# Configure parameters
def configure_parameters(opts, directories):
    # assert opts.reference_assembly is not None, "Must include --reference_assembly"
    assert opts.checkm_tree in {"full", "reduced"}, "--checkm_tree must be either 'reduced' or 'full'"

    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -f <scaffolds.fasta> -b <mapped.sorted.bam> -n <name> -o <output_directory>".format(__program__)

    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-f","--fasta", type=str, required=True, help = "path/to/scaffolds.fasta")
    parser_io.add_argument("-b","--bam", type=str, required=True, help = "path/to/mapped.sorted.bam")
    parser_io.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    parser_io.add_argument("-o","--project_directory", type=str, default="veba_output/binning", help = "path/to/project_directory [Default: veba_output/binning]")
    parser_io.add_argument("-M", "--mode", type=str, default="all", help = "{prokaryotic, viral, eukaryotic, all} [Default: all]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')

    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
    parser_utility.add_argument("-I", "--iterations", type=int, default=-1, help = "Binning iterations for non-viral MAGs. -1 to go until completion. [Default: -1]")

    # parser_utility.add_argument("-c", "--CONDA_PREFIX", type=str, default=None, help = "Set a conda environment")


    # Binning
    parser_binning = parser.add_argument_group('Binning arguments')
    parser_binning.add_argument("-m", "--minimum_contig_length", type=int, default=1500, help="Minimum contig length.  Anything under 2500 will default to 2500 for MetaBat2 [Default: 1500] ")
    parser_binning.add_argument("--concoct_fragment_length", type=int, default=10000, help="CONCOCT | Fragment length [Default: 10000] ")
    parser_binning.add_argument("--concoct_overlap_length", type=int, default=0, help="CONCOCT | Fragment overlap length [Default: 0] ")
    parser_binning.add_argument("--maxbin2_options", type=str, default="", help="MaxBin2 | More options (e.g. --arg 1 ) [Default: ''] | https://sourceforge.net/projects/maxbin/")
    parser_binning.add_argument("--metabat2_options", type=str, default="", help="MetaBat2 | More options (e.g. --arg 1 ) [Default: ''] | https://bitbucket.org/berkeleylab/metabat/src/master/")
    parser_binning.add_argument("--concoct_options", type=str, default="", help="CONCOCT | More options (e.g. --arg 1 ) [Default: '']")
    # parser_binning.add_argument("--vamb_options", type=str, default="--minfasta 200000", help="vamb | More options (e.g. --arg 1 ) [Default: '--minfasta 200000']")

    parser_genemodels = parser.add_argument_group('Gene model arguments')
    parser_genemodels.add_argument("--prodigal_genetic_code", type=str, default=11, help="Prodigal -g translation table [Default: 11]")

    

    parser_evaluation = parser.add_argument_group('Evaluation arguments')
    parser_evaluation.add_argument("--dastool_searchengine", type=str, default="blast", help="DAS_Tool searchengine. [Default: blast] | https://github.com/cmks/DAS_Tool")
    parser_evaluation.add_argument("--dastool_options", type=str, default="", help="DAS_Tool | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/cmks/DAS_Tool")
    parser_evaluation.add_argument("--checkm_tree", type=str, default="reduced", help="CheckM tree type either 'reduced' or 'full' [Default: reduced]")
    parser_evaluation.add_argument("--checkm_completeness", type=float, default=50.0, help="CheckM completeness threshold [Default: 50]")
    parser_evaluation.add_argument("--checkm_contamination", type=float, default=10.0, help="CheckM contamination threshold [Default: 10]")
    parser_evaluation.add_argument("--checkm_strain_heterogeneity", type=float,  help="CheckM strain hetereogeneity threshold")
    parser_evaluation.add_argument("--checkm_options", type=str, default="", help="CheckM lineage_wf | More options (e.g. --arg 1 ) [Default: '']")

    parser_virus = parser.add_argument_group('Virus arguments')
    parser_virus.add_argument("--no_virus", action="store_true", help="Don't run viral detection")
    # parser_virus.add_argument("--include_provirus", action="store_true", help="Include provirus viral detection")
    # parser_virus.add_argument("--virsorter2_database", type=str, default="/usr/local/scratch/CORE/jespinoz/db/virsorter/v3/", help="VirSorter2 | More options (e.g. --arg 1 ) [Default: '/usr/local/scratch/CORE/jespinoz/db/virsorter/v3/']")
    # parser_virus.add_argument("--virsorter2_groups", type=str, default="dsDNAphage,NCLDV,ssDNA,lavidaviridae", help="VirSorter2 | More options (e.g. --arg 1 ) [Default: 'dsDNAphage,NCLDV,ssDNA,lavidaviridae']")
    # parser_virus.add_argument("--virsorter2_options", type=str, default="", help="VirSorter2 | More options (e.g. --arg 1 ) [Default: '']")
    # parser_virus.add_argument("--pfam_database", type=str, default="/usr/local/scratch/CORE/jespinoz/db/pfam/v33.1/Pfam-A.hmm", help="PFAM | More options (e.g. --arg 1 ) [Default: '/usr/local/scratch/CORE/jespinoz/db/pfam/v33.1/Pfam-A.hmm']")
    # parser_virus.add_argument("--viralverify_options", type=str, default="", help="ViralVerify | More options (e.g. --arg 1 ) [Default: '']")
    parser_virus.add_argument("--virfinder_pvalue", type=float, default=0.05 help="VirFinder p-value threshold [Default: 0.05]")
    parser_virus.add_argument("--checkv_database", type=str, default="/usr/local/scratch/CORE/jespinoz/db/checkv/checkv-db-v1.0", help="CheckV | More options (e.g. --arg 1 ) [Default: '/usr/local/scratch/CORE/jespinoz/db/checkv/checkv-db-v1.0']")
    parser_virus.add_argument("--checkv_options", type=str, default="", help="CheckV | More options (e.g. --arg 1 ) [Default: '']")
    parser_virus.add_argument("--multiplier_viral_to_host_genes", type=int, default=5, help = "Minimum number of viral genes [Default: 5]")
    parser_virus.add_argument("--checkv_completeness", type=float, default=50.0, help = "Minimum completeness [Default: 50.0]")
    parser_virus.add_argument("--checkv_quality", type=str, default="High-quality,Medium-quality,Complete", help = "Comma-separated string of acceptable arguments between {High-quality,Medium-quality,Complete} [Default: High-quality,Medium-quality,Complete]")
    parser_virus.add_argument("--miuvig_quality", type=str, default="High-quality,Medium-quality,Complete", help = "Comma-separated string of acceptable arguments between {High-quality,Medium-quality,Complete} [Default: High-quality,Medium-quality,Complete]")

    # Options
    opts = parser.parse_args()

    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    

    # Directories
    directories = dict()
    directories["project"] = create_directory(opts.project_directory)
    directories["sample"] = create_directory(os.path.join(directories["project"], opts.name))
    # directories["preprocessing"] = create_directory(os.path.join(directories["sample"], "preprocessing"))
    directories["output"] = create_directory(os.path.join(directories["sample"], "output"))
    directories["log"] = create_directory(os.path.join(directories["sample"], "log"))
    directories["tmp"] = create_directory(os.path.join(directories["sample"], "tmp"))
    directories["checkpoints"] = create_directory(os.path.join(directories["sample"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["sample"], "intermediate"))


    # Info
    print(format_header(__program__, "="), file=sys.stdout)
    print(format_header("Configuration:", "-"), file=sys.stdout)
    print(format_header("Name: {}".format(opts.name), "."), file=sys.stdout)
    print("Python version:", sys.version.replace("\n"," "), file=sys.stdout)
    print("Python path:", sys.executable, file=sys.stdout) #sys.path[2]
    print("Script version:", __version__, file=sys.stdout)
    print("Moment:", get_timestamp(), file=sys.stdout)
    print("Directory:", os.getcwd(), file=sys.stdout)
    print("Commands:", list(filter(bool,sys.argv)),  sep="\n", file=sys.stdout)
    configure_parameters(opts, directories)
    sys.stdout.flush()

    # if opts.CONDA_PREFIX:


    # Run pipeline
    with open(os.path.join(directories["sample"], "commands.sh"), "w") as f_cmds:
        pipeline = create_pipeline(
                    opts=opts,
                    directories=directories,
                    f_cmds=f_cmds,
        )
        pipeline.compile()
        pipeline.execute(restart_from_checkpoint=opts.restart_from_checkpoint)
   

if __name__ == "__main__":
    main()
