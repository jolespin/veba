#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, random
from collections import OrderedDict, defaultdict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.3.15"

# Assembly
def get_coverage_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [


    # For pipeline purposes
    "cat",
    opts.fasta,
    "|",
    os.environ["seqkit"],
    "seq",
    "-m {}".format(opts.minimum_contig_length),
    ">" ,
    os.path.join(directories["output"], "unbinned.fasta"),

    "&&",

    # Coverage for MetaBat2 
    os.environ["coverm"],
    "contig",
    "--threads {}".format(opts.n_jobs),
    "--methods metabat",
    "--bam-files",
    " ".join(opts.bam),
    ">",
    output_filepaths[0],

    "&&",

    # # Coverage for MaxBin2
    # "cut -f1,4",
    # output_filepaths[0],
    # "|",
    # "tail -n +2",
    # ">",
    # output_filepaths[1],
    """
    python -c "import pandas as pd; df = pd.read_csv('{}', sep='\t', index_col=0); df.loc[:,df.columns.map(lambda x: x.startswith('mapped.sorted.bam') and (not '-var' in x))].to_csv('{}', sep='\t', header=None)"
    """.format( 
            output_filepaths[0],
            output_filepaths[1],
    )
    ]

    return cmd


# Prodigal
def get_prodigal_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [
        "cat",
        input_filepaths[0],
        "|",
        os.environ["seqkit"],
        "seq",
        "-m {}".format(opts.minimum_contig_length),
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

    ]
    return cmd

# Metabat2
def get_metabat2_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, prefix, seed):

    cmd = [
        os.environ["binning_wrapper.py"],
        "-a metabat2",
        "-f {}".format(input_filepaths[0]), # scaffolds.fasta
        "-c {}".format(input_filepaths[1]),
        "-o {}".format(output_directory),
        "-m {}".format(max(opts.minimum_contig_length, 1500)), # mininimum contig length must be >= 1500 nts for MetaBat2
        "-s {}".format(opts.minimum_genome_length), 
        "--n_jobs {}".format(opts.n_jobs),
        "--random_state {}".format(seed),
        "--bin_prefix {}".format(prefix),
        "--remove_bins",
        "--remove_intermediate_files",
        "--metabat2_options {}".format(opts.metabat2_options) if bool(opts.metabat2_options) else "",
    ]
    return cmd


# MaxBin2
def get_maxbin2_107_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, prefix, seed):

    cmd = [
        os.environ["binning_wrapper.py"],
        "-a maxbin2",
        "--maxbin2_markerset 107",
        "-f {}".format(input_filepaths[0]), # scaffolds.fasta
        "-c {}".format(input_filepaths[1]),
        "-o {}".format(output_directory),
        "-m {}".format(opts.minimum_contig_length),
        "-s {}".format(opts.minimum_genome_length), 
        "--n_jobs {}".format(opts.n_jobs),
        "--random_state {}".format(seed),
        "--bin_prefix {}".format(prefix),
        "--remove_bins",
        "--remove_intermediate_files",
        "--maxbin2_options {}".format(opts.maxbin2_options) if bool(opts.maxbin2_options) else "",

    ]
    return cmd


def get_maxbin2_40_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, prefix, seed):

    cmd = [
        os.environ["binning_wrapper.py"],
        "-a maxbin2",
        "--maxbin2_markerset 40",
        "-f {}".format(input_filepaths[0]), # scaffolds.fasta
        "-c {}".format(input_filepaths[1]),
        "-o {}".format(output_directory),
        "-m {}".format(opts.minimum_contig_length),
        "-s {}".format(opts.minimum_genome_length), 
        "--n_jobs {}".format(opts.n_jobs),
        "--random_state {}".format(seed),
        "--bin_prefix {}".format(prefix),
        "--remove_bins",
        "--remove_intermediate_files",
        "--maxbin2_options {}".format(opts.maxbin2_options) if bool(opts.maxbin2_options) else "",
    ]
    return cmd

def get_maxbin2_null_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, prefix, seed):
    cmd = [
        # Create dummy scaffolds_to_bins.tsv to overwrite later. This makes DAS_Tool easier to run
        "echo '' > {}".format(os.path.join(output_directory, "scaffolds_to_bins.tsv")),
        "&&",
        "echo 'Skipping MaxBin2'",
    ]
    return cmd

# CONCOCT
def get_concoct_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, prefix, seed):

    cmd = [
        os.environ["binning_wrapper.py"],
        "--concoct_fragment_length {}".format(opts.concoct_fragment_length),
        "--concoct_overlap_length {}".format(opts.concoct_overlap_length),
        "-a concoct",
        "-f {}".format(input_filepaths[0]), # scaffolds.fasta
        "-b {}".format(" ".join(opts.bam)),
        "-o {}".format(output_directory),
        "-m {}".format(opts.minimum_contig_length),
        "-s {}".format(opts.minimum_genome_length), 
        "--n_jobs {}".format(opts.n_jobs),
        "--random_state {}".format(seed),
        "--bin_prefix {}".format(prefix),
        "--remove_bins",
        "--remove_intermediate_files",
        "--concoct_options {}".format(opts.concoct_options) if bool(opts.concoct_options) else "",

    ]
    return cmd

def get_concoct_null_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, prefix, seed):
    cmd = [
        # Create dummy scaffolds_to_bins.tsv to overwrite later. This makes DAS_Tool easier to run
        "echo '' > {}".format(os.path.join(output_directory, "scaffolds_to_bins.tsv")),
        "&&",
        "echo 'Skipping CONCOCT'",
    ]
    return cmd


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
        os.environ["DAS_Tool"],
        "--bins ${S2B_ARRAY[0]}", # scaffolds.fasta
        "--contigs {}".format(input_filepaths[4]),
        "--outputbasename {}".format(os.path.join(output_directory, "_")),
        "--labels ${S2B_ARRAY[1]}",
        "--search_engine {}".format(opts.dastool_searchengine),
        "--score_threshold {}".format(opts.dastool_minimum_score),
        "--write_bins 1",
        "--create_plots 0",
        "--threads {}".format(opts.n_jobs),
        "--proteins {}".format(os.path.join(directories[("intermediate",  "2__prodigal")], "gene_models.faa")),
        "--debug",
        opts.dastool_options,
        # Eukaryotic
        "&&",
        "cat",
        os.path.join(output_directory, "__DASTool_bins", "*.fa"),
        "|",
        os.environ["seqkit"],
        "seq",
        "-m {}".format(opts.tiara_minimum_length),
        ">" ,
        os.path.join(directories["tmp"], "scaffolds.binned.gte{}.fasta".format(opts.tiara_minimum_length)),


        # Tiara (Now that CheckM2 is being used, is this necessary?)
        "&&",
        "mkdir -p {}".format(os.path.join(output_directory, "consensus_domain_classification")),
        "&&",
        "mkdir -p {}".format(os.path.join(output_directory, "__DASTool_bins", "eukaryota")),
        "&&",
        os.environ["tiara"],
        "-i {}".format(os.path.join(directories["tmp"], "scaffolds.binned.gte{}.fasta".format(opts.tiara_minimum_length))),
        "-o {}".format(os.path.join(output_directory, "consensus_domain_classification", "tiara_output.tsv")),
        "--probabilities",
        "-m {}".format(opts.tiara_minimum_length),
        "-t {}".format(opts.n_jobs),
        opts.tiara_options,

        # Predict domain
        "&&",
        os.environ["consensus_domain_classification.py"],
        "-i {}".format(os.path.join(output_directory, "__DASTool_scaffolds2bin.txt")),
        "-t {}".format(os.path.join(output_directory, "consensus_domain_classification", "tiara_output.tsv")),
        "-o {}".format(os.path.join(output_directory, "consensus_domain_classification")),
        "--logit_transform {}".format(opts.logit_transform),

        # Move eukaryota
        "&&",
        "for ID_GENOME in $(cat %s); do mv %s/${ID_GENOME}.* %s; done"%(
            os.path.join(output_directory, "consensus_domain_classification", "eukaryota.list"),
            os.path.join(output_directory, "__DASTool_bins"),
            os.path.join(output_directory, "__DASTool_bins", "eukaryota"),
            ),
        # Eukaryotic scaffolds
        "&&",
        "cat",
        os.path.join(output_directory, "__DASTool_bins", "eukaryota", "*.fa"),
        "|",
        "grep",
        '"^>"',
        "|",
        "cut -c2-",
        ">",
        os.path.join(output_directory, "__DASTool_bins", "eukaryota", "eukaryota.scaffolds.list"),

        # Remove eukaryotic scaffolds
        "&&",
        "cat",
        os.path.join(output_directory, "__DASTool_scaffolds2bin.txt"),
        "|",
        "grep",
        "-v",
        "-f {}".format(os.path.join(output_directory, "__DASTool_bins", "eukaryota", "eukaryota.scaffolds.list")),
        ">",
        os.path.join(output_directory, "__DASTool_scaffolds2bin.no_eukaryota.txt"),
    
        # Partition
        "&&",
        "cut",
        "-f1",
        os.path.join(output_directory, "__DASTool_scaffolds2bin.no_eukaryota.txt"),
        ">",
        os.path.join(output_directory, "binned.list"),
        "&&",
        os.environ["partition_gene_models.py"],
        "-i {}".format(os.path.join(output_directory, "__DASTool_scaffolds2bin.no_eukaryota.txt")),
        "-g {}".format(os.path.join(directories[("intermediate",  "2__prodigal")], "gene_models.gff")),
        "-d {}".format(os.path.join(directories[("intermediate",  "2__prodigal")], "gene_models.ffn")),
        "-a {}".format(os.path.join(directories[("intermediate",  "2__prodigal")], "gene_models.faa")),
        "-o {}".format(os.path.join(output_directory, "__DASTool_bins")),
        "--use_mag_as_description",

        "&&",
        "rm -rf {} {}".format(
            os.path.join(output_directory, "_.seqlength"), 
            os.path.join(directories["tmp"], "scaffolds.binned.gte{}.fasta".format(opts.tiara_minimum_length)),
        ),

    ]
    return cmd



# CheckM2
def get_checkm2_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    # checkm2 predict -i *.faa -o checkm2_output --genes -x faa --threads 16 --force
    cmd = [
        "mkdir -p {}".format(os.path.join(directories["tmp"], "checkm2")),
          "&&",

        os.environ["checkm2"],
        "predict",
        "-i {}".format(os.path.join(input_filepaths[0], "*.faa")),
        "--genes", 
        "-o {}".format(output_directory),
        "-t {}".format(opts.n_jobs),
        "--force",
        "-x faa",
        # "--tmpdir {}".format(os.environ["TMPDIR"] if "TMPDIR" in os.environ else directories["tmp"]), # Hack around: OSError: AF_UNIX path too long
        "--tmpdir {}".format(os.path.join(directories["tmp"], "checkm2")),
        "--database_path {}".format(os.path.join(os.environ["VEBA_DATABASE"], "Classify", "CheckM2", "uniref100.KO.1.dmnd")),
        opts.checkm2_options,
            "&&",

        "gzip {}".format(os.path.join(output_directory, "diamond_output", "*.tsv")),

            "&&",

        os.environ["filter_checkm2_results.py"],
        "-i {}".format(os.path.join(output_directory, "quality_report.tsv")),
        "-b {}".format(input_filepaths[0]),
        "-o {}".format(os.path.join(output_directory, "filtered")),
        "-f {}".format(opts.fasta),
        "-m {}".format(opts.minimum_contig_length),
        # "--unbinned",
        "--completeness {}".format(opts.checkm2_completeness),
        "--contamination {}".format(opts.checkm2_contamination),
        "-x fa",
            "&&",

        os.environ["scaffolds_to_bins.py"],
        "-x fa",
        "-i {}".format(os.path.join(output_directory, "filtered", "genomes")),
        ">",
        os.path.join(output_directory, "filtered", "scaffolds_to_bins.tsv"),
            "&&",

        "cat",
        input_filepaths[1],
        "|",
        os.environ["seqkit"],
        "seq",
        "-m {}".format(opts.minimum_contig_length),
        "|",
        os.environ["seqkit"],
        "grep",
        "--pattern-file {}".format(os.path.join(output_directory, "filtered", "unbinned.list")),
        ">",
        output_filepaths[-1],
            "&&",

        "rm -rf {}".format(os.path.join(output_directory, "protein_files")),

    ]
    return cmd


def get_featurecounts_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # ORF-Level Counts
    cmd = [
    "mkdir -p {}".format(os.path.join(directories["tmp"], "featurecounts")),
    "&&",
        os.environ["featureCounts"],
        "-G {}".format(input_filepaths[0]),
        "-a {}".format(input_filepaths[1]),
        "-o {}".format(os.path.join(output_directory, "featurecounts.orfs.tsv")),
        "-F GTF",
        "--tmpDir {}".format(os.path.join(directories["tmp"], "featurecounts")),
        "-T {}".format(opts.n_jobs),
        "-g gene_id",
        "-t CDS",
        opts.featurecounts_options,
        " ".join(opts.bam),
    "&&",
    "gzip -f {}".format(os.path.join(output_directory, "featurecounts.orfs.tsv")),
    ]
    return cmd


def get_consolidate_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
        # os.path.join(directories["intermediate"], "*__checkm2",  "filtered", "scaffolds_to_bins.tsv"),
        # os.path.join(directories["intermediate"], "*__checkm2",  "filtered", "bins.list"),
        # os.path.join(directories["intermediate"], "*__checkm2",  "filtered", "binned.list"),
        # os.path.join(directories["intermediate"], "*__checkm2",  "filtered", "checkm2_results.filtered.tsv"),
        # os.path.join(directories["intermediate"], "*__checkm2", "filtered", "genomes"),
        # os.path.join(directories[("intermediate", "{}__featurecounts".format(step-1))], "featurecounts.orfs.tsv.gz"),
    cmd = [
        "rm -rf {}".format(os.path.join(output_directory, "*")),
    ]

    cmd = [ 
        "mkdir -p {}".format(os.path.join(output_directory, "genomes")),

            "&&",   

        # scaffolds_to_bins.tsv
        "cat",
        input_filepaths[0],
        ">",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),

            "&&",
        # bins.list
        "cat",
        input_filepaths[1],
        ">",
        os.path.join(output_directory, "bins.list"),

            "&&",

        # binned.list
        "cat",
        input_filepaths[2],
        ">",
        os.path.join(output_directory, "binned.list"),

            "&&",

        # checkm2_results.filtered.tsv
        os.environ["concatenate_dataframes.py"],
        "-a 0",
        # "-e",
        input_filepaths[3],
        ">",
        os.path.join(output_directory, "checkm2_results.filtered.tsv"),
    ]

    # # Genomes
    # for fp in glob.glob(input_filepaths[4]):
    #     fn = os.path.split(fp)[1]
    #     src = os.path.realpath(os.path.join(fp))
    #     dst = os.path.join(output_directory,"genomes", fn)
    #     cmd += [ 
    #         "&&",
    #         "ln -sf",
    #         src,
    #         dst,
    #     ]
    
        # Genomes

    cmd += [ 
        "&&",
        # "cp -rf",
        "ln -sfr",
        input_filepaths[4],
        os.path.join(output_directory,"genomes"),
    ]
    
    # # Featurecounts
    # for fp in [input_filepaths[5]]:
    #     fn = os.path.split(fp)[1]
    #     src = os.path.realpath(os.path.join(fp))
    #     dst = os.path.join(output_directory, fn)
    #     cmd += [ 
    #         "&&",
    #         "ln -sf",
    #         src,
    #         dst,
    #     ]

    # Featurecounts
    cmd += [
            "&&",
            # "ln -sfr",
            "cp -rf",
            input_filepaths[5],
            output_directory,
        ]
        
    # SeqKit
    cmd += [ 
            "&&",
        # Statistics
        os.environ["seqkit"],
        "stats",
        "-a",
        "-b",
        "-T",
        "-j {}".format(opts.n_jobs),
        os.path.join(output_directory, "genomes", "*.fa"),

        "|",

        """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: x[:-3]); df.to_csv(sys.stdout, sep="\t")'"""
        ">",
        os.path.join(output_directory,"genome_statistics.tsv"),

            "&&",

        # Statistics
        "cat",
        opts.fasta,
        "|",
        os.environ["seqkit"],
        "grep",
        "-v",
        "-f {}".format(os.path.join(output_directory, "binned.list")),
        "-j {}".format(opts.n_jobs),
        "|",
        os.environ["seqkit"],
        "seq",
        "-j {}".format(opts.n_jobs),
        "-m {}".format(opts.minimum_contig_length),
        ">",
        os.path.join(output_directory,"unbinned.fasta"),
    

        # "&&",
        # "rm -rf {}".format(os.path.join(directories["tmp"],"*")),

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
    description = "Calculating coverage for assembly via CoverM"

    input_filepaths = [
            *opts.bam,
        ]

    output_filenames = ["coverage_metabat.tsv", "coverage_noheader.tsv"]

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
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
                errors_ok=False,
                log_prefix=program_label,

    )

    # ==========
    # Prodigal
    # ==========
    step = 2

    program = "prodigal"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Gene calls via Prodigal"
    # i/o
    input_filepaths = [
        opts.fasta,
    ]

    output_filenames = [
        "gene_models.gff",
        "gene_models.faa",
        "gene_models.ffn",
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
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
                errors_ok=False,
                log_prefix=program_label,

    )

    # Iterative mode
    input_fasta = opts.fasta
    for iteration in range(1, opts.n_iter+1):
        steps = dict()

        seed = opts.random_state
        if opts.random_state == -1:
            seed = random.choice(range(1,10000))
        if opts.random_state == 0:
            seed = iteration

        # ==========
        # MetaBat2
        # ==========
        step  += 1

        program = "binning_metabat2"
        program_label = "{}__{}".format(step, program)
        
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


        # Info
        description = "Binning via MetaBat2 [Iteration={}]".format(iteration)

        # i/o
        input_filepaths = [
            input_fasta,
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
            "prefix":"{}__METABAT2__{}.{}__".format(opts.name, "P", iteration),
            "seed":seed,
        }



        cmd = get_metabat2_cmd(**params)
        pipeline.add_step(
                    id=program_label,
                    description = description,
                    step=step,
                    cmd=cmd,
                    input_filepaths = input_filepaths,
                    output_filepaths = output_filepaths,
                    validate_inputs=False,
                    validate_outputs=False,
                    errors_ok=False,
                    acceptable_returncodes={0,1,2},
                    log_prefix=program_label,

        )

        steps[program] = step


        # =============
        # MaxBin2 (107)
        # =============
        step += 1

        program = "binning_maxbin2-107"

        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


        # Info
        if opts.skip_maxbin2:
            description = "[Skipping] Binning via MaxBin2 [Marker Set=107] [Iteration={}]".format(iteration)
        else:
            description = "Binning via MaxBin2 [Marker Set=107] [Iteration={}]".format(iteration)

        # i/o
        input_filepaths = [
            input_fasta,
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
            "prefix":"{}__MAXBIN2-107__{}.{}__".format(opts.name, "P", iteration),
            "seed":seed,

        }

        if opts.skip_maxbin2:
            cmd = get_maxbin2_null_cmd(**params)
        else:
            cmd = get_maxbin2_107_cmd(**params)
        pipeline.add_step(
                    id=program_label,
                    description = description,
                    step=step,
                    cmd=cmd,
                    input_filepaths = input_filepaths,
                    output_filepaths = output_filepaths,
                    validate_inputs=False,
                    validate_outputs=False,
                    errors_ok=False,
                    acceptable_returncodes={0,1,2},                    
                    log_prefix=program_label,

        )

        steps[program] = step

        # ============
        # MaxBin2 (40)
        # ============
        step  += 1

        program = "binning_maxbin2-40"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


        # Info
        if opts.skip_maxbin2:
            description = "[Skipping] Binning via MaxBin2 [Marker Set=40] [Iteration={}]".format(iteration)

        else:
            description = "Binning via MaxBin2 [Marker Set=40] [Iteration={}]".format(iteration)

        # i/o
        input_filepaths = [
            input_fasta,
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
            "prefix":"{}__MAXBIN2-40__{}.{}__".format(opts.name, "P", iteration),
            "seed":seed,

        }
        if opts.skip_maxbin2:
            cmd = get_maxbin2_null_cmd(**params)
        else:
            cmd = get_maxbin2_40_cmd(**params)
        pipeline.add_step(
                    id=program_label,
                    description = description,
                    step=step,
                    cmd=cmd,
                    input_filepaths = input_filepaths,
                    output_filepaths = output_filepaths,
                    validate_inputs=False,
                    validate_outputs=False,
                    errors_ok=False,
                    acceptable_returncodes={0,1,2},                    
                    log_prefix=program_label,

        )

        steps[program] = step

        # ============
        # CONCOCT
        # ============
        step  += 1

        program = "binning_concoct"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


        # Info
        if opts.skip_concoct:
            description = "[Skipping] Binning via CONCOCT [Iteration={}]".format(iteration)
        else:
            description = "Binning via CONCOCT [Iteration={}]".format(iteration)

        # i/o
        input_filepaths = [
            input_fasta,
            *opts.bam,
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
            "prefix":"{}__CONCOCT__{}.{}__".format(opts.name, "P", iteration),
            "seed":seed,
        }

        if opts.skip_concoct:
            cmd = get_concoct_null_cmd(**params)
        else:
            cmd = get_concoct_cmd(**params)
        pipeline.add_step(
                    id=program_label,
                    description = description,
                    step=step,
                    cmd=cmd,
                    input_filepaths = input_filepaths,
                    output_filepaths = output_filepaths,
                    validate_inputs=False,
                    validate_outputs=False,
                    errors_ok=False,
                    acceptable_returncodes={0,1,2},                    
                    log_prefix=program_label,

        )

        steps[program] = step

        # ==========
        # DAS_Tool
        # ==========
        step  += 1

        program = "dastool"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


        # Info
        description = "Evaluation via DAS_Tool [Iteration={}]".format(iteration)
        # i/o
        input_filepaths = [
            os.path.join(directories[("intermediate",  "{}__binning_metabat2".format(steps["binning_metabat2"]))], "scaffolds_to_bins.tsv"),
            os.path.join(directories[("intermediate",  "{}__binning_maxbin2-107".format(steps["binning_maxbin2-107"]))], "scaffolds_to_bins.tsv"),
            os.path.join(directories[("intermediate",  "{}__binning_maxbin2-40".format(steps["binning_maxbin2-40"]))], "scaffolds_to_bins.tsv"),
            os.path.join(directories[("intermediate",  "{}__binning_concoct".format(steps["binning_concoct"]))], "scaffolds_to_bins.tsv"),
            input_fasta,
            # os.path.join(directories[("intermediate",  "2__prodigal")], "gene_models.faa"),
        ]

        output_filenames = [
            "__DASTool_bins",
            "__DASTool_summary.txt",
            "__DASTool_scaffolds2bin.no_eukaryota.txt",
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
                    id=program_label,
                    description = description,
                    step=step,
                    cmd=cmd,
                    input_filepaths = input_filepaths,
                    output_filepaths = output_filepaths,
                    validate_inputs=False,
                    validate_outputs=False,
                    errors_ok=False,
                    acceptable_returncodes={0,1,2},                    
                    log_prefix=program_label,
                    # acceptable_returncodes= {0,1},

        )


        steps[program] = step

        # ==========
        # CheckM2
        # ==========
        step  += 1

        program = "checkm2"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


        # Info 
        description = "Evaluation via CheckM2 [Iteration={}]".format(iteration)
        # i/o
        input_filepaths = [
            os.path.join(directories[("intermediate",  "{}__dastool".format(steps["dastool"]))], "__DASTool_bins"),
            input_fasta,
        ]

        output_filenames = [
            "filtered/checkm2_results.filtered.tsv",
            "filtered/*.list",
        ]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))
        output_filepaths += [ 
            os.path.join(directories["tmp"], "unbinned_{}.fasta".format(iteration)),
        ]

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_checkm2_cmd(**params)
        pipeline.add_step(
                    id=program_label,
                    description = description,
                    step=step,
                    cmd=cmd,
                    input_filepaths = input_filepaths,
                    output_filepaths = output_filepaths,
                    validate_inputs=False,
                    validate_outputs=False,
                    errors_ok=False,
                    acceptable_returncodes={0,1,2},                    
                    log_prefix=program_label,
                    # acceptable_returncodes= {0,1},

        )

        steps[program] = step

        # Reset input_fasta
        input_fasta = os.path.join(directories["tmp"], "unbinned_{}.fasta".format(iteration))


    # ==========
    # featureCounts
    # ==========
    step += 1

    # Info
    program = "featurecounts"
    program_label = "{}__{}".format(step, program)
    description = "Counting reads"

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # i/o
    input_filepaths = [ 
        opts.fasta,
        os.path.join(directories[("intermediate",  "2__prodigal")], "gene_models.gff"),
        *opts.bam,
    ]

    output_filenames = ["featurecounts.orfs.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_featurecounts_cmd(**params)
    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
                errors_ok=False,
                acceptable_returncodes={0,1,2},                    
                log_prefix=program_label,
    )
    
    
    # =============
    # Consolidate
    # =============
    step += 1

    program = "consolidate"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories["output"]

    # Info

    description = "Consolidate output files"

    # i/o
    input_filepaths = [
        os.path.join(directories["intermediate"], "*__checkm2",  "filtered", "scaffolds_to_bins.tsv"),
        os.path.join(directories["intermediate"], "*__checkm2",  "filtered", "bins.list"),
        os.path.join(directories["intermediate"], "*__checkm2",  "filtered", "binned.list"),
        os.path.join(directories["intermediate"], "*__checkm2",  "filtered", "checkm2_results.filtered.tsv"),
        os.path.join(directories["intermediate"], "*__checkm2", "filtered", "genomes", "*"),
        os.path.join(directories[("intermediate", "{}__featurecounts".format(step-1))], "featurecounts.orfs.tsv.gz"),
    ]

    output_filenames =  [
        "scaffolds_to_bins.tsv", 
        "bins.list",
        "binned.list",
        "unbinned.fasta",
        "genomes",
        "checkm2_results.filtered.tsv",
        "featurecounts.orfs.tsv.gz",
        "genome_statistics.tsv",
    ]

    output_filepaths = list(map(lambda fn:os.path.join(directories["output"], fn), output_filenames))

    
    params = {
    "input_filepaths":input_filepaths,
    "output_filepaths":output_filepaths,
    "output_directory":output_directory,
    "opts":opts,
    "directories":directories,
    }

    cmd = get_consolidate_cmd(**params)
    pipeline.add_step(
            id=program_label,
            description = description,
            step=step,
            cmd=cmd,
            input_filepaths = input_filepaths,
            output_filepaths = output_filepaths,
            validate_inputs=False,
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
                "binning_wrapper.py",
                "scaffolds_to_bins.py",
                "check_scaffolds_to_bins.py",
                "partition_gene_models.py",
                "append_geneid_to_prodigal_gff.py",
                "filter_checkm2_results.py",
                "consensus_domain_classification.py",
                "concatenate_dataframes.py",
                "subset_table.py",
                }

    required_executables={
                # 1
                "coverm",
                # 2
                "metabat2",
                # 6
                # 6
                "DAS_Tool",
                # 7
                "prodigal",
                # 8
                "checkm2",
                # 9 
                "seqkit",
                # 10 
                "featureCounts",
                # 11
                "tiara",
 
     } | accessory_scripts

    if not opts.skip_maxbin2:
        # 3-4
        required_executables |= {
            "run_MaxBin.pl",
            }
    if not opts.skip_concoct:
        # 5
        required_executables |= { 
            "cut_up_fasta.py",
            "concoct_coverage_table.py",
            "concoct",
            "merge_cutup_clustering.py",
            "extract_fasta_bins.py",
            }

    if opts.path_config == "CONDA_PREFIX":
        executables = dict()
        for name in sorted(required_executables):
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
        if name.endswith(".py"):
            executables[name] = "python " + os.path.join(opts.script_directory, "scripts", name)
        else:
            executables[name] =  os.path.join(opts.script_directory, "scripts", name)


    print(format_header( "Adding executables to path from the following source: {}".format(opts.path_config), "-"), file=sys.stdout)
    for name, executable in executables.items():
        if name in required_executables:
            print(name, executable, sep = " --> ", file=sys.stdout)
            os.environ[name] = executable.strip()
    print("", file=sys.stdout)


# Configure parameters
def configure_parameters(opts, directories):
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
    parser_io.add_argument("-b","--bam", type=str, nargs="+", required=True, help = "path/to/mapped.sorted.bam files separated by spaces.")
    parser_io.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    parser_io.add_argument("-o","--project_directory", type=str, default="veba_output/binning/prokaryotic/", help = "path/to/project_directory [Default: veba_output/binning/prokaryotic/]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')

    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("-I", "--n_iter", type=int, default=3, help = "Number of iterations to run binning. Use -1 for brute force (-1 isn't ready yet) [Default: 3]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Use -1 for completely random. Use 0 for consecutive random states.  Use any other positive integer for the same random state for all iterations [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
    parser_utility.add_argument("--tmpdir", type=str, help="Set temporary directory") 

    # Databases
    parser_databases = parser.add_argument_group('Database arguments')
    parser_databases.add_argument("--veba_database", type=str,  help=f"VEBA database location.  [Default: $VEBA_DATABASE environment variable]")

    # Binning
    parser_binning = parser.add_argument_group('Binning arguments')
    parser_binning.add_argument("-m", "--minimum_contig_length", type=int, default=1500, help="Minimum contig length.  Anything under 2500 will default to 2500 for MetaBat2 [Default: 1500] ")
    parser_binning.add_argument("-s", "--minimum_genome_length", type=int, default=150000, help="Minimum genome length.  [Default: 150000]")
    parser_binning.add_argument("--concoct_fragment_length", type=int, default=10000, help="CONCOCT | Fragment length [Default: 10000] ")
    parser_binning.add_argument("--concoct_overlap_length", type=int, default=0, help="CONCOCT | Fragment overlap length [Default: 0] ")
    parser_binning.add_argument("--skip_maxbin2", action="store_true", help="MaxBin2 | Skip MaxBin2. Useful for large datasets")
    parser_binning.add_argument("--skip_concoct", action="store_true", help="CONCOCT | Skip CONCOCT. Useful when there's a lot of samples")
    parser_binning.add_argument("--maxbin2_options", type=str, default="", help="MaxBin2 | More options (e.g. --arg 1 ) [Default: ''] | https://sourceforge.net/projects/maxbin/")
    parser_binning.add_argument("--metabat2_options", type=str, default="", help="MetaBat2 | More options (e.g. --arg 1 ) [Default: ''] | https://bitbucket.org/berkeleylab/metabat/src/master/")
    parser_binning.add_argument("--concoct_options", type=str, default="", help="CONCOCT | More options (e.g. --arg 1 ) [Default: '']")
    # parser_binning.add_argument("--vamb_options", type=str, default="", help="vamb | More options (e.g. --arg 1 ) [Default: '']")
    # parser_binning.add_argument("--metacoag_options", type=str, default="", help="metacoag | More options (e.g. --arg 1 ) [Default: '']")

    parser_genemodels = parser.add_argument_group('Gene model arguments')
    parser_genemodels.add_argument("--prodigal_genetic_code", type=str, default=11, help="Prodigal -g translation table [Default: 11]")

    parser_evaluation = parser.add_argument_group('Evaluation arguments')
    parser_evaluation.add_argument("--dastool_searchengine", type=str, default="diamond", help="DAS_Tool searchengine. [Default: diamond] | https://github.com/cmks/DAS_Tool")
    parser_evaluation.add_argument("--dastool_minimum_score", type=float, default=0.1, help="DAS_Tool score_threshold. Score threshold until selection algorithm will keep selecting bins. This is set to a relaxed setting because CheckM2 is run post hoc. [Default: 0.1] | https://github.com/cmks/DAS_Tool")
    parser_evaluation.add_argument("--dastool_options", type=str, default="", help="DAS_Tool | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/cmks/DAS_Tool")
    parser_evaluation.add_argument("--checkm2_completeness", type=float, default=50.0, help="CheckM2 completeness threshold [Default: 50.0]")
    parser_evaluation.add_argument("--checkm2_contamination", type=float, default=10.0, help="CheckM2 contamination threshold [Default: 10.0]")
    parser_evaluation.add_argument("--checkm2_options", type=str, default="", help="CheckM lineage_wf | More options (e.g. --arg 1 ) [Default: '']")

    # featureCounts
    parser_featurecounts = parser.add_argument_group('featureCounts arguments')
    parser_featurecounts.add_argument("--featurecounts_options", type=str, default="", help="featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/")

    # Tiara
    parser_domain = parser.add_argument_group('Domain classification arguments')
    parser_domain.add_argument("--logit_transform", type=str, default="softmax", help = " Transformation for consensus_domain_classification: {softmax, tss} [Default: softmax")
    parser_domain.add_argument("--tiara_minimum_length", type=int, default=3000, help="Tiara | Minimum contig length. Anything lower than 3000 is not recommended. [Default: 3000]")
    parser_domain.add_argument("--tiara_options", type=str, default="", help="Tiara | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/ibe-uw/tiara")


    # Options
    opts = parser.parse_args()

    # Directories
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Threads
    if opts.n_jobs == -1:
        from multiprocessing import cpu_count 
        opts.n_jobs = cpu_count()
    assert opts.n_jobs >= 1, "--n_jobs must be ≥ 1.  To select all available threads, use -1."

    # Database
    if opts.veba_database is None:
        assert "VEBA_DATABASE" in os.environ, "Please set the following environment variable 'export VEBA_DATABASE=/path/to/veba_database' or provide path to --veba_database"
        opts.veba_database = os.environ["VEBA_DATABASE"]

    # Directories
    directories = dict()

    directories["project"] = create_directory(opts.project_directory)
    directories["sample"] = create_directory(os.path.join(directories["project"], opts.name))
    directories["output"] = create_directory(os.path.join(directories["sample"], "output"))
    directories["log"] = create_directory(os.path.join(directories["sample"], "log"))
    if not opts.tmpdir:
        opts.tmpdir = os.path.join(directories["sample"], "tmp")
    directories["tmp"] = create_directory(opts.tmpdir)
    directories["checkpoints"] = create_directory(os.path.join(directories["sample"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["sample"], "intermediate"))


    # Info
    print(format_header(__program__, "="), file=sys.stdout)
    print(format_header("Configuration:", "-"), file=sys.stdout)
    print(format_header("Name: {}".format(opts.name), "."), file=sys.stdout)
    print("Python version:", sys.version.replace("\n"," "), file=sys.stdout)
    print("Python path:", sys.executable, file=sys.stdout) #sys.path[2]
    print("Script version:", __version__, file=sys.stdout)
    print("VEBA Database:", opts.veba_database, file=sys.stdout)
    print("Moment:", get_timestamp(), file=sys.stdout)
    print("Directory:", os.getcwd(), file=sys.stdout)
    print("Commands:", list(filter(bool,sys.argv)),  sep="\n", file=sys.stdout)
    configure_parameters(opts, directories)
    sys.stdout.flush()

    # opts.bam = " ".join(opts.bam)


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
