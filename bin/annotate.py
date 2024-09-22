#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob,  shutil, gzip
from collections import OrderedDict, defaultdict

import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Soothsayer Ecosystem
from genopype import *
from genopype import __version__ as genopype_version
from soothsayer_utils import *
from soothsayer_utils.soothsayer_utils import assert_acceptable_arguments

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.9.21"

DIAMOND_HEADER_FIELDS = "qseqid sseqid stitle pident evalue bitscore qcovhsp scovhsp"


def get_preprocess_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, program):
    cmd = [
        "cat",
        " ".join(input_filepaths),
        "|",
        os.environ["seqkit"],
        "seq",
        "-M",
        opts.maximum_protein_length,
        "--only-id",
        ">",
        os.path.join(directories["tmp"], "proteins.faa"),
    ]
    return cmd

def get_diamond_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, program):
    tmp = os.path.join(directories["tmp"], program)
    # Command
    cmd = [

        "mkdir -p {}".format(tmp),

            "&&",

        os.environ["diamond"],
        "blastp",
        "--db {}".format(input_filepaths[1]),
        "--query {}".format(input_filepaths[0]),
        "--threads {}".format(opts.n_jobs),
        "-f 6 {}".format(DIAMOND_HEADER_FIELDS),
        "--evalue {}".format(opts.diamond_evalue),
        "-o {}".format(os.path.join(output_directory, "output.tsv")),
        "--max-target-seqs 1",
        "--tmpdir {}".format(tmp),
    ]
    if bool(opts.diamond_sensitivity):
        cmd += [ 
            "--{}".format(opts.diamond_sensitivity),
        ]
    cmd += [ 
        opts.diamond_options,
    ]

    cmd += [ 
            "&&",
            
        "pigz",
        "-f",
        "-p {}".format(opts.n_jobs),
        os.path.join(output_directory, "output.tsv"),

            "&&",

        "rm -rf {}".format(tmp),
    ]           

    return cmd


# HMMER
# def get_hmmsearch_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

#     # Command
#     cmd = [
#         "(",
#         os.environ["hmmsearch"],
#         "--tblout {}".format(os.path.join(output_directory, "output.tsv")),
#         "--cut_ga",
#         "--cpu {}".format(opts.n_jobs),
#         "--seed {}".format(opts.random_state + 1),
#         input_filepaths[1],
#         input_filepaths[0],
#         ">",
#         "/dev/null",
#         ")",

#             "&&",

#         "pigz",
#         "-f",
#         "-p {}".format(opts.n_jobs),
#         os.path.join(output_directory, "output.tsv"),
        
#     ]    
#     return cmd

# PyHMMSearch
def get_pyhmmsearch_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [
        os.environ["pyhmmsearch.py"],
        "-i {}".format(input_filepaths[0]),
        "-d {}".format(input_filepaths[1]),
        "-m gathering",
        "--n_jobs {}".format(opts.n_jobs),
        "|",
        os.environ["reformat_pyhmmsearch.py"],
        "-o {}".format(os.path.join(output_directory, "output.tsv.gz")),

    ]    
    return cmd

# # KofamScan
# def get_kofamscan_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

#     # Command
#     cmd = [
#         "(",


#         os.environ["exec_annotation"],
#         "--cpu {}".format(opts.n_jobs),
#         "-f detail-tsv",
#         "-p {}".format(input_filepaths[1]),
#         "-k {}".format(input_filepaths[2]),
#         "--tmp-dir {}".format(os.path.join(directories["tmp"], "kofamscan")),
#         opts.kofamscan_options,
#         input_filepaths[0],
#         "|",
#         'grep "*"',
#         ">",
#         os.path.join(output_directory, "output.tsv"),
#         ")",

#             "&&",

#         "pigz",
#         "-f",
#         "-p {}".format(opts.n_jobs),
#         os.path.join(output_directory, "output.tsv"),
#     ]    
#     return cmd

# PyKofamSearch
def get_pykofamsearch_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [
        
        os.environ["pykofamsearch.py"],
        "-i {}".format(input_filepaths[0]),
        "-d {}".format(input_filepaths[1]),
        "--n_jobs {}".format(opts.n_jobs),
        "|",
        os.environ["reformat_pykofamsearch.py"],
        "-o {}".format(os.path.join(output_directory, "output.tsv.gz")),

    ]    
    return cmd

def get_merge_annotations_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [
        os.environ["merge_annotations.py"],
        "--diamond_uniref {}".format(input_filepaths[0]),
        "--diamond_mibig {}".format(input_filepaths[1]),
        "--diamond_vfdb {}".format(input_filepaths[2]),
        "--diamond_cazy {}".format(input_filepaths[3]),
        "--pyhmmsearch_pfam {}".format(input_filepaths[4]),
        "--pyhmmsearch_amr {}".format(input_filepaths[5]),
        "--pyhmmsearch_antifam {}".format(input_filepaths[6]),
        "--pykofamsearch {}".format(input_filepaths[7]),
        "--pfam_clans {}".format(input_filepaths[8]),
        # "--veba_database {}".format(opts.veba_database),
        "--fasta {}".format(opts.proteins),
        "-o {}".format(output_directory),
        '--composite_name_joiner="{}"'.format(opts.composite_name_joiner),
    ]
    if opts.identifier_mapping:
        cmd += [ 
            "-i {}".format(opts.identifier_mapping),
            # "--genome_cluster_column_label {}".format(opts.genome_cluster_column_label),
            # "--protein_cluster_column_label {}".format(opts.protein_cluster_column_label),
        ]

    return cmd


def get_mcr_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [
        
        # Genomes
        os.environ["compile_ko_from_annotations.py"],
        "-i {}".format(input_filepaths[0]),
        "-o {}".format(os.path.join(output_directory, "kos.genomes.tsv")),
        "-l genome",

            "&&",
        
        os.environ["profile-pathway-coverage.py"],
        "-i {}".format(os.path.join(output_directory, "kos.genomes.tsv")),
        "-o {}".format(output_directory),
        "-d {}".format(os.path.join(opts.veba_database, "Annotate", "KEGG-Pathway-Profiler", "database.pkl.gz")),
        "--index_name id_genome",
        
            "&&",
            
        "mv",
        os.path.join(output_directory, "pathway_coverage.tsv.gz"),
        os.path.join(output_directory, "module_completion_ratios.genomes.tsv.gz"),
        
            "&&",
            
        "mv",
        os.path.join(output_directory, "pathway_output.pkl.gz"),
        os.path.join(output_directory, "module_completion_ratios.genomes.pkl.gz"),

            "&&",
            
        # Genome Clusters
        os.environ["compile_ko_from_annotations.py"],
        "-i {}".format(input_filepaths[0]),
        "-o {}".format(os.path.join(output_directory, "kos.genome_clusters.tsv")),
        "-l genome_cluster",

            "&&",
            
        os.environ["profile-pathway-coverage.py"],
        "-i {}".format(os.path.join(output_directory, "kos.genome_clusters.tsv")),
        "-o {}".format(output_directory),
        "-d {}".format(os.path.join(opts.veba_database, "Annotate", "KEGG-Pathway-Profiler", "database.pkl.gz")),
        "--index_name id_genome_cluster",
        
            "&&",
            
        "mv",
        os.path.join(output_directory, "pathway_coverage.tsv.gz"),
        os.path.join(output_directory, "module_completion_ratios.genome_clusters.tsv.gz"),
        
            "&&",
            
        "mv",
        os.path.join(output_directory, "pathway_output.pkl.gz"),
        os.path.join(output_directory, "module_completion_ratios.genome_clusters.pkl.gz"),
    ]
        
    return cmd

# def get_propagate_annotations_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

#     # Command
#     cmd = [
#         os.environ["propagate_annotations_from_representatives.py"],
#         "-i {}".format(input_filepaths[0]),
#         "-c {}".format(opts.protein_clusters),
#         "-o {}".format(output_filepaths[0]),
#     ]

#     return cmd

# ============
# Run Pipeline
# ============
# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = {
                "merge_annotations.py",
                # "propagate_annotations_from_representatives.py",
                "compile_ko_from_annotations.py",
                # "module_completion_ratios.py",
                "merge_annotations.py",
                "compile_custom_humann_database_from_annotations.py",
                }

    required_executables={
                # 1
                "diamond",
                # 2 
                "pyhmmsearch.py",
                "reformat_pyhmmsearch.py",
                # 3
                "pykofamsearch.py",
                "reformat_pykofamsearch.py",
                "profile-pathway-coverage.py",
                

     } | accessory_scripts

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
        executables[name] = "'{}'".format(os.path.join(opts.script_directory, "scripts", name)) # Can handle spaces in path

    print(format_header( "Adding executables to path from the following source: {}".format(opts.path_config), "-"), file=sys.stdout)
    for name, executable in executables.items():
        if name in required_executables:
            print(name, executable, sep = " --> ", file=sys.stdout)
            os.environ[name] = executable.strip()
    print("", file=sys.stdout)

# Pipeline
def create_pipeline(opts, directories, f_cmds):

    # .................................................................
    # Primordial
    # .................................................................
    # Commands file
    pipeline = ExecutablePipeline(name=__program__,  f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])


    # # ==========
    # # Diamond [nr]
    # # ==========
    # step = 0

    # program = "preprocess"
    # program_label = "{}__{}".format(step, program)
    # # Add to directories

    # # Info
    # description = "Preprocess sequences"
    # if os.path.isdir(opts.proteins):
    #     input_filepaths = glob.glob(os.path.join(opts.proteins, "*.{}".format(opts.extension)))
    # else:
    #     input_filepaths = [opts.proteins]
    # opts.protein_files = input_filepaths

    # # i/o

    # output_filepaths = [os.path.join(directories["tmp"], "proteins.faa")]

    # params = {
    #     "input_filepaths":input_filepaths,
    #     "output_filepaths":output_filepaths,
    #     "output_directory":output_directory,
    #     "opts":opts,
    #     "directories":directories,
    #     "program":program,
    # }

    # cmd = get_preprocess_cmd(**params)

    # pipeline.add_step(
    #             id=program,
    #             description = description,
    #             step=step,
    #             cmd=cmd,
    #             input_filepaths = input_filepaths,
    #             output_filepaths = output_filepaths,
    #             validate_inputs=True,
    #             validate_outputs=True,
    # )

    # ==========
    # Diamond [nr]
    # ==========
    step = 1

    program = "diamond-uniref"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Diamond [UniRef]"

    # i/o
    input_filepaths = [
        opts.proteins, 
        os.path.join(opts.veba_database, "Annotate", "UniRef", "{}.dmnd".format(opts.uniref)),
        ]
    output_filenames = ["output.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        "program":program,
    }

    cmd = get_diamond_cmd(**params)

    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
    )

    # ==========
    # Diamond
    # ==========
    step = 2

    program = "diamond-mibig"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Diamond [MIBiG]"

    # i/o
    input_filepaths = [
        opts.proteins, 
         os.path.join(opts.veba_database, "Annotate", "MIBiG", "mibig_v3.1.dmnd"),
        ]
    output_filenames = ["output.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        "program":program,
    }

    cmd = get_diamond_cmd(**params)

    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
    )

    # ==========
    # Diamond
    # ==========
    step = 3

    program = "diamond-vfdb"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Diamond [VFDB]"

    # i/o
    input_filepaths = [
        opts.proteins, 
         os.path.join(opts.veba_database, "Annotate", "VFDB", "VFDB_setA_pro.dmnd"),
        ]
    output_filenames = ["output.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        "program":program,
    }

    cmd = get_diamond_cmd(**params)

    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
    )

    # ==========
    # Diamond
    # ==========
    step = 4

    program = "diamond-cazy"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Diamond [CAZy]"

    # i/o
    input_filepaths = [
        opts.proteins, 
         os.path.join(opts.veba_database, "Annotate", "CAZy", "CAZyDB.07262023.dmnd"),
        ]
    output_filenames = ["output.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        "program":program,
    }

    cmd = get_diamond_cmd(**params)

    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
    )

    # ==========
    # HMMSearch-Pfam
    # ==========
    step = 5

    program = "pyhmmsearch-pfam"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "PyHMMSearch [Pfam]"

    # i/o
    input_filepaths = [
        opts.proteins, 
        os.path.join(opts.veba_database, "Annotate", "Pfam", "Pfam-A.hmm.gz"),
        ]
    output_filenames = ["output.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_pyhmmsearch_cmd(**params)

    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
    )

    # =============
    # HMMSearch-AMR
    # =============
    step = 6

    program = "pyhmmsearch-amr"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "PyHMMSearch [NCBIfam-AMR]"

    # i/o
    input_filepaths = [
        opts.proteins, 
        os.path.join(opts.veba_database, "Annotate", "NCBIfam-AMRFinder", "NCBIfam-AMRFinder.hmm.gz"),
        ]
    output_filenames = ["output.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_pyhmmsearch_cmd(**params)

    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
    )

    # =================
    # HMMSearch-AntiFam
    # =================
    step = 7

    program = "pyhmmsearch-antifam"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "PyHMMSearch [AntiFam]"

    # i/o
    input_filepaths = [
        opts.proteins, 
        os.path.join(opts.veba_database, "Contamination", "AntiFam", "AntiFam.hmm.gz"),
        ]
    output_filenames = ["output.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_pyhmmsearch_cmd(**params)

    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
    )
    # ==========
    # KOFAMSCAN
    # ==========
    step = 8

    program = "pykofamsearch"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "PyKofamSearch [KOfam]"
    # i/o
    input_filepaths = [
        opts.proteins, 
        os.path.join(opts.veba_database, "Annotate", "KOfam"),
        ]
    output_filenames = ["output.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_pykofamsearch_cmd(**params)

    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
    )

    # ==========
    # Merge annotations
    # ==========
    step = 9

    program = "merge_annotations"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories["output"]

    # Info
    description = "Merging annotation results"

    # i/o
    input_filepaths = [
        os.path.join(directories[("intermediate",  "1__diamond-uniref")], "output.tsv.gz"),
        os.path.join(directories[("intermediate",  "2__diamond-mibig")], "output.tsv.gz"),
        os.path.join(directories[("intermediate",  "3__diamond-vfdb")], "output.tsv.gz"),
        os.path.join(directories[("intermediate",  "4__diamond-cazy")], "output.tsv.gz"),
        os.path.join(directories[("intermediate",  "5__pyhmmsearch-pfam")], "output.tsv.gz"),
        os.path.join(directories[("intermediate",  "6__pyhmmsearch-amr")], "output.tsv.gz"),
        os.path.join(directories[("intermediate",  "7__pyhmmsearch-antifam")], "output.tsv.gz"),
        os.path.join(directories[("intermediate",  "8__pykofamsearch")], "output.tsv.gz"),
        os.path.join(opts.veba_database, "Annotate", "Pfam", "Pfam-A.clans.tsv.gz"),
    ]
    output_filenames = ["annotations.proteins.tsv.gz"]

    if opts.identifier_mapping:
        output_filenames += ["annotations.protein_clusters.tsv.gz"]

    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_merge_annotations_cmd(**params)

    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
    )

    # ==========
    # Propagate annotations
    # ==========
    if opts.identifier_mapping:
        step = 10

        program = "module_completion_ratios"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories["output"]

        # Info
        description = "Calculate module completion ratios"

        # i/o
        input_filepaths = [
            os.path.join(directories["output"], "annotations.proteins.tsv.gz"),
        ]
        output_filenames = ["kos.genomes.tsv", "kos.genome_clusters.tsv", "module_completion_ratios.genomes.tsv.gz", "module_completion_ratios.genome_clusters.tsv.gz"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_mcr_cmd(**params)

        pipeline.add_step(
                    id=program,
                    description = description,
                    step=step,
                    cmd=cmd,
                    input_filepaths = input_filepaths,
                    output_filepaths = output_filepaths,
                    validate_inputs=True,
                    validate_outputs=True,
        )

    return pipeline



# Configure parameters
def configure_parameters(opts, directories):

    assert_acceptable_arguments(opts.diamond_sensitivity, {"", "fast", "mid-sensitive", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"})
    assert_acceptable_arguments(opts.uniref, {"uniref90","uniref50"})

    # if opts.identifier_mapping or opts.protein_cluster_column_label:
        # assert bool(opts.identifier_mapping) == bool(opts.protein_cluster_column_label), "--identifier_mapping must be provided with --protein_cluster_column_label"

    assert opts.maximum_protein_length < 100000, "Must be < 100k because sequences this length or greater will cause HMMSearch to crash"
    # Set environment variables
    add_executables_to_environment(opts=opts)

def check_protein_sequence_lengths(path, maximum_protein_length):
    if path.endswith(".gz"):
        f = gzip.open(path, "rt")
    else:
        f = open(path, "r")

    failed_identifiers = list()
    for header, seq in pv(SimpleFastaParser(f), "Checking sequence lengths"):
        n = len(seq)
        if n > maximum_protein_length:
            failed_identifiers.append(header.split(" ")[0])
    assert len(failed_identifiers) == 0, "Please remove sequences longer than {}\nYou can do this by running `seqkit seq -M {} {}\n\nThe following identifiers are too long: {}".format(maximum_protein_length, maximum_protein_length, path, "\n\t".join(failed_identifiers))



def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -a <proteins> -o <output_directory> |Optional: -i <identifier_mapping]".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-a","--proteins", type=str, required=True, help = "path/to/proteins.faa fasta to annotate")
    parser_io.add_argument("-o","--output_directory", type=str, default="veba_output/annotation", help = "path/to/project_directory [Default: veba_output/annotation]")
    parser_io.add_argument("-i","--identifier_mapping", type=str, required=False,  help = "Tab-seperated value table (identifier_mapping.proteins.tsv created by cluster.py)  Requirements: 1) contain headers; 2) first column must be protein identifiers; and 3) contains these columns to the right in any order. Format: [id_protein]<tab>[organism_type]<tab>[id_genome]<tab>[sample_of_origin]<tab>[id_genome_cluster]<tab>[id_protein_cluster] with headers [Optional]")
    # parser_io.add_argument("-r","--representatives", type=str, required=False, help = "Tab-seperated value table of [id_protein]<tab>[id_protein_cluster].  Use this if you only want to annotate the representatives.")
    # parser_io.add_argument("--genome_column_label", type=str, default="id_genome", help = "--genome_column_label must be in --identifier_mapping header [Default: id_genome]")
    # parser_io.add_argument("--genome_cluster_column_label", type=str, default="id_genome_cluster", help = "--genome_cluster_column_label must be in --identifier_mapping header [Default: id_genome_cluster]")
    # parser_io.add_argument("--protein_cluster_column_label", type=str, default="id_protein_cluster", help = "--protein_cluster_column_label must be in --identifier_mapping header [Default: id_protein_cluster]")
    # parser_io.add_argument("--no_genome_clusters", action="store_true", help = "No genome clusters in --identifier_mapping header")
    # parser_io.add_argument("--no_protein_clusters", action="store_true", help = "No protein clusters in --identifier_mapping header")
    parser_io.add_argument("-x", "--extension", type=str, default="faa", help = "Fasta file extension for proteins if a directory is provided for --proteins (Can be gzipped) [Default: faa]")
    parser_io.add_argument("-M", "--maximum_protein_length", type=int, default=99999, help = "Proteins ≥ 100k will cause HMMSearch to crash and will use more resources [Default: 99999]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("--keep_temporary_directory", action="store_true",  help = "Keep temporary directory [Default is to remove]")
    parser_utility.add_argument("--no_check_protein_lengths", action="store_true",  help = "Do not check protein sequence lengths.  Not recommended.  Sequences must be < 100k or else PyHMMSearch will fail.")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Databases
    parser_databases = parser.add_argument_group('Database arguments')
    parser_databases.add_argument("-u", "--uniref", type=str,  default="uniref50", help="UniRef database to use {uniref90, uniref50}.  uniref90 receommended for well-characterized systems and uniref50 for less characterized systems [Default: uniref50]")
    parser_databases.add_argument("--veba_database", type=str,  help=f"VEBA database location.  [Default: $VEBA_DATABASE environment variable]")
    
    # Diamond
    parser_diamond = parser.add_argument_group('Diamond arguments')
    parser_diamond.add_argument("--diamond_sensitivity", type=str, default="", help="Diamond | Sensitivity [Default:  '']")
    parser_diamond.add_argument("--diamond_evalue", type=float, default=0.001, help="Diamond | E-Value [Default: 0.001]")
    parser_diamond.add_argument("--diamond_options", type=str, default="", help="Diamond | More options (e.g. --arg 1 ) [Default: '']")

    # Composite
    parser_composite = parser.add_argument_group('Composite arguments')
    parser_composite.add_argument("-j", "--composite_name_joiner", type=str, required=False,  default=";", help = "Composite label separator [Default: ; ]")

    # Options
    opts = parser.parse_args()
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
    directories["project"] = create_directory(opts.output_directory)
    directories["output"] = create_directory(os.path.join(directories["project"], "output"))
    directories["log"] = create_directory(os.path.join(directories["project"], "log"))
    directories["tmp"] = create_directory(os.path.join(directories["project"], "tmp"))
    directories["checkpoints"] = create_directory(os.path.join(directories["project"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["project"], "intermediate"))
    os.environ["TMPDIR"] = directories["tmp"]
    
    # Info
    print(format_header(__program__, "="), file=sys.stdout)
    print(format_header("Configuration:", "-"), file=sys.stdout)
    print("Python version:", sys.version.replace("\n"," "), file=sys.stdout)
    print("Python path:", sys.executable, file=sys.stdout) #sys.path[2]
    print("GenoPype version:", genopype_version, file=sys.stdout) #sys.path[2]
    print("Script version:", __version__, file=sys.stdout)
    print("VEBA Database:", opts.veba_database, file=sys.stdout)
    print("Moment:", get_timestamp(), file=sys.stdout)
    print("Directory:", os.getcwd(), file=sys.stdout)
    if "TMPDIR" in os.environ: print(os.environ["TMPDIR"], file=sys.stdout)
    print("Commands:", list(filter(bool,sys.argv)),  sep="\n", file=sys.stdout)
    configure_parameters(opts, directories)

    # Check sequence lengths
    if not opts.no_check_protein_lengths:
        check_protein_sequence_lengths(opts.proteins, opts.maximum_protein_length)

    sys.stdout.flush()

    # Run pipeline
    with open(os.path.join(directories["project"], "commands.sh"), "w") as f_cmds:
        pipeline = create_pipeline(
                     opts=opts,
                     directories=directories,
                     f_cmds=f_cmds,
        )
        pipeline.compile()
        pipeline.execute(restart_from_checkpoint=opts.restart_from_checkpoint)

    if not opts.keep_temporary_directory: # Add this to genopype
        for fp in glob.glob(os.path.join(directories["tmp"],"*")):
            if os.path.isdir(fp):
                shutil.rmtree(fp, ignore_errors=True)
            else:
                os.remove(fp)

if __name__ == "__main__":
    main()
