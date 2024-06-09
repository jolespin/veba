#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, shutil
from collections import OrderedDict, defaultdict

import pandas as pd
import numpy as np

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.11.07"



# Set up directories
def get_partition_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [
        
        # Make temporary directories
        "mkdir -p {}".format(os.path.join(output_directory,"genomes")),
            "&&",
        "mkdir -p {}".format(os.path.join(output_directory,"proteins")),
            "&&",
        "mkdir -p {}".format(os.path.join(output_directory,"proteins", "cpr_bacteria")),
            "&&",
        "mkdir -p {}".format(os.path.join(output_directory, "proteins", "non-cpr_prokaryotes")),
            "&&",

        # Merge scaffolds_to_bins.tsv
        "cat",
        " ".join(opts.scaffolds_to_bins), # os.path.join(directories["intermediate"], "*__checkm", "filtered", "scaffolds_to_bins.tsv"),
        ">",
        os.path.join(output_directory,"scaffolds_to_bins.tsv"),

            "&&",

        # Partition gene models
        os.environ["partition_gene_models.py"],
        "-i {}".format(os.path.join(output_directory,"scaffolds_to_bins.tsv")),
        "-f {}".format(opts.fasta),
        "-g {}".format(opts.gff), 
        "-d {}".format(opts.cds), 
        "-a {}".format(opts.protein), 
        "-o {}".format(os.path.join(output_directory, "genomes")),
        "--use_mag_as_description",
    ]
    return cmd

def get_gtdbtk_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

        # Run GTDB-Tk
        cmd = [
        "export GTDBTK_DATA_PATH={}".format(opts.gtdbtk_database),

            "&&",

        os.environ["gtdbtk"],
        "classify_wf",
        "--genome_dir {}".format(os.path.join(directories[("intermediate", "1__partition")],"genomes")),
        "--pplacer_cpus {}".format(opts.pplacer_threads),
        "--out_dir {}".format(output_directory),
        "-x fa",
        "--cpus {}".format(opts.n_jobs),
        "--tmpdir {}".format(opts.tmpdir),
        opts.gtdbtk_options,

            "&&",

        os.environ["concatenate_dataframes.py"],
        "--axis 0",
        "--allow_empty_or_missing_files",
        os.path.join(output_directory,  "classify", "gtdbtk.ar122.summary.tsv"),
        os.path.join(output_directory,  "classify", "gtdbtk.bac120.summary.tsv"),
        ">",
        os.path.join(output_directory,   "gtdbtk.summary.tsv"),
        ]

        return cmd

def get_split_cpr_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [
        "(",
            # CPR and bacterial genomes
            "cat",
            os.path.join(directories[("intermediate", "2__gtdbtk")], "classify", "gtdbtk.bac120.summary.tsv"),
            "|",
            "grep 'Patescibacteria'", 
            "|",
            "cut -f1",
            ">",
            os.path.join(directories[("intermediate", "1__partition")], "genomes", "cpr_bacteria.list"),

                "&&",

            # Non-CPR genomes
            "cat",
            os.path.join(directories[("intermediate", "2__gtdbtk")],  "classify", "gtdbtk.bac120.summary.tsv"),
            "|",
            "grep -v 'Patescibacteria'",
            "|",
            "cut -f1",
            "|",
            "tail -n +2",
            ">",
            os.path.join(directories[("intermediate", "1__partition")], "genomes", "non-cpr_bacteria.list"),
        ")",
        "2> /dev/null",
        "||",
        "(",
            "echo 'No bacterial genomes'",
            "&&",
            ">",
            os.path.join(directories[("intermediate", "1__partition")], "genomes", "cpr_bacteria.list"),
            "&&",
            ">",
            os.path.join(directories[("intermediate", "1__partition")], "genomes", "non-cpr_bacteria.list"),
        ")",
        "&&",

        # Archaea
        "(",
            "cat",
            os.path.join(directories[("intermediate", "2__gtdbtk")], "classify", "gtdbtk.ar122.summary.tsv"),
            "|",
            "cut -f1",
            "|",
            "tail -n +2",
            ">",
            os.path.join(directories[("intermediate", "1__partition")],  "genomes", "archaea.list"),
            ")",
            "2> /dev/null",
            "||",
            "(",
            "echo 'No archael genomes'",
            "&&",
            ">",
            os.path.join(directories[("intermediate", "1__partition")],  "genomes", "archaea.list"),
        ")",
        "&&",

        # Organize non-cpr genomes
        "for ID in $(cat %s %s); do cp %s %s; done"%( 
            os.path.join(directories[("intermediate", "1__partition")], "genomes", "non-cpr_bacteria.list"),
            os.path.join(directories[("intermediate", "1__partition")], "genomes", "archaea.list"),
            os.path.join(directories[("intermediate", "1__partition")], "genomes", "${ID}.faa"),
            os.path.join(directories[("intermediate", "1__partition")], "proteins", "non-cpr_prokaryotes"),
        ),
        "&&",

        # Organize CPR genomes
        "for ID in $(cat %s); do cp %s %s; done"%( 
            os.path.join(directories[("intermediate", "1__partition")], "genomes", "cpr_bacteria.list"),
            os.path.join(directories[("intermediate", "1__partition")], "genomes", "${ID}.faa"),
            os.path.join(directories[("intermediate", "1__partition")], "proteins", "cpr_bacteria"),
        ),
    ]
    return cmd


def get_checkm_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [
       # CheckM directories
        # "mkdir -p {}".format(output_directory),
        # "&&",
        "mkdir -p {}".format(os.path.join(directories[("intermediate", "1__partition")], "non-cpr_prokaryotes")),
        "&&",
        "mkdir -p {}".format(os.path.join(directories[("intermediate", "1__partition")],  "cpr_bacteria")),
        "&&",

        "(",
            "(",
                os.environ["checkm"],
                "lineage_wf",
                {"full":"", "reduced":"-r"}[opts.checkm_tree],
                "-g", 
                "--tab_table",
                "-f {}".format(os.path.join(output_directory,  "non-cpr_prokaryotes", "output.tsv")),
                "--pplacer_threads {}".format(opts.pplacer_threads),
                "-t {}".format(opts.n_jobs),
                "-x faa",
                "--tmpdir {}".format(opts.tmpdir), # Hack around: OSError: AF_UNIX path too long
                # "--tmpdir {}".format(os.path.join(directories["tmp"], "checkm")),
                os.path.join(directories[("intermediate", "1__partition")], "proteins", "non-cpr_prokaryotes"),
                os.path.join(output_directory,"non-cpr_prokaryotes"),
                    "&&",
                "rm -rf {} {}".format(
                    os.path.join(output_directory, "non-cpr_prokaryotes", "bins"),
                    os.path.join(output_directory, "non-cpr_prokaryotes", "storage"),
                ),
            ")",
            # "2>/dev/null",
            "||",
            ">",
            os.path.join(output_directory, "non-cpr_prokaryotes", "output.tsv"),
            "&&",

            # CPR Bactiera
            "(",
                os.environ["checkm"],
                "analyze",
                "-g",
                "-x faa",
                "-t {}".format(opts.n_jobs),
                "--tmpdir {}".format(opts.tmpdir), # Hack around: OSError: AF_UNIX path too long
                opts.cpr_database,
                os.path.join(directories[("intermediate", "1__partition")], "proteins", "cpr_bacteria"),
                os.path.join(output_directory, "cpr_bacteria"),

                    "&&",

                os.environ["checkm"],
                "qa",
                "--tab_table",
                "-f {}".format(os.path.join(output_directory, "cpr_bacteria", "output.tsv")),
                opts.cpr_database,
                os.path.join(output_directory, "cpr_bacteria"),

                    "&&",

                "rm -rf {} {}".format(
                    os.path.join(output_directory, "cpr_bacteria", "bins"),
                    os.path.join(output_directory, "cpr_bacteria", "storage"),
                ),
            ")",
            "||",
            ">",
            os.path.join(output_directory, "cpr_bacteria", "output.tsv"),
        ")",

            "&&",

        os.environ["concatenate_dataframes.py"],
        "--axis 0",
        "--allow_empty_or_missing_files",
        os.path.join(output_directory,  "non-cpr_prokaryotes", "output.tsv"),
        os.path.join(output_directory, "cpr_bacteria", "output.tsv"),
        ">",
        os.path.join(output_directory,  "output.concatenated.tsv"),

            "&&",
        
        os.environ["filter_checkm_results.py"],
        "-i {}".format(os.path.join(output_directory,  "output.concatenated.tsv")),
        "-b {}".format(os.path.join(directories[("intermediate", "1__partition")], "genomes")),
        "-o {}".format(directories["output"]),
        "-f {}".format(opts.fasta),
        "-m {}".format(opts.minimum_contig_length),
        "--unbinned",
        "--completeness {}".format(opts.checkm_completeness),
        "--contamination {}".format(opts.checkm_contamination),
        {True:"--strain_heterogeneity {}".format(opts.checkm_strain_heterogeneity),False:""}[bool(opts.checkm_strain_heterogeneity)],
        "-x fa",

    ]
    return cmd
    



def get_output_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [
        # Partion gene models
        os.environ["partition_gene_models.py"],
        "-i {}".format(os.path.join(output_directory, "scaffolds_to_bins.tsv")),
        # "-f {}".format(opts.fasta),
        "-g {}".format(opts.gff), 
        "-d {}".format(opts.cds), 
        "-a {}".format(opts.protein), 
        "-o {}".format(os.path.join(output_directory, "genomes")),
        "--use_mag_as_description",

            "&&",

        # Subset the GTDB-Tk output
        os.environ["subset_table.py"],
        "-t {}".format(os.path.join(directories[("intermediate", "2__gtdbtk")], "gtdbtk.summary.tsv")),
        "-i {}".format(os.path.join(output_directory,  "bins.list")),
        "--axis 0",
        ">",
        os.path.join(output_directory, "gtdbtk_output.filtered.tsv"),


        # "&&",
        # "rm -rf {}".format(
        #     os.path.join(directories["tmp"],"genomes"),
        #     os.path.join(directories["tmp"],"proteins"),
        # ),

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
    pipeline = ExecutablePipeline(name=__program__, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])

    # ==========
    # Partition
    # ==========
    step = 1

    # Info
    program = "partition"
    program_label = "{}__{}".format(step, program)
    description = "Create subdirectories and partition genomes"

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # i/o
    input_filepaths = [
        opts.fasta,
        opts.gff,
        opts.cds,
        opts.protein,
        *opts.scaffolds_to_bins,
    ]

    output_filenames =  [
        "scaffolds_to_bins.tsv", 
        "genomes/*.fa",
    ]

    output_filepaths = list(map(lambda fn:os.path.join(output_directory, fn), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_partition_cmd(**params)
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
                acceptable_returncodes={0},                    
                log_prefix=program_label,
    )

    # ==========
    # GTDB-Tk
    # ==========
    step = 2

    # Info
    program = "gtdbtk"
    program_label = "{}__{}".format(step, program)
    description = "Classify genomes and determine CPR taxa"

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # i/o
    input_filepaths = [
        "genomes/*.fa",
    ]

    output_filenames =  [
        "gtdbtk.summary.tsv",
    ]

    output_filepaths = list(map(lambda fn:os.path.join(output_directory, fn), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_gtdbtk_cmd(**params)
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
                acceptable_returncodes={0},                    
                log_prefix=program_label,
    )

    # ==========
    # Split CPR
    # ==========
    step = 3

    # Info
    program = "split_cpr"
    program_label = "{}__{}".format(step, program)
    description = "Split CPR bacteria from other prokaryotes"

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # i/o
    input_filepaths = output_filepaths

    output_filenames =  [
    ]

    output_filepaths = list(map(lambda fn:os.path.join(output_directory, fn), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_split_cpr_cmd(**params)
    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=False,
                errors_ok=False,
                acceptable_returncodes={0},                    
                log_prefix=program_label,
    )

    # ==========
    # CheckM
    # ==========
    step = 4

    # Info
    program = "checkm"
    program_label = "{}__{}".format(step, program)
    description = "Run CheckM CPR marker set for CPR bacteria and lineage_wf for other prokaryotes"

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # i/o
    input_filepaths = [
    ]

    output_filenames =  [
        "checkm_output.filtered.tsv",
        "bins.list",
        "binned.list",
        "scaffolds_to_bins.tsv",
    ]

    output_filepaths = list(map(lambda fn:os.path.join(directories["output"], fn), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_checkm_cmd(**params)
    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=False,
                validate_outputs=True,
                errors_ok=False,
                acceptable_returncodes={0},                    
                log_prefix=program_label,
    )


    # ==========
    # Output
    # ==========
    step = 5

    # Info
    program = "output"
    program_label = "{}__{}".format(step, program)
    description = "Output files"

    # Add to directories
    output_directory = directories["output"]
    # i/o
    input_filepaths = [
        os.path.join(directories[("intermediate", "2__gtdbtk")], "gtdbtk.summary.tsv"),
        os.path.join(directories["output"], "bins.list"),
        os.path.join(directories["output"], "scaffolds_to_bins.tsv"),
        opts.fasta,
        opts.gff,
        opts.cds,
        opts.protein,
    ]

    output_filenames =  [
        "gtdbtk_output.filtered.tsv",
    ]

    output_filepaths = list(map(lambda fn:os.path.join(output_directory, fn), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_output_cmd(**params)
    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=False,
                validate_outputs=True,
                errors_ok=False,
                acceptable_returncodes={0},                    
                log_prefix=program_label,
    )


    return pipeline



# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = {
        "partition_gene_models.py",
        "concatenate_dataframes.py",
        "subset_table.py",
        "filter_checkm_results.py",
    }

    required_executables={
        "gtdbtk",
        "checkm",
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
        if name.endswith(".py"):
            executables[name] = "python " + os.path.join(opts.script_directory, name)
        else: 
            executables[name] = os.path.join(opts.script_directory, name)


    print(format_header( "Adding executables to path from the following source: {}".format(opts.path_config), "-"), file=sys.stdout)
    for name, executable in executables.items():
        if name in required_executables:
            print(name, executable, sep = " --> ", file=sys.stdout)
            os.environ[name] = executable.strip()
    print("", file=sys.stdout)


# Configure parameters
def configure_parameters(opts, directories):
    updated_s2b = list()
    for fp in opts.scaffolds_to_bins:
        if  os.stat(fp).st_size == 0:
            print("[-] Empty file. Removing {}".format(fp), file=sys.stdout)
        else:
            print("[=] Keeping {}".format(fp), file=sys.stdout)
            updated_s2b.append(fp)
    opts.scaffolds_to_bins = updated_s2b


    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -g <gene_models.gff> -d <gene_models.ffn> -a <gene_models.faa> -f <scaffolds.fasta> -i <s2b_1.tsv s2b_2.tsv ... s2b_n.tsv >-o <output_directory>".format(__program__)

    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-i","--scaffolds_to_bins", type=str, nargs="+", required=True,  help = "One or more path/to/scaffolds_to_bins.tsv files separated by spaces, Format: [id_scaffold]<tab>[id_bin], No header")
    parser_io.add_argument("-f","--fasta", type=str, required=False,  help = "path/to/scaffolds.fasta")
    parser_io.add_argument("-g","--gff", type=str, required=True, help = "path/to/gene_models.gff")
    parser_io.add_argument("-d","--cds", type=str,  required=True, help = "path/to/gene_models.ffn")
    parser_io.add_argument("-a","--protein", type=str, required=False,  help = "path/to/gene_models.faa")
    parser_io.add_argument("-o","--output_directory", type=str, default="cpr_adjustment_output", help = "path/to/project_directory [Default: cpr_adjustment_output]")
    parser_io.add_argument("-m", "--minimum_contig_length", type=int, default=1500, help="Minimum contig length.  Anything under 2500 will default to 2500 for MetaBat2 [Default: 1500] ")
    parser_io.add_argument("-s", "--minimum_genome_length", type=int, default=150000, help="Minimum genome length.  [Default: 150000]")
    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    # parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
    parser_utility.add_argument("--tmpdir", type=str, help="Set temporary directory")  #site-packges in future
    parser_utility.add_argument("--remove_intermediate_directory", action="store_true", help="Remove intermediate directory")

    # CheckM
    parser_checkm = parser.add_argument_group('CheckM arguments')
    parser_checkm.add_argument("--pplacer_threads", type=int, default=1, help = "Number of threads used for pplacer. Multithreaded uses a lot memory so don't use unless you have the resources [Default: 1]")
    parser_checkm.add_argument("--checkm_tree", type=str, default="reduced", help="CheckM tree type either 'reduced' or 'full' [Default: reduced]")
    parser_checkm.add_argument("--checkm_completeness", type=float, default=50.0, help="CheckM completeness threshold [Default: 50]")
    parser_checkm.add_argument("--checkm_contamination", type=float, default=10.0, help="CheckM contamination threshold [Default: 10]")
    parser_checkm.add_argument("--checkm_strain_heterogeneity", type=float,  help="CheckM strain hetereogeneity threshold")
    parser_checkm.add_argument("--checkm_options", type=str, default="", help="CheckM lineage_wf | More options (e.g. --arg 1 ) [Default: '']")
    parser_checkm.add_argument("--cpr_database", type=str, required=True, help="CheckM | path/to/cpr.hmm (e.g. --arg 1 )")

    # GTDBTk
    parser_gtdbtk = parser.add_argument_group('GTDB-Tk arguments')
    parser_gtdbtk.add_argument("--gtdbtk_database", type=str, required=True, help="GTDB-Tk | path/to/gtdbtk_database (e.g. --arg 1 )")
    parser_gtdbtk.add_argument("--gtdbtk_options", type=str, default="", help="GTDB-Tk | classify_wf options (e.g. --arg 1 ) [Default: '']")

    # Options
    opts = parser.parse_args()

    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Directories
    directories = dict()
    directories["output"] = create_directory(opts.output_directory)
    directories["intermediate"] = create_directory(os.path.join(directories["output"], "intermediate"))
    directories["log"] = create_directory(os.path.join(directories["output"], "log"))
    if not opts.tmpdir:
        opts.tmpdir = create_directory(os.path.join(directories["output"],"intermediate", "tmp"))
    directories["tmp"] = create_directory(opts.tmpdir)
    directories["checkpoints"] = create_directory(os.path.join(directories["output"],"intermediate", "checkpoints"))
    os.environ["TMPDIR"] = directories["tmp"]

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

    if opts.remove_intermediate_directory:
        print("[-] Removing intermediate directory: {}".format(directories["intermediate"]), file=sys.stdout)
        shutil.rmtree(directories["intermediate"])
    
   


if __name__ == "__main__":
    main(sys.argv[1:])
