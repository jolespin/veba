#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, shutil, time
from multiprocessing import cpu_count
from collections import OrderedDict, defaultdict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.9.15"

# Check
def get_check_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command

    # Command
    cmd = [
        os.environ["check_fasta_duplicates.py"],
        opts.fasta,
        ]

    return cmd

def get_mmseqs2_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [
            os.environ["mmseqs"],
            opts.algorithm,
            opts.fasta,
            os.path.join(output_directory, "mmseqs2"),
            directories["tmp"],
            "--threads {}".format(opts.n_jobs),
            "--min-seq-id {}".format(opts.minimum_identity_threshold/100),
            "-c {}".format(opts.minimum_coverage_threshold),
            "--cov-mode 1",
            opts.mmseqs2_options,

                "&&",

            "gzip",
            os.path.join(output_directory, "mmseqs2_rep_seq.fasta"),

                "&&",

            "rm -rf",
            os.path.join(output_directory, "mmseqs2_all_seqs.fasta"),
            os.path.join(directories["tmp"], "*"),
        ]

    return cmd

# Compile
def get_compile_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    
    # Command
    cmd = [

            os.environ["edgelist_to_clusters.py"],
            "-i {}".format(input_filepaths[0]),
            "--no_singletons" if bool(opts.no_singletons) else "",
            "--cluster_prefix {}".format(opts.cluster_prefix) if bool(opts.cluster_prefix) else "",
            "--cluster_suffix {}".format(opts.cluster_suffix) if bool(opts.cluster_suffix) else "",
            "--cluster_prefix_zfill {}".format(opts.cluster_prefix_zfill),
            "-o {}".format(os.path.join(output_directory, "{}.tsv".format(opts.basename))),
            # "-g {}".format(os.path.join(output_directory, "{}.networkx_graph.pkl".format(opts.basename))),
            # "-d {}".format(os.path.join(output_directory, "{}.dict.pkl".format(opts.basename))),
            "--identifiers {}".format(opts.identifiers) if bool(opts.identifiers) else "",
            
                "&&",

            os.environ["reformat_representative_sequences.py"],
            "-c {}".format(os.path.join(output_directory, "{}.tsv".format(opts.basename))),
            "-i {}".format(input_filepaths[1]),
            "-f {}".format(opts.representative_output_format),
            "-o {}".format(output_filepaths[1]),
    ]

    if opts.no_sequences_and_header:
        cmd += [ 
            "--no_sequences",
            "--no_header",
        ]

    return cmd

# ============
# Run Pipeline
# ============
# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = set([ 
        "check_fasta_duplicates.py",
        "edgelist_to_clusters.py",
        "reformat_representative_sequences.py",
    ])

    
    required_executables={
                "mmseqs",

     } | accessory_scripts

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
    for name in sorted(accessory_scripts):
        executables[name] = "'{}'".format(os.path.join(opts.script_directory, name)) # Can handle spaces in path
        
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

    
    # ==========
    # Preprocessing
    # ==========
    
    program = "check"
    # Add to directories
    output_directory = directories["tmp"] 

    # Info
    step = 0
    description = "Check sequences for duplicates"

    # i/o
    input_filepaths = [opts.fasta]
    output_filepaths = [
    ]

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_check_cmd(**params)

    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=False,
    )

    # ==========
    # MMSEQS2
    # ==========
    
    program = "mmseqs2"
    # Add to directories
    output_directory = directories["intermediate"] 

    # Info
    step = 1
    description = "Cluster sequences via MMSEQS2"

    # i/o
    input_filepaths = [opts.fasta]
    output_filenames = [
        "mmseqs2_cluster.tsv",
        "mmseqs2_rep_seq.fasta.gz",
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_mmseqs2_cmd(**params)

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
    # Compile
    # ==========
    
    program = "compile"
    # Add to directories
    output_directory = directories["output"] 

    # Info
    step = 2
    description = "Compile MMSEQS2 results"

    # i/o
    input_filepaths = output_filepaths
    output_filenames = [
        "{}.tsv".format(opts.basename),
    ]
    if opts.representative_output_format == "table":
        output_filenames += ["representative_sequences.tsv.gz"]
    if opts.representative_output_format == "fasta":
        output_filenames += ["representative_sequences.fasta.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_compile_cmd(**params)

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

    assert_acceptable_arguments(opts.algorithm, {"easy-cluster", "easy-linclust"})
    assert_acceptable_arguments(opts.representative_output_format, {"table", "fasta"})
    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <sequences.fasta> -o <output_directory>".format(__program__)

    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i","--fasta", type=str, help = "Fasta file")
    parser_io.add_argument("-o","--output_directory", type=str, default="mmseqs2_output", help = "path/to/project_directory [Default: mmseqs2_output]")
    parser_io.add_argument("-e", "--no_singletons", action="store_true", help="Exclude singletons")
    parser_io.add_argument("-b", "--basename", type=str, default="clusters", help="Basename for clustering files [Default: clusters]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
    # parser_utility.add_argument("--verbose", action='store_true')

    # MMSEQS2
    parser_mmseqs2 = parser.add_argument_group('MMSEQS2 arguments')
    parser_mmseqs2.add_argument("-a", "--algorithm", type=str, default="easy-cluster", help="MMSEQS2 | {easy-cluster, easy-linclust} [Default: easy-cluster]")
    parser_mmseqs2.add_argument("-t", "--minimum_identity_threshold", type=float, default=50.0, help="MMSEQS2 | Percent identity threshold (Range (0.0, 100.0]) [Default: 50.0]")
    parser_mmseqs2.add_argument("-c", "--minimum_coverage_threshold", type=float, default=0.8, help="MMSEQS2 | Coverage threshold (Range (0.0, 1.0]) [Default: 0.8]")
    parser_mmseqs2.add_argument("--cluster_prefix", type=str, default="SC-", help="Sequence cluster prefix [Default: 'SC-]")
    parser_mmseqs2.add_argument("--cluster_suffix", type=str, default="", help="Sequence cluster suffix [Default: '']")
    parser_mmseqs2.add_argument("--cluster_prefix_zfill", type=int, default=0, help="Sequence cluster prefix zfill. Use 7 to match identifiers from OrthoFinder.  Use 0 to add no zfill. [Default: 0]") #7
    parser_mmseqs2.add_argument("--mmseqs2_options", type=str, default="", help="MMSEQS2 | More options (e.g. --arg 1 ) [Default: '']")
    parser_mmseqs2.add_argument("--identifiers", type=str, help = "Identifiers to include for `edgelist_to_clusters.py`.  If missing identifiers and singletons are allowed, then they will be included as singleton clusters with weight of np.inf")
    parser_mmseqs2.add_argument("--no_sequences_and_header", action="store_true", help = "Don't include sequences or header in table.  Useful for concatenation and reduced redundancy of sequences")
    parser_mmseqs2.add_argument("-f","--representative_output_format", type=str, default="fasta", help = "Format of output for representative sequences: {table, fasta} [Default: fasta]") # Should fasta be the new default?

    # Options
    opts = parser.parse_args()

    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    # Threads
    if opts.n_jobs == -1:
        opts.n_jobs = cpu_count()
    assert opts.n_jobs >= 1, "--n_jobs must be ≥ 1 (or -1 to use all available threads)"

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
    print("Script version:", __version__, file=sys.stdout)
    print("Moment:", get_timestamp(), file=sys.stdout)
    print("Directory:", os.getcwd(), file=sys.stdout)
    print("Commands:", list(filter(bool,sys.argv)),  sep="\n", file=sys.stdout)
    configure_parameters(opts, directories)
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

if __name__ == "__main__":
    main(sys.argv[1:])


