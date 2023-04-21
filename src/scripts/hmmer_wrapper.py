#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
from collections import OrderedDict, defaultdict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.4.17"

# Check
def get_check_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command

    # Command
    cmd = [
        os.environ["check_fasta_duplicates.py"],
        opts.fasta,
        ]

    return cmd


# HMMSearch
def get_hmmer_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command

    # Command
    cmd = [
        os.environ[opts.algorithm],
        "--tblout {}".format(os.path.join(directories["intermediate"], "tblout.tsv")),
        "--domtblout {}".format(os.path.join(directories["intermediate"], "domtblout.tsv")),
        "--cpu {}".format(opts.n_jobs),
        "--seed {}".format(opts.random_state + 1),
        "--{}".format(opts.hmmsearch_threshold) if opts.hmmsearch_threshold == "e" else "-E {}".format(opts.hmmsearch_evalue),
        opts.database_hmm,
        opts.fasta,
        ">",
        "/dev/null",

            "&&",

        "gzip",
        os.path.join(directories["intermediate"], "*.tsv"),

    ]

    return cmd

# Compile
def get_compile_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    
    # Command
    cmd = [
        os.environ["filter_hmmsearch_results.py"],
        "-i {}".format(os.path.join(directories["intermediate"], "tblout.tsv.gz")),
        "-o {}".format(os.path.join(directories["output"], "output.tsv.gz")),
        "-f {}".format(opts.hmm_marker_field),
        "--synopsis {}".format(os.path.join(directories["output"], "synopsis.tsv.gz")),
        "--region {}".format(opts.region),
    ]

    if opts.scores_cutoff:
        cmd += [ 
            "-s {}".format(opts.scores_cutoff),
        ]
    
    cmd += [ 
            "&&",

        "gzip -d -c",
        os.path.join(directories["output"], "output.tsv.gz"),
        "|",
        "cut -f1",
        "|",
        "tail -n +3",
        ">",
        os.path.join(directories["output"], "identifiers.list"),

            "&&",

        "gzip -d -c",
        os.path.join(directories["output"], "output.tsv.gz"),
        "|",
        "cut -f3",
        "|",
        "tail -n +3",
        "|",
        "sort -u",
        ">",
        os.path.join(directories["output"], "hmms_detected.list"),

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
        "filter_hmmsearch_results.py",
        "check_fasta_duplicates.py",
    ])

    
    required_executables={
                # 1
                opts.algorithm,

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
        executables[name] = "python " + os.path.join(opts.script_directory, name)
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
    # HMMSearch
    # ==========
    step = 1

    program = "hmmsearch"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories["intermediate"] 
    # Info
    description = "HMMSearch for proteins"


    # i/o
    input_filepaths = [
        opts.fasta,
        ]
    output_filenames = ["tblout.tsv.gz", "domtblout.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_hmmer_cmd(**params)

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
    description = "Compile HMMER results"

    # i/o
    input_filepaths = output_filepaths
    output_filenames = [
        "output.tsv.gz",
        "identifiers.list",
        "hmms_detected.list",
        "synopsis.tsv.gz",
    ]
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
    if opts.hmmsearch_threshold.lower() == "e":
        opts.hmmsearch_threshold = None
    assert_acceptable_arguments(opts.algorithm, {"hmmsearch", "hmmscan"})
    
    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <fasta> -d <database_hmms>  -o <output_directory> |Optional: -s <scores_cutoff.tsv>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i","--fasta", type=str, required=True, help = "Protein fasta file")
    parser_io.add_argument("-d", "--database_hmm", type=str,  required=True, help=f"path/to/HMM database of markers")
    parser_io.add_argument("-o","--output_directory", type=str, default="hmmer_output", help = "path/to/project_directory [Default: hmmer_output]")
    parser_io.add_argument("-s","--scores_cutoff", type=str, help = "path/to/scores_cutoff.tsv. No header. [id_hmm]<tab>[score]")


    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Use -1 for completely random. Use 0 for consecutive random states.  Use any other positive integer for the same random state for all iterations [Default: 0]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # HMMER
    parser_hmmer = parser.add_argument_group('HMMER arguments')
    parser_hmmer.add_argument("-a", "--algorithm", type=str, default="hmmsearch", help="HMMER | {hmmsearch, hmmscan} [Default: hmmsearch]")
    parser_hmmer.add_argument("-f", "--hmm_marker_field", default="accession", type=str, help="HMM reference type (accession, name) [Default: accession]")
    parser_hmmer.add_argument("--region",  default="full_sequence", type=str, help="{full_sequence, best_domain} [Default: full_sequence]")
    parser_hmmer.add_argument("--hmmsearch_threshold", type=str, default="e", help="HMMER | Threshold {cut_ga, cut_nc, gut_tc, e} [Default:  e]")
    parser_hmmer.add_argument("--hmmsearch_evalue", type=float, default=10.0, help="HMMER | E-Value [Default: 10.0]")
    parser_hmmer.add_argument("--hmmsearch_options", type=str, default="", help="HMMER | More options (e.g. --arg 1 ) [Default: '']")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Threads
    if opts.n_jobs == -1:
        from multiprocessing import cpu_count 
        opts.n_jobs = cpu_count()
    assert opts.n_jobs >= 1, "--n_jobs must be ≥ 1.  To select all available threads, use -1."

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

    # Symlink genomes

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
    main()
