#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
from collections import OrderedDict, defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

# DATABASE_CHECKV="/usr/local/scratch/CORE/jespinoz/db/checkv/checkv-db-v1.0/"

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.7.8"


# Preprocess
def get_preprocess_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [ 
"""
python -c "import sys, glob, pandas as pd; pd.concat(map(lambda fp: pd.read_csv(fp, sep='\t', index_col=0), glob.glob('{}')), axis=0).to_csv('{}', sep='\t')"
""".format(
    os.path.join(opts.viral_binning_directory,"*", "output", "quality_summary.filtered.tsv"),
    os.path.join(directories["preprocessing"], "all_samples.quality_summary.filtered.tsv"),
    ),

    ]
    return cmd

def get_compile_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [ 
        os.environ["compile_viral_classifications.py"],
        "-i {}".format(input_filepaths[0]),
        "--checkv_database {}".format(opts.checkv_database),
        "-o {}".format(output_filepaths[0]),
    ]
    if opts.clusters:
        cmd += ["-c {}".format(opts.clusters)]

    return cmd

def get_consensus_genome_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [ 

        "cat",
        input_filepaths[0],
        "|",
        "cut -f1,2,3,4",
        "|",
        "tail -n +2",
        "|",
        os.environ["consensus_genome_classification.py"],
        "--leniency 1",
        "-o {}".format(output_filepaths[0]),
        "--simple",
    ]
    return cmd

def get_consensus_source_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [ 

        "cat",
        input_filepaths[0],
        "|",
        "cut -f1,2,6",
        "|",
        "tail -n +2",
        ">",
        os.path.join(directories["tmp"], "source.tsv"),
        "&&",
        os.environ["consensus_orthogroup_annotation.py"],
        "-i {}".format(os.path.join(directories["tmp"], "source.tsv")),
        "--similarity_threshold {}".format(opts.similarity_threshold),
        "--retain_unannotated {}".format(opts.retain_unannotated),
        "--unannotated_weight {}".format(opts.unannotated_weight),
        "--representative_threshold {}".format(opts.representative_threshold),
        "-o {}".format(os.path.join(output_directory, "unifunc_output")),
        "--orthogroup_label id_cluster",
        "&&",
        "mv {} {}".format(
            os.path.join(output_directory, "unifunc_output", "consensus_annotations.tsv"),
            os.path.join(output_directory,  "consensus_source_classification.tsv"),
        ),
        "&&",
        os.path.join(directories["tmp"], "source.tsv"),

    ]
    return cmd



# Symlink
# def get_symlink_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    
#     # Command
#     cmd = ["("]
#     for filepath in input_filepaths:
#         # cmd.append("ln -f -s {} {}".format(os.path.realpath(filepath), os.path.realpath(output_directory)))
#         cmd.append("ln -f -s {} {}".format(os.path.realpath(filepath), output_directory))
#         cmd.append("&&")
#     cmd[-1] = ")"
#     return cmd

# ============
# Run Pipeline
# ============
# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = set([ 
        "compile_viral_classifications.py",
        "consensus_genome_classification.py",
        "consensus_orthogroup_annotation.py",
    ])

    
    required_executables={
        "unifunc",
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
        executables[name] = "python " + os.path.join(opts.script_directory, "scripts", name)
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

    if os.path.isdir(opts.viral_binning_directory):
        # ==========
        # Preprocess
        # ==========
        step = 0

        program = "preprocess"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories["preprocessing"]

        # Info
        description = "Merge CheckV results"


        # i/o
        input_filepaths = [
            opts.viral_binning_directory,
            ]
        output_filenames = ["all_samples.quality_summary.filtered.tsv"]
        output_filepaths = list(map(lambda filename: os.path.join(directories["preprocessing"], filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_preprocess_cmd(**params)

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

        opts.viral_binning_directory = os.path.join(directories["preprocessing"], "all_samples.quality_summary.filtered.tsv")

  
    # ==========
    # Get classifications
    # ==========
    step = 1

    program = "compile"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories["output"] 

    # Info
    description = "Compiile CheckV classifications"


    # i/o
    input_filepaths = [
        opts.viral_binning_directory,
        ]
    if opts.clusters:
        input_filepaths += [opts.clusters]
    output_filenames = ["viral_taxonomy.tsv"]
    output_filepaths = list(map(lambda filename: os.path.join(directories["output"], filename), output_filenames))

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

    if opts.clusters:
        # ==========
        # consensus genome classification
        # ==========
        step = 2

        program = "consensus_genome_classification"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories["output"]# = create_directory(os.path.join(directories["intermediate"], program_label))

        # Info
        description = "Consensus genome classification"


        # i/o
        input_filepaths = [
            os.path.join(directories["output"], "viral_taxonomy.tsv"),
            opts.clusters,
            ]
        output_filenames = ["consensus_genome_classification.tsv"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_consensus_genome_cmd(**params)

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
        # consensus source classification
        # ==========
        step = 3

        program = "consensus_source_habitat"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories["output"]# = create_directory(os.path.join(directories["intermediate"], program_label))

        # Info
        description = "Consensus source habitat"


        # i/o
        input_filepaths = [
            os.path.join(directories["output"], "viral_taxonomy.tsv"),
            opts.clusters,
            ]
        output_filenames = ["consensus_source_classification.tsv"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_consensus_source_cmd(**params)

        pipeline.add_step(
                    id=program,
                    description = description,
                    step=step,
                    cmd=cmd,
                    input_filepaths = input_filepaths,
                    output_filepaths = output_filepaths,
                    validate_inputs=True,
                    validate_outputs=True,
                    acceptable_returncodes={0,126},
        )


    return pipeline

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
    usage = "{} -i <viral_binning_directory>  -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i","--viral_binning_directory", type=str, required=True, help = "Either: path/to/checkv/quality_summary.tsv or directory of veba_output/binning/viral")
    parser_io.add_argument("-c","--clusters", type=str, help = "path/to/clusters.tsv, Format: [id_mag]<tab>[id_cluster], No header.")
    parser_io.add_argument("-o","--output_directory", type=str, default="veba_output/classify/viral", help = "path/to/output_directory [Default: veba_output/classify/viral]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Databases
    parser_databases = parser.add_argument_group('Database arguments')
    parser_databases.add_argument("--veba_database", type=str,  help=f"VEBA database location.  [Default: $VEBA_DATABASE environment variable]")

    # # CheckV
    # parser_checkv = parser.add_argument_group('CheckV arguments')
    # parser_checkv.add_argument("--checkv_database", type=str, default=DATABASE_CHECKV, help="CheckV | path/to/gtdbtk_database (e.g. --arg 1 ) [Default: {}]".format(DATABASE_CHECKV))

    # Consensus genome classification
    parser_consensus_source = parser.add_argument_group('Consensus habitat/isolation source arguments')
    parser_consensus_source.add_argument("--similarity_threshold", type=float, default=0.8, help = "Threshold for similarity analysis [Default: 0.8]")
    parser_consensus_source.add_argument("--retain_unannotated", type=int, default=1, help = "Consider unannotations (i.e., blank functions) in the scording system [Default: 1]")
    parser_consensus_source.add_argument("--unannotated_weight", type=float, default=0.382, help = "Weight for unannotations (i.e., blank functions) in the scording system? [Default: 0.382]")
    parser_consensus_source.add_argument("--representative_threshold", type=float, default=0.618, help = "Score to consider as representative [Default: 0.618]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Database
    if opts.veba_database is None:
        assert "VEBA_DATABASE" in os.environ, "Please set the following environment variable 'export VEBA_DATABASE=/path/to/veba_database' or provide path to --veba_database"
        opts.veba_database = os.environ["VEBA_DATABASE"]
    opts.checkv_database = os.path.join(opts.veba_database, "Classify", "CheckV")

    # Directories
    directories = dict()
    directories["project"] = create_directory(opts.output_directory)
    directories["preprocessing"] = create_directory(os.path.join(directories["project"], "preprocessing"))
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
    print("VEBA Database:", opts.veba_database, file=sys.stdout)
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
    main()
