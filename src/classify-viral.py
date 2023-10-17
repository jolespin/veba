#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
from collections import OrderedDict, defaultdict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from genopype import __version__ as genopype_version
from soothsayer_utils import *


pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.10.16"


def get_concatenate_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [
        os.environ["concatenate_dataframes.py"],
        "-a 0",
        os.path.join(opts.viral_binning_directory, "*", "output", "checkv_results.filtered.tsv"),

        "|",

        os.environ["cut_table_by_column_labels.py"], #[id_mag]<tab>[id_slc]<tab>[classification]|OPTIONAL:<tab>[weight]
        "--columns lineage,agreement,taxid,provirus",
        ">",
        os.path.join(directories["output"], "taxonomy.tsv"),
    ]
    return cmd


def get_genomad_taxonomy_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [ 
        os.environ["scaffolds_to_bins.py"],
        "-g {}".format(opts.genomes),
        "-x {}".format(opts.extension),
        ">",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),

            "&&",

        os.environ["seqkit"],
        "seq",
        "$(cat {})".format(opts.genomes),
        ">",
        os.path.join(directories["tmp"], "concatenated_genomes.fa"),

            "&&",

        os.environ["genomad_taxonomy_wrapper.py"],
        "-i {}".format(os.path.join(output_directory, "scaffolds_to_bins.tsv")),
        "-f {}".format(os.path.join(directories["tmp"], "concatenated_genomes.fa")),
        "-o {}".format(output_directory),
        "--veba_database {}".format(os.path.join(opts.veba_database)),
        "--n_jobs {}".format(opts.n_jobs),

            "&&",

        os.environ["cut_table_by_column_labels.py"], #[id_mag]<tab>[id_slc]<tab>[classification]|OPTIONAL:<tab>[weight]
        "--columns lineage,agreement,taxid",
        os.path.join(output_directory,  "taxonomy.tsv"),
        ">",
        os.path.join(directories["output"], "taxonomy.tsv"),

    ]

    return cmd

def get_consensus_genome_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [ 
        "cat",
        input_filepaths[0],
        "|",
        os.environ["cut_table_by_column_labels.py"], #[id_mag]<tab>[id_slc]<tab>[classification]|OPTIONAL:<tab>[weight]
        "--columns lineage,agreement",
        "|",
        "tail -n +2",
        "|",
        os.environ["insert_column_to_table.py"],
        "-c {}".format(opts.clusters),
        "-n id_genome_cluster",
        "-i 0",
        "|",
        os.environ["consensus_genome_classification_unranked.py"],
        "-o {}".format(output_filepaths[0]),
        "-t {}".format(opts.threshold),
        "--unclassified_label 'Unclassified virus'",
        "--unclassified_taxonid -1",
        "--veba_database {}".format(opts.veba_database),
        "--verbose",
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
        "insert_column_to_table.py",
        "concatenate_dataframes.py",
        "consensus_genome_classification_unranked.py",
        "cut_table_by_column_labels.py",
        "genomad_taxonomy_wrapper.py",
        "scaffolds_to_bins.py",

    ])

    
    required_executables=set(
        [
            "genomad",
            "seqkit",
        ]
     ) | accessory_scripts

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

    step = 1

    if opts.viral_binning_directory:
        program = "concatenate"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories["output"]

        # Info
        description = "Concatenating taxonomy results from geNomad"

        # i/o
        input_filepaths = [
            opts.viral_binning_directory,
            ]
        output_filenames = ["taxonomy.tsv"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_concatenate_cmd(**params)

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

    if opts.genomes:
        program = "genomad"
        program_label = "{}__{}".format(step, program)

        # Add to directories
        output_directory = create_directory(os.path.join(directories["intermediate"], "genomad_taxonomy"))

        # Info
        description = "Classifying viral contigs with geNomad"

        # i/o
        input_filepaths = [
            opts.genomes,
            ]
        output_filenames = ["taxonomy.tsv"]
        output_filepaths = list(map(lambda filename: os.path.join(directories["output"], filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_genomad_taxonomy_cmd(**params)

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
            os.path.join(directories["output"], "taxonomy.tsv"),
            opts.clusters,
            ]
        output_filenames = ["taxonomy.clusters.tsv"]
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



    return pipeline

# Configure parameters
def configure_parameters(opts, directories):
    assert bool(opts.viral_binning_directory) != bool(opts.genomes), "Must choose either --viral_binning_directory or --genomes, not both."

    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <viral_binning_directory>|-g <genomes.list>  -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i","--viral_binning_directory", type=str, help = "Either: path/to/checkv/quality_summary.tsv or directory of veba_output/binning/viral")
    parser_io.add_argument("-g","--genomes", type=str, help = "path/to/genomes.list where each line is a path to a genome.fasta [Cannot be ued with --viral_binning_directory]")
    parser_io.add_argument("-c","--clusters", type=str, help = "path/to/clusters.tsv, Format: [id_mag]<tab>[id_cluster], No header.")
    parser_io.add_argument("-o","--output_directory", type=str, default="veba_output/classify/viral", help = "path/to/output_directory [Default: veba_output/classify/viral]")
    parser_io.add_argument("-x","--extension", type=str, default="fa", help = "path/to/output_directory.  Does not support gzipped. [Default: fa]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Databases
    parser_databases = parser.add_argument_group('Database arguments')
    parser_databases.add_argument("--veba_database", type=str,  help=f"VEBA database location.  [Default: $VEBA_DATABASE environment variable]")

    # Consensus genome classification
    parser_consensus = parser.add_argument_group('Consensus genome classification arguments')
    parser_consensus.add_argument("-t","--threshold", default=0.5, type=float, help = "Fraction of classifications for consensus [Default: 0.5]")

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
