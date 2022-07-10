#!/usr/bin/env python
from __future__ import print_function, division
from pdb import Restart
import sys, os, argparse, glob
from collections import OrderedDict, defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

# DATABASE_GTDBTK="/usr/local/scratch/CORE/jespinoz/db/gtdbtk/release202/"

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.06.07"

# .............................................................................
# Notes
# .............................................................................

# GTDB-Tk
def get_concatenate_gtdbtk_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [ 
# """
# python -c "import sys, glob, pandas as pd; pd.concat(map(lambda fp: pd.read_csv(fp, sep='\t', index_col=0), glob.glob('{}')), axis=0).to_csv('{}', sep='\t')"
# """.format(
#     os.path.join(opts.prokaryotic_binning_directory,"*", "output", "gtdbtk_output.filtered.tsv"),
#     os.path.join(directories["preprocessing"], "prokaryotic_taxonomy.tsv"),
#     ),

        os.environ["concatenate_dataframes.py"],
        "--axis 0",
        "--allow_empty_or_missing_files",
        os.path.join(opts.prokaryotic_binning_directory,"*", "output", "gtdbtk_output.filtered.tsv"),
        ">",
        os.path.join(directories["output"], "prokaryotic_taxonomy.tsv"),
    ]
    return cmd



# GTDB-Tk
def get_gtdbtk_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
#gtdbtk classify_wf --genome_dir genomes --out_dir gtdbtk_output_batch1-3 -x fa --cpus 32"

    # Command
    cmd = [
        "for FP in $(cat %s); do cp $FP %s; done"%(
            opts.prokaryotic_binning_directory,
            directories["tmp"],
            ),
        "&&",
    ]
    if opts.gtdbtk_database:
        cmd += [ 
        "export GTDBTK_DATA_PATH={}".format(opts.gtdbtk_database),
        "&&",
        ]

    cmd += [
        # "mkdir -p {}".format(os.path.join(directories["tmp"], "gtdbtk")),
        # "&&",
        os.environ["gtdbtk"],
        "classify_wf",
        "--genome_dir {}".format(os.path.join(directories["tmp"])),
        "--out_dir {}".format(output_directory),
        "-x {}".format(opts.extension),
        "--cpus {}".format(opts.n_jobs),
        "--tmpdir {}".format(opts.tmpdir),
        opts.gtdbtk_options,
        "&&",

        os.environ["concatenate_dataframes.py"],
        "--axis 0",
        "--allow_empty_or_missing_files",
        os.path.join(output_directory, "classify", "gtdbtk.ar122.summary.tsv"),
        os.path.join(output_directory, "classify", "gtdbtk.bac120.summary.tsv"),
        ">",
        os.path.join(directories["output"], "prokaryotic_taxonomy.tsv"),
        "&&",
        "rm -rf {}".format(os.path.join(directories["tmp"], "*")),
    ]

    return cmd

def get_consensus_cluster_classification_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [ 

        # cat gtdbtk.summary.tsv | cut -f1,3,18 | tail -n +2 |
        "cat",
        input_filepaths[0],
        "|",
        "cut -f1,2,17",
        "|",
        "tail -n +2",
        "|",
        os.environ["insert_column_to_table.py"],
        "-c {}".format(opts.clusters),
        "-n id_slc",
        "-i 0",
        "|",
        os.environ["consensus_genome_classification.py"],
        "--leniency {}".format(opts.leniency),
        "-o {}".format(output_filepaths[0]),


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
        "concatenate_dataframes.py",
        "consensus_genome_classification.py",
        "insert_column_to_table.py",
    ])

    
    required_executables={
                # 1
                "gtdbtk",
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

    if os.path.isdir(opts.prokaryotic_binning_directory):
        # ==========
        # Preprocess
        # ==========
        step = 1

        program = "concatenate_gtdbtk"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories["output"]

        # Info
        description = "Merge GTDB-Tk results"


        # i/o
        input_filepaths = [
            opts.prokaryotic_binning_directory,
            ]
        output_filenames = ["prokaryotic_taxonomy.tsv"]
        output_filepaths = list(map(lambda filename: os.path.join(directories["output"], filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_concatenate_gtdbtk_cmd(**params)

        pipeline.add_step(
                    id=program_label,
                    description = description,
                    step=step,
                    cmd=cmd,
                    input_filepaths = input_filepaths,
                    output_filepaths = output_filepaths,
                    validate_inputs=True,
                    validate_outputs=True,
        )


    else:
  
        # ==========
        # GTDB-Tk
        # ==========
        step = 1

        program = "gtdbtk"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

        # Info
        description = "GTDB-Tk classification"


        # i/o
        input_filepaths = [
            opts.prokaryotic_binning_directory,
            ]
        output_filenames = ["prokaryotic_taxonomy.tsv"]
        output_filepaths = list(map(lambda filename: os.path.join(directories["output"], filename), output_filenames))

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
        )

    if opts.clusters:
        # ==========
        # consensus genome classification
        # ==========
        step = 2

        program = "consensus_cluster_classification"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories["output"]# = create_directory(os.path.join(directories["intermediate"], program_label))

        # Info
        description = "Consensus genome classification"


        # i/o
        input_filepaths = [
            os.path.join(directories["output"], "prokaryotic_taxonomy.tsv"),
            opts.clusters,
            ]
        output_filenames = ["prokaryotic_taxonomy.clusters.tsv"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_consensus_cluster_classification_cmd(**params)

        pipeline.add_step(
                    id=program_label,
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
    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <prokaryotic_binning_directory> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    # parser_io.add_argument("-i","--genomes", required=False, type=str, help = "path/to/genomes.list")
    parser_io.add_argument("-i","--prokaryotic_binning_directory", type=str, required=True, help = "path/to/prokaryotic_binning_directory")

    parser_io.add_argument("-c","--clusters", type=str, help = "path/to/clusters.tsv, Format: [id_mag]<tab>[id_cluster], No header.")
    parser_io.add_argument("-o","--output_directory", type=str, default="veba_output/classify/prokaryotic", help = "path/to/output_directory [Default: veba_output/classify/prokaryotic]")
    parser_io.add_argument("-x", "--extension", type=str, default="fa", help = "Fasta file extension for genomes if a list is provided [Default: fa]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("--tmpdir", type=str,  default=os.environ["TMPDIR"], help="path/to/TMPDIR [Default: {}]".format(os.environ["TMPDIR"]))  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Databases
    parser_databases = parser.add_argument_group('Database arguments')
    parser_databases.add_argument("--veba_database", type=str,  help=f"VEBA database location.  [Default: $VEBA_DATABASE environment variable]")

    # GTDB-Tk
    parser_gtdbtk = parser.add_argument_group('GTDB-Tk arguments')
    # parser_gtdbtk.add_argument("--gtdbtk_database", type=str, default=DATABASE_GTDBTK, help="GTDB-Tk | path/to/gtdbtk_database (e.g. --arg 1 ) [Default: {}]".format(DATABASE_GTDBTK))
    parser_gtdbtk.add_argument("--gtdbtk_options", type=str, default="", help="GTDB-Tk | classify_wf options (e.g. --arg 1 ) [Default: '']")


    # Consensus genome classification
    parser_consensus = parser.add_argument_group('Consensus genome classification arguments')

    parser_consensus.add_argument("-l","--leniency", default=1.382, type=float, help = "Leniency parameter. Lower value means more conservative weighting. A value of 1 indiciates no weight bias. A value greater than 1 puts higher weight on higher level taxonomic assignments. A value less than 1 puts lower weights on higher level taxonomic assignments.  [Default: 1.382]")
    # parser_consensus.add_argument("-r", "--rank_prefixes", type=str, default=RANK_PREFIXES, help = "Rank prefixes separated by , delimiter'\n[Default: {}]".format(RANK_PREFIXES))
    # parser_consensus.add_argument("-d", "--delimiter", type=str, default=";", help = "Taxonomic delimiter [Default: ; ]")
    # parser_consensus.add_argument("-s", "--simple", action="store_true", help = "Simple classification that does not use lineage information from --rank_prefixes")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Database
    if opts.veba_database is None:
        assert "VEBA_DATABASE" in os.environ, "Please set the following environment variable 'export VEBA_DATABASE=/path/to/veba_database' or provide path to --veba_database"
        opts.veba_database = os.environ["VEBA_DATABASE"]

    opts.gtdbtk_database = os.path.join(opts.veba_database, "Classify", "GTDBTk")


    # Directories
    directories = dict()
    directories["project"] = create_directory(opts.output_directory)
    directories["output"] = create_directory(os.path.join(directories["project"], "output"))
    directories["log"] = create_directory(os.path.join(directories["project"], "log"))
    directories["tmp"] = create_directory(os.path.join(directories["project"], "tmp"))
    directories["checkpoints"] = create_directory(os.path.join(directories["project"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["project"], "intermediate"))
    # os.environ["TMPDIR"] = directories["tmp"]

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
