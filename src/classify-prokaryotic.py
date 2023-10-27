#!/usr/bin/env python
from __future__ import print_function, division
from pdb import Restart
import sys, os, argparse, glob
from collections import OrderedDict, defaultdict
# from Bio.SeqIO.FastaIO import SimpleFastaParser

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from genopype import __version__ as genopype_version
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.10.16"

# GTDB-Tk
def get_gtdbtk_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
#gtdbtk classify_wf --genome_dir genomes --out_dir gtdbtk_output_batch1-3 -x fa --cpus 32"

    # Command
    cmd = [ 
        "export GTDBTK_DATA_PATH={}".format(opts.gtdbtk_database),

        "&&",

        "mkdir -p {}".format(os.path.join(directories["tmp"], "gtdbtk")),

        "&&",

        "mkdir -p {}".format(os.path.join(directories["tmp"], "genomes")),
    ]

    if opts.prokaryotic_binning_directory:
        cmd += [ 
                "&&",

            "ls {}".format(os.path.join(opts.prokaryotic_binning_directory, "*", "output", "genomes", "*.fa")),
            ">",
            os.path.join(directories["tmp"], "genomes.list"),
        ]

    if opts.genomes:
        cmd += [ 
                "&&",

            "cp",
            opts.genomes,
            os.path.join(directories["tmp"], "genomes.list"),
        ]

    # GTDB-Tk
    cmd += [

            "&&",

        "(for FP in $(cat %s); do cp $FP %s; done)"%(
            os.path.join(directories["tmp"], "genomes.list"),
            os.path.join(directories["tmp"],"genomes"),
            ),

            "&&",

        os.environ["gtdbtk"],
        "classify_wf",
        "--genome_dir {}".format(os.path.join(directories["tmp"], "genomes")),
        "--out_dir {}".format(output_directory),
        "-x {}".format(opts.extension),
        "--cpus {}".format(opts.n_jobs),
        "--tmpdir {}".format(os.path.join(directories["tmp"], "gtdbtk")),
        opts.gtdbtk_options,
    ]

    if opts.skip_ani_screen:
        cmd += [
        "--skip_ani_screen",
        ]

    else:
        cmd += [
        "--mash_db {}".format(os.path.join(opts.gtdbtk_database, "mash","gtdb_r214.msh")),
        ]
        

    # Concatenation
    cmd += [

            "&&",

        os.environ["concatenate_dataframes.py"],
        "--axis 0",
        "--allow_empty_or_missing_files",
        os.path.join(output_directory, "classify", "gtdbtk.ar122.summary.tsv"),
        os.path.join(output_directory, "classify", "gtdbtk.bac120.summary.tsv"),
        ">",
        os.path.join(directories["output"], "taxonomy.tsv"),

            "&&",

        "rm -rf {}".format(os.path.join(directories["tmp"], "*")),
    ]

    return cmd

def get_krona_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [ 
        os.environ["compile_krona.py"],
        "-i {}".format(input_filepaths[0]),
        "-m prokaryotic",
        "-o {}".format(os.path.join(output_directory, "krona.tsv")),

            "&&",

        os.environ["ktImportText"],
        "-n='Prokaryotic Taxonomy'",
        "-o {}".format(os.path.join(output_directory, "krona.html")),
        os.path.join(output_directory, "krona.tsv"),

            "&&",

        "SRC={}; DST={}; SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST".format(
        os.path.join(output_directory, "krona.html"),
        directories["output"],
        )
    ]
    return cmd

def get_consensus_cluster_classification_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [ 
        os.environ["compile_prokaryotic_genome_cluster_classification_scores_table.py"],
        "-i {}".format(input_filepaths[0]),
        "-c {}".format(input_filepaths[1]),
        "|",
        os.environ["consensus_genome_classification.py"],
        "--leniency {}".format(opts.leniency),
        "-o {}".format(output_filepaths[0]),
        "-u 'Unclassified prokaryote'",
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
        "compile_prokaryotic_genome_cluster_classification_scores_table.py",
        # "cut_table_by_column_labels.py",
        "concatenate_dataframes.py",
        "consensus_genome_classification.py",
        # "insert_column_to_table.py",
        "compile_krona.py",

    ])

    
    required_executables={
                # 1
                "gtdbtk",

                # Krona
                "ktImportText",

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
    if opts.prokaryotic_binning_directory:
        input_filepaths = [
            opts.prokaryotic_binning_directory,
            ]
    if opts.genomes:
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

    cmd = get_gtdbtk_cmd(**params)

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
    # Krona
    # ==========
    step = 2

    program = "krona"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Krona graph for prokaryotic classification"


    # i/o
    input_filepaths = output_filepaths
    output_filenames = ["krona.tsv", "krona.html"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_krona_cmd(**params)

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
        step = 3

        program = "consensus_cluster_classification"
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

        cmd = get_consensus_cluster_classification_cmd(**params)

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
    assert bool(opts.prokaryotic_binning_directory) != bool(opts.genomes), "Must choose either --prokaryotic_binning_directory or --genomes, not both."

    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <prokaryotic_binning_directory>|-g <genomes.list> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i","--prokaryotic_binning_directory", type=str, required=False, help = "path/to/prokaryotic_binning_directory [Cannot be used with --genomes]")
    parser_io.add_argument("-g","--genomes", required=False, type=str, help = "path/to/genomes.list [Cannot be ued with --prokaryotic_binning_directory]")
    parser_io.add_argument("-c","--clusters", type=str, help = "path/to/clusters.tsv, Format: [id_mag]<tab>[id_cluster], No header.")
    parser_io.add_argument("-o","--output_directory", type=str, default="veba_output/classify/prokaryotic", help = "path/to/output_directory [Default: veba_output/classify/prokaryotic]")
    parser_io.add_argument("-x", "--extension", type=str, default="fa", help = "Fasta file extension for genomes if a list is provided [Default: fa]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("--tmpdir", type=str, help="path/to/TMPDIR")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Databases
    parser_databases = parser.add_argument_group('Database arguments')
    parser_databases.add_argument("--veba_database", type=str,  help=f"VEBA database location.  [Default: $VEBA_DATABASE environment variable]")

    # GTDB-Tk
    parser_gtdbtk = parser.add_argument_group('GTDB-Tk arguments')
    parser_gtdbtk.add_argument("--skip_ani_screen", action="store_true", help = "Skip ANI screen [Default: Don't skip ANI screen]")
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
    directories["checkpoints"] = create_directory(os.path.join(directories["project"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["project"], "intermediate"))
    if not opts.tmpdir:
        opts.tmpdir = create_directory(os.path.join(directories["project"], "tmp"))
    directories["tmp"] = opts.tmpdir
    os.environ["TMPDIR"] = directories["tmp"]

    # Database
    if opts.veba_database is None:
        assert "VEBA_DATABASE" in os.environ, "Please set the following environment variable 'export VEBA_DATABASE=/path/to/veba_database' or provide path to --veba_database"
        opts.veba_database = os.environ["VEBA_DATABASE"]

    opts.gtdbtk_database = os.path.join(opts.veba_database, "Classify", "GTDB")

    # if opts.mash_database is None:
    #     opts.mash_database = os.path.join(directories["tmp"], "gtdbtk.msh")
    

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
