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
__version__ = "2023.5.8"


# MetaEuk
def get_metaeuk_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [

        # Run MetaEuk
        os.environ["metaeuk"],
        "easy-predict",
        "--threads {}".format(opts.n_jobs),
        "-s {}".format(opts.metaeuk_sensitivity),
        "-e {}".format(opts.metaeuk_evalue),
        opts.metaeuk_options,
        opts.fasta, # contigs
        opts.metaeuk_database, # db
        os.path.join(output_directory, "metaeuk"), # output prefix
        os.path.join(directories["tmp"],"metaeuk"),

        # Convert MetaEuk identifiers
        "&&",
        os.environ["compile_metaeuk_identifiers.py"],
        "--cds {}".format(os.path.join(output_directory, "metaeuk.codon.fas")),
        "--protein {}".format(os.path.join(output_directory, "metaeuk.fas")),
        "-o {}".format(output_directory),
        "-b {}".format(opts.basename),
    ]

    # Remove temporary files
    cmd += [
        "&&",
        "rm -rf {} {} {}".format(
            os.path.join(output_directory, "*.fas"), # output prefix
            os.path.join(output_directory, "metaeuk.gff"), # output prefix
            os.path.join(directories["tmp"],"metaeuk", "*"),
        ),
    ]

    if opts.scaffolds_to_bins:
        cmd += [
        # Partition the gene models and genomes
        "&&",
        os.environ["partition_gene_models.py"],
        "-i {}".format(opts.scaffolds_to_bins),
        "-f {}".format(opts.fasta),
        "-g {}".format(os.path.join(output_directory, "{}.gff".format(opts.basename))),
        "-d {}".format(os.path.join(output_directory, "{}.ffn".format(opts.basename))),
        "-a {}".format(os.path.join(output_directory, "{}.faa".format(opts.basename))),
        "-o {}".format(os.path.join(output_directory, "genomes")),
        "--use_mag_as_description",

        "&&",
        "rm -rf {}".format(os.path.join(output_directory, "{}.*".format(opts.basename))),
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
    # Binning
    # ==========
    step = 1

    program = "metaeuk"

    program_label = "{}__{}".format(step, program)

    # Add to directories


    # Info
    description = "Modeling genes with via MetaEuk"
    
    # i/o
    input_filepaths = [
        opts.fasta, 
        ]

    output_directory = directories["output"]
    output_filepaths = [ 
        os.path.join(output_directory, "identifier_mapping.metaeuk.tsv"),
        os.path.join(output_directory, "metaeuk_to_simple.tsv"),
    ]
    if not opts.scaffolds_to_bins:
        output_filepaths += [
            os.path.join(output_directory, "{}.gff".format(opts.basename)),
            os.path.join(output_directory, "{}.ffn".format(opts.basename)),
            os.path.join(output_directory, "{}.faa".format(opts.basename)),
        ]
    else:
        output_filepaths += [
            os.path.join(output_directory, "genomes","*.gff"),
            os.path.join(output_directory, "genomes","*.ffn"),
            os.path.join(output_directory, "genomes","*.faa"),
        ]

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_metaeuk_cmd(**params)


    # Add step to pipeline
    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
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
        "partition_gene_models.py",
        "compile_metaeuk_identifiers.py",
    }

    required_executables={
        "metaeuk",
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
        executables[name] = "'{}'".format(os.path.join(opts.script_directory, name)) # Can handle spaces in path



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
    usage = "{} -f <scaffolds.fasta> -d <metaeuk_database> -i <scaffolds_to_bins.tsv>  -o <output_directory>".format(__program__)

    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-f","--fasta", type=str, required=True, help = "path/to/scaffolds.fasta")
    parser_io.add_argument("-d", "--metaeuk_database", type=str,  required=True, help=f"MetaEuk/MMSEQS2 database")
    parser_io.add_argument("-o","--output_directory", type=str, default="metaeuk_output", help = "path/to/project_directory [Default: metaeuk_output]")
    parser_io.add_argument("-i","--scaffolds_to_bins", type=str, required=False,  help = "path/to/scaffolds_to_bins.tsv, [Optional] Format: [id_scaffold]<tab>[id_bin], No header")
    parser_io.add_argument("-b","--basename", type=str, default="gene_models", required=False,  help = "Basename for output files. [Default: gene_models]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    # parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # MetaEuk
    parser_metaeuk = parser.add_argument_group('MetaEuk arguments')
    parser_metaeuk.add_argument("--metaeuk_sensitivity", type=float, default=4.0, help="MetaEuk | Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive  [Default: 4.0]")
    parser_metaeuk.add_argument("--metaeuk_evalue", type=float, default=0.01, help="MetaEuk | List matches below this E-value (range 0.0-inf) [Default: 0.01]")
    parser_metaeuk.add_argument("--metaeuk_options", type=str, default="", help="MetaEuk | More options (e.g. --arg 1 ) [Default: ''] https://github.com/soedinglab/metaeuk")

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
    directories["output"] = create_directory(opts.output_directory)
    directories["intermediate"] = create_directory(os.path.join(directories["output"], "intermediate"))
    directories["log"] = create_directory(os.path.join(directories["output"], "log"))
    directories["tmp"] = create_directory(os.path.join(directories["output"],"intermediate", "tmp"))
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

    shutil.rmtree(directories["intermediate"])
    
   


if __name__ == "__main__":
    main(sys.argv[1:])
