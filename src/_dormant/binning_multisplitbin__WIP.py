#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
from collections import OrderedDict, defaultdict

import pandas as pd
import numpy as np

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.07.30"

# .............................................................................
# Primordial
# .............................................................................
# Bowtie2
def get_vamb_multisplit_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]
    # Command
    cmd = [
   
    # VAMB
    "(",
    "rm -r {}".format(output_directory), # There can't be an existing directory for some reason
    "&&",
    os.environ["vamb"],
    "--fasta {}".format(input_filepaths[0]),
    "--jgi {}".format(input_filepaths[1]),
    "--minfasta {}".format(opts.minimum_genome_size),
    "-m {}".format(opts.minimum_contig_size),
    "-o __".format(opts.id_separator),
    "--outdir {}".format(output_directory),
    opts.vamb_options,
    ")",
    ]
    return cmd



# # Symlink
# def get_symlink_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
#     # Command
#     cmd = ["("]
#     for filepath in input_filepaths:
#         cmd.append("ln -f -s {} {}".format(os.path.realpath(filepath), os.path.realpath(output_directory)))
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

    required_executables={
                "vamb",
     }

    if opts.path_config == "CONDA_PREFIX":
        executables = dict()
        for name in required_executables:
            executables[name] = os.path.join(os.environ["CONDA_PREFIX"], "bin", name)
    else:
        if opts.path_config is None:
            opts.path_config = os.path.join(opts.script_directory, "metagenomics_pipeline_config.tsv")
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
    accessory_scripts = []
    for name in accessory_scripts:
        executables[name] = "python " + os.path.join(opts.script_directory, name)
    print(format_header( "Adding executables to path from the following source: {}".format(opts.path_config), "-"), file=sys.stdout)
    for name, executable in executables.items():
        if name in required_executables:
            print(name, executable, sep = " --> ", file=sys.stdout)
            os.environ[name] = executable.strip()
    print("", file=sys.stdout)



def create_pipeline(opts, directories, f_cmds):

    # .................................................................
    # Primordial
    # .................................................................
    # Commands file
    pipeline = ExecutablePipeline(name="binning.py", description=opts.name, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])

    # ==========
    # VAMB
    # ==========
    program = "vamb_multi-split"
    # Add to directories
    output_directory = directories[("intermediate",  program)] = create_directory(os.path.join(directories["intermediate"], "{}_output".format(program)))

    # Info
    step = 1
    description = "Multibinning"

    # i/o
    input_filepaths = [opts.fasta, opts.coverage]

    output_filenames = ["clusters.tsv", "bins"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_vamb_multisplit_cmd(**params)
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

    

    # # =============
    # # Symlink
    # # =============
    # program = "symlink"
    # # Add to directories
    # output_directory = directories["output"]

    # # Info
    # step = 3
    # description = "Symlinking relevant output files"

    # # i/o
    # input_filepaths = [
    #     os.path.join(directories[("intermediate", "bowtie2")], "mapped.sorted.bam"),
    #     os.path.join(directories[("intermediate", "coverage")], "coverage.tsv"),

    # ]

    # output_filenames =  map(lambda fp: fp.split("/")[-1], input_filepaths)
    # output_filepaths = list(map(lambda fn:os.path.join(directories["output"], fn), output_filenames))
    #     # Prodigal
    #     # os.path.join(directories["output"], "*"),
    
    # params = {
    # "input_filepaths":input_filepaths,
    # "output_filepaths":output_filepaths,
    # "output_directory":output_directory,
    # "opts":opts,
    # "directories":directories,
    # }

    # cmd = get_symlink_cmd(**params)
    # pipeline.add_step(
    #         id=program,
    #         description = description,
    #         step=step,
    #         cmd=cmd,
    #         input_filepaths = input_filepaths,
    #         output_filepaths = output_filepaths,
    #         validate_inputs=True,
    #         validate_outputs=False,
    # )

    return pipeline

# Configure parameters
def configure_parameters(opts, directories):
    # assert opts.reference_assembly is not None, "Must include --reference_assembly"

    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -f <scaffolds.fasta> -c <coverage.tsv -o <output_directory>".format(__program__)

    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-o","--project_directory", type=str, default="./binning_output", help = "path/to/project_directory [Default: ./binning_output]")

    # VAMB
    parser_vamb = parser.add_argument_group('VAMB arguments')
    parser_vamb.add_argument("-f","--fasta", type=str, required=True, help = "path/to/scaffolds_catalogue.fasta")
    parser_vamb.add_argument("-c","--coverage", type=str, required=True, help = "path/to/jgi_coverage.tsv")
    parser_vamb.add_argument("-m", "--minimum_contig_size", type=int, default=1500, help = "Minimum contig size [Default: 1500" )
    parser_vamb.add_argument("-g", "--minimum_genome_size", type=int, default=200000, help = "Minimum genome size of bins [Default: 200000")
    parser_vamb.add_argument("-s", "--id_separator", type=str, default="__", help = "Identifier separator for scaffolds.  e.g., <sample>__<contig> [Default: __" )
    parser_vamb.add_argument("--vamb_options", type=str, default="", help="vamb | More options (e.g. --arg 1 ) [Default: '']")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    # parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    # parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Options
    opts = parser.parse_args()

    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    opts.name = "catalogue"

    # Directories
    directories = dict()
    directories["project"] = create_directory(opts.project_directory)
    directories["sample"] = create_directory(os.path.join(directories["project"], opts.name))
    # directories["preprocessing"] = create_directory(os.path.join(directories["sample"], "preprocessing"))
    directories["output"] = create_directory(os.path.join(directories["sample"], "output"))
    directories["log"] = create_directory(os.path.join(directories["sample"], "log"))
    directories["tmp"] = create_directory(os.path.join(directories["sample"], "tmp"))
    directories["checkpoints"] = create_directory(os.path.join(directories["sample"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["sample"], "intermediate"))


    # Info
    print(format_header("binning.py", "="), file=sys.stdout)
    print(format_header("Configuration:", "-"), file=sys.stdout)
    print(format_header("Name: {}".format(opts.name), "."), file=sys.stdout)
    print("Python version:", sys.version.replace("\n"," "), file=sys.stdout)
    print("Python path:", sys.executable, file=sys.stdout) #sys.path[2]
    print("Script version:", __version__, file=sys.stdout)
    print("Moment:", get_timestamp(), file=sys.stdout)
    print("Directory:", os.getcwd(), file=sys.stdout)
    print("Commands:", list(filter(bool,sys.argv)),  sep="\n", file=sys.stdout)
    configure_parameters(opts, directories)
    sys.stdout.flush()

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
