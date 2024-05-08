#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, shutil, time, warnings
from multiprocessing import cpu_count
from collections import OrderedDict, defaultdict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from genopype import __version__ as genopype_version
from soothsayer_utils import *

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.12.15"

# ============
# Run Pipeline
# ============
# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = set([ 

    ])

    required_executables={
                "sylph",

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
    usage = "{} -i <input.tsv> -o <output_directory>".format(__program__)

    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i","--input", type=str, default="stdin", help = "path/to/input.tsv. Format: Must include the following columns (No header)[organism_type]<tab>[path/to/genome.fa]. You can get this from `cut -f1,4 veba_output/misc/genomes_table.tsv` [Default: stdin]")
    parser_io.add_argument("-o","--output_directory", type=str, default="veba_output/profiling/databases", help = "path/to/output_directory for databases [Default: veba_output/profiling/databases]")
    parser_io.add_argument("--viral_tag", type=str, default="viral", help = "[Not case sensitive] Tag/Label of viral organisms in first column of --input (e.g., viral, virus, viron) [Default: viral]")


    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
    # parser_utility.add_argument("--verbose", action='store_true')

    # Sylph
    parser_sylph = parser.add_argument_group('Sylph sketch arguments')
    parser_sylph.add_argument("-k", "--sylph_k", type=int, choices={21,31}, default=31,  help="Sylph |  Value of k. Only k = 21, 31 are currently supported. [Default: 31]")
    parser_sylph.add_argument("-s", "--sylph_minimum_spacing", type=int,  default=30,  help="Sylph |  Minimum spacing between selected k-mers on the genomes [Default: 30]")

    parser_sylph_nonviral = parser.add_argument_group('[Prokaryotic & Eukaryotic] Sylph sketch arguments')
    parser_sylph_nonviral.add_argument("--sylph_nonviral_subsampling_rate", type=int, default=200,  help="Sylph [Prokaryotic & Eukaryotic]|  Subsampling rate.	[Default: 200]")
    parser_sylph_nonviral.add_argument("--sylph_nonviral_options", type=str, default="", help="Sylph [Prokaryotic & Eukaryotic] | More options for `sylph sketch` (e.g. --arg 1 ) [Default: '']")

    parser_sylph_viral = parser.add_argument_group('[Viral] Sylph sketch arguments')
    parser_sylph_viral.add_argument("--sylph_viral_subsampling_rate", type=int, default=100,  help="Sylph [Viral]|  Subsampling rate. [Default: 100]")
    parser_sylph_viral.add_argument("--sylph_viral_options", type=str, default="", help="Sylph [Viral] | More options for `sylph sketch` (e.g. --arg 1 ) [Default: '']")

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
    directories["output"] = create_directory(opts.output_directory)
    directories["intermediate"] = create_directory(os.path.join(directories["output"], "intermediate"))
    directories["log"] = create_directory(os.path.join(directories["intermediate"], "log"))
    directories["checkpoints"] = create_directory(os.path.join(directories["intermediate"], "checkpoints"))

    # Info
    print(format_header(__program__, "="), file=sys.stdout)
    print(format_header("Configuration:", "-"), file=sys.stdout)
    print("Python version:", sys.version.replace("\n"," "), file=sys.stdout)
    print("Python path:", sys.executable, file=sys.stdout) #sys.path[2]
    print("Script version:", __version__, file=sys.stdout)
    print("GenoPype version:", genopype_version, file=sys.stdout) #sys.path[2]
    print("Moment:", get_timestamp(), file=sys.stdout)
    print("Directory:", os.getcwd(), file=sys.stdout)
    print("Commands:", list(filter(bool,sys.argv)),  sep="\n", file=sys.stdout)
    configure_parameters(opts, directories)
    sys.stdout.flush()

    # Make directories
    t0 = time.time()
    # print(format_header("* ({}) Creating directories:".format(format_duration(t0)), opts.output_directory), file=sys.stdout)
    # os.makedirs(opts.output_directory, exist_ok=True)
    
    # Load input
    if opts.input == "stdin":
        opts.input = sys.stdin
    df_genomes = pd.read_csv(opts.input, sep="\t", header=None)
    assert df_genomes.shape[1] == 2, "Must include the follow columns (No header) [organism_type]<tab>[genome]).  Suggested input is from `compile_genomes_table.py` script using `cut -f1,4` to get the necessary columns."
    df_genomes.columns = ["organism_type", "genome"]

    opts.viral_tag = opts.viral_tag.lower()

    print(format_header("* ({}) Organizing genomes by organism_type".format(format_duration(t0))), file=sys.stdout)
    organism_to_genomes = defaultdict(set)
    for i, (organism_type, genome_filepath) in pv(df_genomes.iterrows(), unit="genomes ", total=df_genomes.shape[0]):
        organism_type = organism_type.lower()
        if organism_type == opts.viral_tag:
            organism_to_genomes["viral"].add(genome_filepath)
        else:
            organism_to_genomes["nonviral"].add(genome_filepath)
    # del df_genomes

    # Commands
    f_cmds = open(os.path.join(directories["intermediate"], "commands.sh"), "w")

    for organism_type, filepaths in organism_to_genomes.items():
        # Write genomes to file
        print(format_header("* ({}) Creating genome database: (N={}) for organism_type='{}'".format(format_duration(t0),len(filepaths), organism_type)), file=sys.stdout)

        genome_filepaths_list = os.path.join(directories["intermediate"], "{}_genomes.list".format(organism_type))
        with open(genome_filepaths_list, "w") as f:
            for fp in sorted(filepaths):
                print(fp, file=f)

        name = "sylph__{}".format(organism_type)
        description = "[Program = sylph sketch] [Organism_Type = {}]".format(organism_type)

        arguments = [
            os.environ["sylph"],
            "sketch",
            "-t {}".format(opts.n_jobs),
            "--gl {}".format(genome_filepaths_list),
            "-o {}".format(os.path.join(opts.output_directory, "genome_database-{}".format(organism_type))),
            "-k {}".format(opts.sylph_k),
            "--min-spacing {}".format(opts.sylph_minimum_spacing),
        ]

        if organism_type == "nonviral":
            arguments += [
            "-c {}".format(opts.sylph_nonviral_subsampling_rate),
            opts.sylph_nonviral_options,
        ]

        else:
            arguments += [
            "-c {}".format(opts.sylph_viral_subsampling_rate),
            opts.sylph_viral_options,
        ]
        print(arguments, file=sys.stdout)                
        cmd = Command(
            arguments,
            name=name, 
            f_cmds=f_cmds,
            )
        
    
        # Run command
        cmd.run(
            checkpoint_message_notexists="[Running ({})] | {}".format(format_duration(t0), description),
            checkpoint_message_exists="[Loading Checkpoint ({})] | {}".format(format_duration(t0), description),
            write_stdout=os.path.join(directories["log"], "{}.o".format(name)),
            write_stderr=os.path.join(directories["log"], "{}.e".format(name)),
            write_returncode=os.path.join(directories["log"], "{}.returncode".format(name)),
            checkpoint=os.path.join(directories["checkpoints"], name),
            )
        
        if hasattr(cmd, "returncode_"):
            if cmd.returncode_ != 0:
                print("[Error] | {}".format(description), file=sys.stdout)
                print("Check the following files:\ncat {}".format(os.path.join(directories["log"], "{}.*".format(name))), file=sys.stdout)
                sys.exit(cmd.returncode_)
            else:
                output_filepath = os.path.join(opts.output_directory, "genome_database-{}.syldb".format(organism_type))
                size_bytes = os.path.getsize(output_filepath)
                size_mb = size_bytes >> 20
                if size_mb < 1:
                    print("Output Database:", output_filepath, "({} bytes)".format(size_bytes), file=sys.stdout)
                else:
                    print("Output Database:", output_filepath, "({} MB)".format(size_mb), file=sys.stdout)

    f_cmds.close()

if __name__ == "__main__":
    main(sys.argv[1:])


