#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
import pandas as pd
from genopype import * 
from soothsayer_utils import *

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.2.1"

# ==============
# Agostic commands
# ==============

# concatenate fasta
def get_concatenate_fasta_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]
    # Command 
    cmd = [ 
        os.environ["concatenate_fasta.py"],
        "-i {}".format(input_filepaths[0]),
        "-o {}".format(output_directory),
        "-m {}".format(opts.minimum_contig_length),
        "-x {}".format("fa.gz"),
        "-b reference",
        "-M {}".format(opts.mode),


    ]
    return cmd

# concatenate gene models
def get_concatenate_gff_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]
    # Command 
    cmd = [ 
        os.environ["concatenate_gff.py"],
        "-i {}".format(input_filepaths[0]),
        "-o {}".format(output_directory),
        "-x {}".format("gff"),
        "-b reference",
        "-M {}".format(opts.mode),

    ]
    return cmd

# ==============
# Local commands
# ==============
def get_bowtie2_local_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]
    # Command 
    cmd = [
"""

for ID_SAMPLE in $(cut -f1 %s); 
    do %s --threads %d --seed %d %s/${ID_SAMPLE}/reference.fa.gz %s/${ID_SAMPLE}/reference.fa.gz
    done
"""%(
    opts.references,
    os.environ["bowtie2-build"],
    opts.n_jobs,
    opts.random_state,
    output_directory,
    output_directory,
    ),
    ]

    return cmd

# ==============
# Global commands
# ==============
def get_bowtie2_global_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]
    # Command 
    cmd = [ 
        os.environ["bowtie2-build"],
        "--threads {}".format(opts.n_jobs),
        "--seed {}".format(opts.random_state),
        opts.bowtie2_build_options,
        input_filepaths[0],
        input_filepaths[0],
    ]
    return cmd

# Pipeline
def create_local_pipeline(opts, directories, f_cmds):

    # .................................................................
    # Primordial
    # .................................................................
    # Commands file
    pipeline = ExecutablePipeline(name=__program__,  f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])
    

    # ==========
    # concatenate fasta
    # ==========
    
    step = 1

    program = "concatenate_fasta"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories["output"]# = create_directory(os.path.join(directories["intermediate"], "concatenated"))


    # Info
    description = "Concatenate fasta files"
    # i/o
    input_filepaths = [
        opts.references,
    ]

    output_filenames = [
        "*/reference.fa.gz",
        "*/reference.saf",

    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    } 

    cmd = get_concatenate_fasta_cmd(**params)
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
                log_prefix=program_label,
                # acceptable_returncodes= {0,1},

    )

    # =======================
    # concatenate gene models
    # =======================
    
    step = 2

    program = "concatenate_gff"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories["output"]# = create_directory(os.path.join(directories["intermediate"], "concatenated"))


    # Info
    description = "Concatenate gene models"
    # i/o
    input_filepaths = [
        opts.gene_models,
    ]

    output_filenames = [
        "*/reference.gff",
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    } 

    cmd = get_concatenate_gff_cmd(**params)
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
                log_prefix=program_label,
                # acceptable_returncodes= {0,1},

    )

    # ==========
    # Bowtie2 index
    # ==========
    step = 3

    program = "bowtie2"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories["output"] #= create_directory(os.path.join(directories["intermediate"], "concatenated"))


    # Info
    description = "Build mapping index"
    # i/o
    input_filepaths = list(
        map(lambda id_sample: os.path.join(directories["output"], id_sample, "reference.fa.gz"), 
        opts.samples,
        ),
    )
    

    output_filepaths = list(
        map(lambda fp: "{}.*.bt2".format(fp), input_filepaths),
    )

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    } 

    cmd = get_bowtie2_local_cmd(**params)
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
                log_prefix=program_label,
                # acceptable_returncodes= {0,1},

    )

 
    return pipeline

def create_global_pipeline(opts, directories, f_cmds):
    # .................................................................
    # Primordial
    # .................................................................
    # Commands file
    pipeline = ExecutablePipeline(name=__program__,  f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])
    

    # ==========
    # concatenate fasta
    # ==========
    
    step = 1

    program = "concatenate_fasta"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories["output"] #= create_directory(os.path.join(directories["intermediate"], "concatenated"))


    # Info
    description = "Concatenate fasta files"
    # i/o
    input_filepaths = [
        opts.references,
    ]

    output_filenames = [
        "reference.fa.gz",
        "reference.saf",
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    } 

    cmd = get_concatenate_fasta_cmd(**params)
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
                log_prefix=program_label,
                # acceptable_returncodes= {0,1},

    )

    # ==========
    # concatenate gff
    # ==========
    
    step = 2

    program = "concatenate_gff"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories["output"] #= create_directory(os.path.join(directories["intermediate"], "concatenated"))


    # Info
    description = "Concatenate gff files"
    # i/o
    input_filepaths = [
        opts.gene_models,
    ]

    output_filenames = [
        "reference.gff",
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    } 

    cmd = get_concatenate_gff_cmd(**params)
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
                log_prefix=program_label,
                # acceptable_returncodes= {0,1},

    )

    # ==========
    # Bowtie2 index
    # ==========
    
    step = 3

    program = "bowtie2"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories["output"] #= create_directory(os.path.join(directories["intermediate"], "concatenated"))


    # Info
    description = "Build mapping index"
    # i/o
    input_filepaths = [
        os.path.join(directories["output"], "reference.fa.gz"),
    ]

    output_filenames = [
        "reference.fa.gz.*.bt2",
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    } 

    cmd = get_bowtie2_global_cmd(**params)
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
                log_prefix=program_label,
                # acceptable_returncodes= {0,1},

    )


    return pipeline

# ============
# Run Pipeline
# ============
# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = set([ 
        "concatenate_fasta.py",
        "concatenate_gff.py",
        # "fasta_to_saf.py",
    ])

    
    required_executables = set([ 
         "bowtie2-build",
    ])| accessory_scripts

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

# Configure parameters
def configure_parameters(opts, directories):
    df_references = pd.read_csv(opts.references, sep="\t", header=None)
    df_gene_models = pd.read_csv(opts.gene_models, sep="\t", header=None)

    assert df_references.shape[1] == df_gene_models.shape[1], "--references and --gene_models must have the same number of columns"

    m = df_references.shape[1]
    if opts.mode == "infer":
        assert_acceptable_arguments(m, {1,2})
        opts.mode = {1:"global", 2:"local"}[m]
        print("Inferring mode is {}".format(opts.mode), file=sys.stderr)

    if opts.mode == "global":
        assert m == 1, "There should be only one column if mode='global'"
    if opts.mode == "local":
        assert m == 2, "There should be two columns if mode='local'"
    if opts.mode == "local":
        opts.samples = sorted(set(df_references.iloc[:,0]))



    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__

    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <mags> -o <output> --heatmap_output <pdf> ".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-r","--references", type=str, required=True, help = "local mode: [id_sample]<tab>[path/to/reference.fa] and global mode: [path/to/reference.fa]")
    parser_io.add_argument("-g","--gene_models", type=str, required=True, help = "local mode: [id_sample]<tab>[path/to/reference.gff] and global mode: [path/to/reference.gff]")
    parser_io.add_argument("-o","--output_directory", type=str, default="veba_output/index", help = "path/to/project_directory [Default: veba_output/index]")
    parser_io.add_argument("-m", "--minimum_contig_length", type=int, default=1500, help="Minimum contig length [Default: 1500]")
    parser_io.add_argument("-M", "--mode", type=str, default="infer", help="Concatenate all references with global and build index or build index for each reference {global, local, infer}")
    # parser_io.add_argument("-c", "--copy_files", action="store_true", help="Copy files instead of symlinking. Only applies to global.")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Use -1 for completely random. Use 0 for consecutive random states.  Use any other positive integer for the same random state for all iterations [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    # parser_utility.add_argument("--skip_concatenation", action="store_true", help="Skip concatenation step. Useful when references are concatenated before hand (e.g., already ran for bowtie2 but want STAR index as well)")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Utility
    parser_bowtie2 = parser.add_argument_group('Bowtie2 Index arguments')
    parser_bowtie2.add_argument("--bowtie2_build_options", type=str, default="", help="bowtie2-build | More options (e.g. --arg 1 ) [Default: '']")

    # parser_star = parser.add_argument_group('STAR arguments')
    # parser_star.add_argument("--read_length", type=int, default=151, help = "Read length [Default: 151]")
    # parser_star.add_argument("--star_index_options", type=str, default="", help="bowtie2-build | More options (e.g. --arg 1 ) [Default: '']")

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
        if opts.mode == "local":
            pipeline = create_local_pipeline(
                        opts=opts,
                        directories=directories,
                        f_cmds=f_cmds,
            )
        if opts.mode == "global":
            pipeline = create_global_pipeline(
                        opts=opts,
                        directories=directories,
                        f_cmds=f_cmds,
            )
        pipeline.compile()
        pipeline.execute(restart_from_checkpoint=opts.restart_from_checkpoint)



if __name__ == "__main__":
    main()
