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
def get_bowtie2_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]

    # Command
    cmd = [
    # Clear temporary directory just in case
    "rm -rf {}".format(os.path.join(directories["tmp"], "*")),
    "&&",
    # Bowtie2
    "(",
    os.environ["bowtie2"],
   "-x {}".format(opts.bowtie2_index),
    "-1 {}".format(input_filepaths[0]),
    "-2 {}".format(input_filepaths[1]),
   "--threads {}".format(opts.n_jobs),
    # "--un-conc-gz {}".format(os.path.join(output_directory, "unmapped_%.{}.gz".format(unmapped_ext))),
    # "--un-gz {}".format(os.path.join(output_directory, "unmapped_singletons_%{}.gz".format(unmapped_ext))),
    "--seed {}".format(opts.random_state),
    opts.bowtie2_options,
    ")",
    # Convert to sorted BAM
    "|",
    "(",
    os.environ["samtools"],
    "sort",
    "--threads {}".format(opts.n_jobs),
    "--reference {}".format(opts.reference_assembly),
    "-T {}".format(os.path.join(directories["tmp"], "samtools_sort")),
    ">",
    output_filepaths[0],
    ")",
    ]
    return cmd

# Bowtie2
def get_coverage_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]

    # Command
    cmd = [
    "(",
    os.environ["jgi_summarize_bam_contig_depths"],
    "--outputDepth",
    output_filepaths[0],
    "--noIntraDepthVariance",
    opts.coverage_options,
    input_filepaths[0],
    ")",
    ]
    return cmd

# Symlink
def get_symlink_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    cmd = ["("]
    for filepath in input_filepaths:
        cmd.append("ln -f -s {} {}".format(os.path.realpath(filepath), os.path.realpath(output_directory)))
        cmd.append("&&")
    cmd[-1] = ")"
    return cmd

# ============
# Run Pipeline
# ============
# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """

    required_executables={
                "bowtie2",
                "samtools",
                "jgi_summarize_bam_contig_depths",
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
    pipeline = ExecutablePipeline(name="coverage.py", description=opts.name, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])

    # ==========
    # Bowtie2
    # ==========
    program = "bowtie2"
    # Add to directories
    output_directory = directories[("intermediate",  program)] = create_directory(os.path.join(directories["intermediate"], "{}_output".format(program)))

    # Info
    step = 1
    description = "Aligning reads to reference"

    # i/o
    input_filepaths = [opts.r1, opts.r2]

    output_filenames = ["mapped.sorted.bam"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_bowtie2_cmd(**params)
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
    # Coverage
    # ==========
    program = "coverage"
    # Add to directories
    output_directory = directories[("intermediate",  program)] = create_directory(os.path.join(directories["intermediate"], "{}_output".format(program)))

    # Info
    step = 2
    description = "Get coverage table"

    # i/o
    input_filepaths = output_filepaths

    output_filenames = ["coverage.tsv"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_coverage_cmd(**params)
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

       

    # =============
    # Symlink
    # =============
    program = "symlink"
    # Add to directories
    output_directory = directories["output"]

    # Info
    step = 3
    description = "Symlinking relevant output files"

    # i/o
    input_filepaths = [
        os.path.join(directories[("intermediate", "bowtie2")], "mapped.sorted.bam"),
        os.path.join(directories[("intermediate", "coverage")], "coverage.tsv"),

    ]

    output_filenames =  map(lambda fp: fp.split("/")[-1], input_filepaths)
    output_filepaths = list(map(lambda fn:os.path.join(directories["output"], fn), output_filenames))
        # Prodigal
        # os.path.join(directories["output"], "*"),
    
    params = {
    "input_filepaths":input_filepaths,
    "output_filepaths":output_filepaths,
    "output_directory":output_directory,
    "opts":opts,
    "directories":directories,
    }

    cmd = get_symlink_cmd(**params)
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

    return pipeline

# Configure parameters
def configure_parameters(opts, directories):
    if not opts.merge:
        assert opts.name is not None, "Must include --name"
        assert opts.r1 is not None, "Must include --r1 reads"
        assert opts.r2 is not None, "Must include --r2 reads"
        assert opts.r1 != opts.r2, "You probably mislabeled the input files because `r1` should not be the same as `r2`: {}".format(opts.r1)
        assert opts.reference_assembly is not None, "Must include --reference_assembly"
        if not opts.bowtie2_index:
            # Do a check that the index files exist
            opts.bowtie2_index = opts.reference_assembly

    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "\n".join([
         "{} -1 <r1.fq> -2 <r2.fq> -R <reference_assembly> -I <bowtie2_index> -o <output_directory>".format(__program__),
         "{} --merge <output_directory from above> > merged_coverage.tsv".format(__program__),
    ])
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_required = parser.add_argument_group('Required I/O arguments')
    parser_required.add_argument("-1","--r1", type=str, help = "path/to/r1.fq")
    parser_required.add_argument("-2","--r2", type=str, help = "path/to/r2.fq")
    parser_required.add_argument("-R", "--reference_assembly", type=str, help = "path/to/reference.fasta" )
    parser_required.add_argument("-I","--bowtie2_index", type=str, help = "path/to/index")
    parser_required.add_argument("-n", "--name", type=str, help="Name of sample")
    parser_required.add_argument("-o","--project_directory", type=str, default="./coverage_output", help = "path/to/project_directory [Default: ./coverage_output]")
    parser_required.add_argument("--merge", type=str, help="path/to/directory to merge coverage files.  Outputs to stdout")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Bowtie2
    parser_bowtie2 = parser.add_argument_group('Bowtie2 arguments')
    parser_bowtie2.add_argument("--bowtie2_options", type=str, default="", help="bowtie2 | More options (e.g. --arg 1 ) [Default: '']")

    # Coverage
    parser_coverage = parser.add_argument_group('Coverage arguments')
    parser_coverage.add_argument("--coverage_options", type=str, default="", help="jgi_summarize_bam_contig_depths | More options (e.g. --arg 1 ) [Default: '']\nhttp://manpages.ubuntu.com/manpages/groovy/man1/jgi_summarize_bam_contig_depths.1.html")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    if not opts.merge:
        # Directories
        directories = dict()
        directories["project"] = create_directory(opts.project_directory)
        directories["sample"] = create_directory(os.path.join(directories["project"], opts.name))
        directories["preprocessing"] = create_directory(os.path.join(directories["sample"], "preprocessing"))
        directories["output"] = create_directory(os.path.join(directories["sample"], "output"))
        directories["log"] = create_directory(os.path.join(directories["sample"], "log"))
        directories["tmp"] = create_directory(os.path.join(directories["sample"], "tmp"))
        directories["checkpoints"] = create_directory(os.path.join(directories["sample"], "checkpoints"))
        directories["intermediate"] = create_directory(os.path.join(directories["sample"], "intermediate"))


        # Info
        print(format_header("coverage.py", "="), file=sys.stdout)
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
    else:
        filepaths = glob.glob(os.path.join(opts.merge, "*", "output", "coverage.tsv"))
        assert len(filepaths) > 0, "Not detecting any `coverage.tsv` files from {}".format(opts.merge)
        coverages = dict()
        lengths = dict()
        for fp in pv(filepaths, "Reading coverage files"):
            id_sample = fp.split("/")[-3]
            df = pd.read_csv(fp, sep="\t", index_col=0)
            coverages[id_sample] = df.iloc[:,-1]
            lengths.update(df["contigLen"].to_dict())
        df_coverage = pd.DataFrame(coverages).sort_index(axis=1)   #.fillna(0)
        average = df_coverage.mean(axis=1).to_frame("totalAvgDepth")
        lengths = pd.Series(lengths)[df_coverage.index].to_frame("contigLen").astype(int)
        df_output = pd.concat([lengths, average, df_coverage], axis=1)
        df_output.index.name = "contigName"
        df_output.to_csv(sys.stdout, sep="\t")

if __name__ == "__main__":
    main()
