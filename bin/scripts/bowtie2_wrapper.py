#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
from collections import OrderedDict, defaultdict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.5.15"


# Bowtie2
def get_bowtie2_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    # Is fasta or fastq?
    ignore_quals = False
    unmapped_ext = "fastq"
    fp, ext = os.path.splitext(input_filepaths[0])
    if ext in {".gz",".bz2"}:
        fp, ext = os.path.splitext(fp)
        ext = ext[1:]
    if ext in {"fa", "fasta", "fna"}:
        unmapped_ext = "fasta"
        ignore_quals = True
    
    # Command
    cmd = [
    # Clear temporary directory just in case
    "rm -rf {}".format(os.path.join(directories["tmp"], "*")),
    "&&",
    # Bowtie2
    "(",
    os.environ["bowtie2"],
   "-x {}".format(opts.reference_index),
   "-1 {}".format(input_filepaths[0]),
   "-2 {}".format(input_filepaths[1]),
   "--threads {}".format(opts.n_jobs),
    "--seed {}".format(opts.random_state),
    "--met-file {}".format(os.path.join(output_directory, "metrics.txt")),
    "--no-unal",
    ]

    # Do something with unpaired reads 
    if opts.retain_unmapped_reads:
        cmd += [
        "--un-conc-gz {}".format(os.path.join(output_directory, "unmapped_%.{}.gz".format(unmapped_ext))),
        # "--un-gz {}".format(os.path.join(output_directory, "unmapped_singletons_%.{}.gz".format(unmapped_ext))),
        ]
    if ignore_quals:
        cmd.append("-f")

    cmd += [
        opts.bowtie2_options,
        ")",
    ]

    # Convert to sorted BAM
    cmd += [
        "|",
        "(",
        os.environ["samtools"],
        "sort",
        "--threads {}".format(opts.n_jobs),
        "--reference {}".format(opts.reference_fasta),
        "-T {}".format(os.path.join(directories["tmp"], "samtools_sort")),
        ">",
        output_filepaths[0],
        ")",

        # Index BAM
        "&&",
        "(",
        os.environ["samtools"],
        "index",
        "-@",
        opts.n_jobs,
        output_filepaths[0],
        ")",

        # Calculate coverage
        "&&",
        "(",
        os.environ["samtools"],
        "coverage",
        output_filepaths[0],
        "|",
        "gzip",
        ">",
        "{}.coverage.tsv.gz".format(output_filepaths[0]),
        ")",
    ]


    cmd += [ 
        "&&",
        "gzip {}".format(os.path.join(output_directory, "metrics.txt")),
    ]
    return cmd



# featureCounts
def get_featurecounts_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command

    # ORF-Level Counts
    cmd = [
    "mkdir -p {}".format(os.path.join(directories["tmp"], "featurecounts")),

        "&&",

    os.environ["featureCounts"],
    "-a {}".format(opts.reference_saf),
    "-o {}".format(os.path.join(output_directory, "featurecounts.tsv")),
    "-F SAF",
    "--tmpDir {}".format(os.path.join(directories["tmp"], "featurecounts")),
    "-T {}".format(opts.n_jobs),
    "-p --countReadPairs",
    opts.featurecounts_options,
    input_filepaths[0],
      
        "&&",

    "tail -n +3",
    os.path.join(output_directory, "featurecounts.tsv"),
    "|",
    "cut -f1,7",
    "|",
    "gzip",
    ">",
    os.path.join(output_directory, "counts.tsv.gz"),
    ]
    if opts.retain_featurecounts:
        cmd += [
            "&&",
        "gzip {}".format(os.path.join(output_directory, "featurecounts.tsv")),
        ]
    else:
        cmd += [
            "&&",
        "rm {}".format(os.path.join(output_directory, "featurecounts.tsv")),
        ]

    return cmd

# Symlink
def get_symlink_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    cmd = [
    "DST={}; (for SRC in {}; do SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST; done)".format(
        output_directory,
        " ".join(input_filepaths), 
        )
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

    accessory_scripts = {
        "fasta_to_saf.py",
        }

    required_executables={
                "bowtie2",
                "featureCounts",
                "samtools",
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

# Pipeline
def create_pipeline(opts, directories, f_cmds):

    # .................................................................
    # Primordial
    # .................................................................
    # Commands file
    pipeline = ExecutablePipeline(name=__program__,  f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])

    
    # ==========
    # Bowtie2
    # ==========
    step = 1

    # Info
    program = "bowtie2"
    program_label = "{}__{}".format(step, program)
    description = "Aligning reads to reference"

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # i/o
    input_filepaths = [opts.forward_reads, opts.reverse_reads]

    output_filenames = ["mapped.sorted.bam", "mapped.sorted.bam.bai", "mapped.sorted.bam.coverage.tsv.gz"]
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
    # featureCounts
    # ==========
    step = 2

    # Info
    program = "featurecounts"
    program_label = "{}__{}".format(step, program)
    description = "Counting reads"

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # i/o
    input_filepaths = output_filepaths

    output_filenames = ["counts.tsv.gz"]

    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_featurecounts_cmd(**params)
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
    step = 3

    # Info
    program = "symlink"
    program_label = "{}__{}".format(step, program)
    description = "Symlinking relevant output files"

    # Add to directories
    output_directory = directories["output"]

    # i/o
    input_filepaths = [
        os.path.join(directories[("intermediate", "1__bowtie2")], "mapped.sorted.bam"),
        os.path.join(directories[("intermediate", "1__bowtie2")], "mapped.sorted.bam.bai"),
        os.path.join(directories[("intermediate", "1__bowtie2")], "mapped.sorted.bam.coverage.tsv.gz"),
    ]
    if opts.retain_unmapped_reads:
         input_filepaths += [
            os.path.join(directories[("intermediate", "1__bowtie2")], "unmapped_*.gz"),
        ]

    input_filepaths += [ 
        os.path.join(directories[("intermediate", "2__featurecounts")], "counts.tsv.gz"),
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
            validate_outputs=True,
    )
    return pipeline

# Configure parameters
def configure_parameters(opts, directories):

    # Set environment variables
    add_executables_to_environment(opts=opts)

    # If --reference_fasta isn't provided then set it to the --reference_index
    if opts.reference_fasta is None:
        opts.reference_fasta =  opts.reference_index
        
    if not opts.reference_saf:
        if os.path.exists(os.path.join(directories["output"], "reference.saf")):
            print(" * No --reference_saf provided.  Loading SAF file from the following location: {}".format(os.path.join(directories["output"], "reference.saf")), file=sys.stdout)
        else:
            print(" * No --reference_saf provided.  Creating one from --reference_fasta in the following location: {}".format(os.path.join(directories["output"], "reference.saf")), file=sys.stdout)
            cmd = [ 
                os.environ["fasta_to_saf.py"],
                "-i {}".format(opts.reference_fasta),
                ">",
                os.path.join(directories["output"], "reference.saf"),
            ]
            Command(cmd).run() 
            opts.reference_saf = os.path.join(directories["output"], "reference.saf")




def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -1 <reads_1.fq[.gz]> -2 <reads_2.fq[.gz]> -n <name> -o <output_directory> -x <reference_directory>"    .format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-1","--forward_reads", type=str, help = "path/to/reads_1.fastq[.gz]", required=True)
    parser_io.add_argument("-2","--reverse_reads", type=str, help = "path/to/reads_2.fastq[.gz]", required=True)
    parser_io.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    parser_io.add_argument("-o","--project_directory", type=str, default="bowtie2_output", help = "path/to/project_directory [Default: bowtie2_output]")

    parser_reference = parser.add_argument_group('Reference arguments')
    parser_reference.add_argument("-x", "--reference_index",type=str, required=True, help="path/to/bowtie2_index")
    parser_reference.add_argument("-r", "--reference_fasta", type=str, required=False, help = "path/to/reference.fasta. If not provided then it is set to the --reference_index" )
    parser_reference.add_argument("-s", "--reference_saf",type=str, required=False, help="path/to/reference.saf. If not provided then it is created from --reference_fasta")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv. Must have at least 2 columns [name, executable] [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=int, help = "Restart from a particular checkpoint")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Bowtie2
    parser_bowtie2 = parser.add_argument_group('Bowtie2 arguments')
    parser_bowtie2.add_argument("--retain_unmapped_reads", default=1, type=int, help = "Retain reads that do not map to reference. 0=No, 1=yes [Default: 1]") 
    parser_bowtie2.add_argument("--bowtie2_options", type=str, default="", help="Bowtie2 | More options (e.g. --arg 1 ) [Default: '']\nhttp://bowtie-bio.sourceforge.net/bowtie2/manual.shtml")

    # featureCounts
    parser_featurecounts = parser.add_argument_group('featureCounts arguments')
    parser_featurecounts.add_argument("--retain_featurecounts", default=0, type=int, help = "Retain feature counts output table (a slimmer version is output regardless). 0=No, 1=yes [Default: 0]") 
    parser_featurecounts.add_argument("--featurecounts_options", type=str, default="", help="featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/")

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

    directories["project"] = create_directory(opts.project_directory)
    directories["sample"] = create_directory(os.path.join(directories["project"], opts.name))
    directories["output"] = create_directory(os.path.join(directories["sample"], "output"))
    directories["log"] = create_directory(os.path.join(directories["sample"], "log"))
    directories["tmp"] = create_directory(os.path.join(directories["sample"], "tmp"))
    directories["checkpoints"] = create_directory(os.path.join(directories["sample"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["sample"], "intermediate"))
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
