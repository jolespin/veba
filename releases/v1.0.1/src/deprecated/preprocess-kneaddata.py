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
__version__ = "2021.12.21"

DB_RIBOKMERS="/usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/bbmap_env/ribokmers.fa.gz"
DB_HUMAN="/usr/local/scratch/CORE/jespinoz/db/genomes/human/GRCh38.p13/"

# .............................................................................
# Notes
# .............................................................................
# * Make batch version that takes in a manifest file
# .............................................................................
# Primordial
# .............................................................................


# Kneaddata
def get_kneaddata_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]

    # Command
    cmd = [
    "(",
    os.environ["kneaddata"],
    "--input {}".format(input_filepaths[0]),
    "--input {}".format(input_filepaths[1]),
    {True:"--reference-db {}".format(opts.kneaddata_contamination_db), False:""}[bool(opts.kneaddata_contamination_db)],
    "--output {}".format(output_directory),
    "--log {}".format(os.path.join(output_directory, "kneaddata.log")),
    "--threads {}".format(opts.n_jobs),
    "--output-prefix kneaddata",
    '--bowtie2-options="--seed {}"'.format(opts.random_state, {True:" {}".format(opts.kneaddatabowtie2_options), False:""}[bool(opts.kneaddatabowtie2_options)]), # Work around to ensure there's not an extra space for the options
    opts.kneaddata_options,
    ")",
    "&&",
    "(",
    os.environ["repair.sh"],
    ]
    # Contamination
    if bool(opts.kneaddata_contamination_db):
        cmd += [
        "in1={}".format(os.path.join(output_directory, "kneaddata_paired_1.fastq")),
        "in2={}".format(os.path.join(output_directory, "kneaddata_paired_2.fastq")),
        ]
    # No Contamination
    else:
        cmd += [
        "in1={}".format(os.path.join(output_directory, "kneaddata.trimmed.1.fastq")),
        "in2={}".format(os.path.join(output_directory, "kneaddata.trimmed.2.fastq")),
        ]

    cmd += [
        "out1={}".format(os.path.join(output_directory,"kneaddata_repaired_1.fastq.gz")),
        "out2={}".format(os.path.join(output_directory, "kneaddata_repaired_2.fastq.gz")),
        "outs={}".format(os.path.join(output_directory, "kneaddata_repaired_singletons.fastq.gz")),
        "overwrite=t",
        ")",
        "&&",
        "rm -f {} {}".format(
            os.path.join(output_directory, "kneaddata_paired_*"),
            os.path.join(output_directory, "kneaddata_unmatched_*"),
        )

    ]

    # Gzip
    cmd += [
        "&&",
        "pigz -f -p {} {}".format(opts.n_jobs, os.path.join(output_directory, "*.fastq")),
    ]

    # Get counts
    cmd += [ 
        "&&",
        os.environ["seqkit"],
        "stats",
        "-T",
        "--threads {}".format(opts.n_jobs),

        os.path.join(output_directory, "*.fastq.gz"),
        ">",
        output_filepaths[2],
    ]

    if opts.remove_contamination:
        cmd += [
            "&&",
            "rm -f {}".format(
                os.path.join(output_directory, "*_contam*"),
            )
        ]
    if opts.remove_trimmed:
        cmd += [
            "&&",
            "rm -f {}".format(
            # os.path.join(output_directory, "kneaddata_paired*"),
            os.path.join(output_directory, "*trimmed*"),
            ),
        ]


    return cmd

# BBDuk
def get_bbduk_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    cmd = [
        os.environ["bbduk.sh"],
        "fastawrap=1000",
        "overwrite=t",
        "threads={}".format(opts.n_jobs),
        "in1={}".format(input_filepaths[0]),
        "in2={}".format(input_filepaths[1]),
        "refstats={}".format(output_filepaths[0]),
        "stats={}".format(output_filepaths[1]),
        "k={}".format(opts.kmer_size),
        "minlen={}".format(opts.minimum_read_length),
        opts.bbduk_options,
    ]

    if bool(opts.bbduk_db):
        cmd += [ 
            "ref={}".format(opts.bbduk_db),
        ]            


    if not opts.report_only:
        cmd +=  [ 
        "out1={}".format(output_filepaths[2]),
        "out2={}".format(output_filepaths[3]),
        "outm1={}".format(output_filepaths[4]),
        "outm2={}".format(output_filepaths[5]),
        ]
    if (not opts.report_only):
        # Get counts
        cmd += [ 
            "&&",
            os.environ["seqkit"],
            "stats",
            "-T",
            "--threads {}".format(opts.n_jobs),
            os.path.join(output_directory, "*.fastq.gz"),
            ">",
            output_filepaths[-1],
        ]
        if opts.remove_contamination:
            cmd += [ 
                "&&",
                "rm {}".format(output_filepaths[4]),
                "&&",
                "rm {}".format(output_filepaths[5]),
            ]
    return cmd


# Symlink
def get_symlink_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    cmd = ["("]
    for filepath in input_filepaths:
        # cmd.append("ln -f -s {} {}".format(os.path.realpath(filepath), os.path.realpath(output_directory)))
        cmd.append("ln -f -s {} {}".format(os.path.realpath(filepath), output_directory))
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
                "repair.sh",
                "bbduk.sh",
                "kneaddata",
                "seqkit",
    }

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
    accessory_scripts = []
    for name in accessory_scripts:
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
    pipeline = ExecutablePipeline(name=__program__, description=opts.name, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])

    # =========
    # Kneaddata
    # =========
    step = 0
    program = "kneaddata"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Quality trimming, removing contaminated reads, and optimizing file compression"
    # i/o
    # if not opts.bypass_decontamination:
    input_filepaths = [opts.r1, opts.r2]
    # else:
        # assert opts.single_reads is not None, "If `bypass_decontamination` then the reads must be single-ended"
        # input_filepaths = [os.path.join(output_directory, "reads.subsampled.clumpify.fastq.gz")]
        # os.path.join(output_directory,"kneaddata_repaired_1.fastq.gz")
    output_filenames = ["kneaddata_repaired_1.fastq.gz", "kneaddata_repaired_2.fastq.gz", "kneaddata_seqkit-stats.tsv"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    # Parameters
    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }
    # Command
    cmd = get_kneaddata_cmd(**params)
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

    # =========
    # BBDuk
    # =========
    step = 1
    program = "bbduk"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))
    description = "Removing contamination via BBDuk"
    # i/o
    # if not opts.bypass_decontamination:
    input_filepaths = output_filepaths
    # else:
        # assert opts.single_reads is not None, "If `bypass_decontamination` then the reads must be single-ended"
        # input_filepaths = [os.path.join(output_directory, "reads.subsampled.clumpify.fastq.gz")]
        # os.path.join(output_directory,"kneaddata_repaired_1.fastq.gz")
    output_filenames = ["bbduk_refstats.txt","bbduk_stats.txt"]
    if not opts.report_only:
        output_filenames += [ 
            "kneaddata_repaired_1.bbduk.fastq.gz",
            "kneaddata_repaired_2.bbduk.fastq.gz",
        ]
    if not opts.remove_contamination:
        output_filenames += [ 
            "contamination_1.fastq.gz",
            "contamination_2.fastq.gz",
        ]
    if not opts.report_only:
        output_filenames += [
        "bbduk_seqkit-stats.tsv",
        ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    # Parameters
    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }
    # Command
    cmd = get_bbduk_cmd(**params)
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
    step = 2
    program = "symlink"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories["output"] 
    description = "Symlinking relevant output files"

    # i/o
    input_filepaths = [
        os.path.join(directories[("intermediate", "0__kneaddata")], "kneaddata_repaired_1.fastq.gz"),
        os.path.join(directories[("intermediate", "0__kneaddata")], "kneaddata_repaired_2.fastq.gz"),
        os.path.join(directories[("intermediate", "0__kneaddata")], "kneaddata_seqkit-stats.tsv"),

        os.path.join(directories[("intermediate", "1__bbduk")], "bbduk_refstats.txt"),
        os.path.join(directories[("intermediate", "1__bbduk")], "bbduk_stats.txt"),


    ]
    
    if not opts.report_only: 
        input_filepaths += [
        os.path.join(directories[("intermediate", "1__bbduk")], "kneaddata_repaired_1.bbduk.fastq.gz"),
        os.path.join(directories[("intermediate", "1__bbduk")], "kneaddata_repaired_2.bbduk.fastq.gz"),
        os.path.join(directories[("intermediate", "1__bbduk")], "bbduk_seqkit-stats.tsv"),

    ]

    output_filenames =  map(lambda fp: fp.split("/")[-1], input_filepaths)
    output_filepaths = list(map(lambda fn:os.path.join(directories["output"], fn), output_filenames))

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

    if bool(opts.r1):
        assert opts.r1 != opts.r2, "You probably mislabeled the input files because `r1` should not be the same as `r2`: {}".format(opts.r1)
        # assert not bool(opts.unpaired_reads), "Cannot have --unpaired_reads if --r1.  Note, this behavior may be changed in the future but it's an adaptation of interleaved reads."

    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -1 <r1.fq> -2 <r2.fq> -n <name> -o <output_directory> --kneaddata_contamination_db <database/>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-1","--r1", type=str, help = "path/to/r1.fq")
    parser_io.add_argument("-2","--r2", type=str, help = "path/to/r2.fq")
    parser_io.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    parser_io.add_argument("-o","--project_directory", type=str, default="preprocess", help = "path/to/project_directory [Default: veba_output/assembly]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--remove_contamination", action="store_true", help = "Remove contaminant fastq (kneaddata|bbduk)") #!
    parser_utility.add_argument("--remove_trimmed", action="store_true", help = "Remove trimmed fastq (kneaddata)") #!
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))


    # Kneaddata
    parser_kneaddata = parser.add_argument_group('Kneaddata arguments')
    parser_kneaddata.add_argument("-d", "--kneaddata_contamination_db", type=str, help="Kneaddata | path/to/contamination_database\nFor human at JCVI, use the following: {}".format(DB_HUMAN))
    parser_kneaddata.add_argument("--kneaddata_options", type=str, default="", help="Kneaddata | More options (e.g. --arg 1 ) [Default: ''] | https://bitbucket.org/biobakery/kneaddata/wiki/Home")
    parser_kneaddata.add_argument("--kneaddatabowtie2_options", type=str, default="", help="Bowtie2 | More options (e.g. --arg 1 ) [Default: '']\nhttp://bowtie-bio.sourceforge.net/bowtie2/manual.shtml")

    # BBDuk
    parser_bbduk = parser.add_argument_group('BBDuk arguments')
    parser_bbduk.add_argument("--bbduk_db", type=str, default=DB_RIBOKMERS, help="BBDuk | path/to/contamination_database [Default: '']\nFor ribokmers at JCVI, use the following: {}".format(DB_RIBOKMERS))
    parser_bbduk.add_argument("-k", "--kmer_size", type=int, default=31, help="BBDuk | k-mer size [Default: 31]")
    parser_bbduk.add_argument("-m", "--minimum_read_length", type=int, default=50, help="BBDuk | Minimum read length [Default: 50]")
    parser_bbduk.add_argument("--report_only", action="store_true", help = "Report only if decontaminated reads should not be output") #!
    parser_bbduk.add_argument("--bbduk_options", type=str, default="", help="BBDuk | More options (e.g., --arg 1) [Default: '']")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

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
    print(format_header(__program__, "="), file=sys.stdout)
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
