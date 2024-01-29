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
__version__ = "2022.03.25"

# .............................................................................
# Notes
# .............................................................................
# * Make batch version that takes in a manifest file
# .............................................................................
# Primordial
# .............................................................................


# Assembly
def get_assembly_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [
        "(",
    os.environ[opts.spades_program],
    ]
    if "restart-from" not in str(opts.spades_options):
        cmd += [
        "-1 {}".format(input_filepaths[0]),
        "-2 {}".format(input_filepaths[1]),
        ]
    cmd += [
    "-o {}".format(output_directory),
    "--tmp-dir {}".format(os.path.join(directories["tmp"], "spades")),
    "--threads {}".format(opts.n_jobs),
    "--memory {}".format(opts.memory),
    opts.spades_options,
    ")",
    ]

    if opts.spades_program == "rnaspades.py":
        cmd += [ 
            "&&",
            os.environ["fasta_to_saf.py"],
            "-i",
            os.path.join(output_directory, "transcripts.fasta"),
            ">",
            os.path.join(output_directory, "transcripts.fasta.saf"),
        ]

    else:
        cmd += [ 
            "&&",
            os.environ["fasta_to_saf.py"],
            "-i",
            os.path.join(output_directory, "scaffolds.fasta"),
            ">",
            os.path.join(output_directory, "scaffolds.fasta.saf"),
        ]


    for fn in [ 
        "before_rr.fasta",
        "K*",
        "misc",
        "corrected",
        "first_pe_contigs.fasta",
        ]:
        cmd += [ 
            "&&",
            "rm -rf {}".format(os.path.join(output_directory, fn)),
        ]
    return cmd

# Bowtie2
def get_alignment_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [
    # Clear temporary directory just in case
    "rm -rf {}".format(os.path.join(directories["tmp"], "*")),
    "&&",

    # Bowtie2 Index
    "(",
    os.environ["bowtie2-build"],
    "--threads {}".format(opts.n_jobs),
    "--seed {}".format(opts.random_state),
    opts.bowtie2_index_options,
    input_filepaths[2], # Reference
    input_filepaths[2], # Index
    ")",

    "&&",

    # Bowtie2
    "(",
    os.environ["bowtie2"],
    "-x {}".format(input_filepaths[2]),
    "-1 {}".format(input_filepaths[0]),
    "-2 {}".format(input_filepaths[1]),
    "--threads {}".format(opts.n_jobs),
    # "--un-conc-gz {}".format(os.path.join(output_directory, "unmapped_%.{}.gz".format(unmapped_ext))),
    # "--un-gz {}".format(os.path.join(output_directory, "unmapped_singletons_%{}.gz".format(unmapped_ext))),
    "--seed {}".format(opts.random_state),
    "--no-unal",

    opts.bowtie2_options,
    ")",
    # Convert to sorted BAM
    "|",
    "(",
    os.environ["samtools"],
    "sort",
    "--threads {}".format(opts.n_jobs),
    "--reference {}".format(input_filepaths[2]),
    "-T {}".format(os.path.join(directories["tmp"], "samtools_sort")),
    ">",
    output_filepaths[0],
    ")",
    "&&",
    "(",
    os.environ["samtools"],
    "index",
    "-@ {}".format(opts.n_jobs),
    output_filepaths[0],
    ")",
    ]

    return cmd


# featureCounts
def get_featurecounts_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command

    # ORF-Level Counts
    cmd = [
    "mkdir -p {}".format(os.path.join(directories["tmp"], "featurecounts")),
    "&&",
    "(",
        os.environ["featureCounts"],
        # "-G {}".format(input_filepaths[0]),
        "-a {}".format(input_filepaths[1]),
        "-o {}".format(os.path.join(output_directory, "featurecounts.tsv")),
        "-F SAF",
        "--tmpDir {}".format(os.path.join(directories["tmp"], "featurecounts")),
        "-T {}".format(opts.n_jobs),
        opts.featurecounts_options,
        input_filepaths[2],
    ")",
        "&&",
    "gzip -f {}".format(os.path.join(output_directory, "featurecounts.tsv")),
        ]
    return cmd

# seqkit
def get_seqkit_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command

    # ORF-Level Counts
    cmd = [

        os.environ["seqkit"],
        "stats",
        "-a", 
        "-j {}".format(opts.n_jobs),
        "-T",
        "-b",
        os.path.join(directories[("intermediate","1__assembly")], "*.fasta"),
        "|",
        "gzip",
        ">",
        output_filepaths[0],
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
    accessory_scripts = {
                "fasta_to_saf.py"
                }

    required_executables={
 
                opts.spades_program,
                "bowtie2-build",
                "bowtie2",
                "samtools",
                "featureCounts",
                "seqkit",
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
    pipeline = ExecutablePipeline(name=__program__, description=opts.name, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])

    # ==========
    # Assembly
    # ==========
    
    step = 1

    # Info
    program = "assembly"
    program_label = "{}__{}".format(step, program)
    description = "Assembling paired-end reads"
    
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # i/o
    input_filepaths = [opts.forward_reads, opts.reverse_reads]
    if opts.spades_program == "rnaspades.py":
        output_filenames = ["transcripts.fasta", "transcripts.fasta.saf"]
    else:
        output_filenames = ["scaffolds.fasta", "scaffolds.fasta.saf"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_assembly_cmd(**params)
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

    # ==========
    # Alignment
    # ==========
    
    step = 2

    # Info
    program = "alignment"
    program_label = "{}__{}".format(step, program)
    description = "Aligning reads to assembly"

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # i/o
    input_filepaths = [
        opts.forward_reads,
        opts.reverse_reads,
    ]

    if opts.spades_program == "rnaspades.py":
        input_filepaths += [ 
            os.path.join(directories[("intermediate", "1__assembly")], "transcripts.fasta"),
        ]
    else:
        input_filepaths += [ 
            os.path.join(directories[("intermediate", "1__assembly")], "scaffolds.fasta"),
        ] 


    output_filenames = ["mapped.sorted.bam"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_alignment_cmd(**params)
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



    # ==========
    # featureCounts
    # ==========
    step = 3 

    # Info
    program = "featurecounts"
    program_label = "{}__{}".format(step, program)
    description = "Counting reads"

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # i/o
    if opts.spades_program == "rnaspades.py":
        input_filepaths = [ 
            os.path.join(directories[("intermediate", "1__assembly")], "transcripts.fasta"),
            os.path.join(directories[("intermediate", "1__assembly")], "transcripts.fasta.saf"),
        ]
    else:
        input_filepaths = [ 
            os.path.join(directories[("intermediate", "1__assembly")], "scaffolds.fasta"),
            os.path.join(directories[("intermediate", "1__assembly")], "scaffolds.fasta.saf"),
        ]

    input_filepaths += [ 
        os.path.join(directories[("intermediate", "2__alignment")], "mapped.sorted.bam"),
    ]

    output_filenames = ["featurecounts.tsv.gz"]
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

    # ==========
    # stats
    # ==========
    
    step = 4

    # Info
    program = "seqkit"
    program_label = "{}__{}".format(step, program)
    description = "Assembly statistics"

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # i/o
    input_filepaths = [
                    os.path.join(directories[("intermediate", "1__assembly")], "*.fasta"),

    ]



    output_filenames = ["seqkit_stats.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_seqkit_cmd(**params)
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


    # =============
    # Symlink
    # =============
    step = 5

    # Info
    program = "symlink"
    program_label = "{}__{}".format(step, program)
    description = "Symlinking relevant output files"

    # Add to directories
    output_directory = directories["output"]

    # i/o
    if opts.spades_program == "rnaspades.py":
        input_filepaths = [ 
            os.path.join(directories[("intermediate", "1__assembly")], "transcripts.fasta"),
            os.path.join(directories[("intermediate", "1__assembly")], "transcripts.fasta.*"),
        ]
    else:
        input_filepaths = [ 
            os.path.join(directories[("intermediate", "1__assembly")], "scaffolds.fasta"),
            os.path.join(directories[("intermediate", "1__assembly")], "scaffolds.fasta.*"),
        ] 
    input_filepaths += [
            os.path.join(directories[("intermediate", "2__alignment")], "mapped.sorted.bam"),
            os.path.join(directories[("intermediate", "2__alignment")], "mapped.sorted.bam.bai"),
            os.path.join(directories[("intermediate", "3__featurecounts")], "featurecounts.tsv.gz"),
            os.path.join(directories[("intermediate", "4__seqkit")], "seqkit_stats.tsv.gz"),
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

# Configure parameters
def configure_parameters(opts, directories):
    # os.environ[]

    assert opts.forward_reads != opts.reverse_reads, "You probably mislabeled the input files because `forward_reads` should not be the same as `reverse_reads`: {}".format(opts.forward_reads)
        # assert not bool(opts.unpaired_reads), "Cannot have --unpaired_reads if --forward_reads.  Note, this behavior may be changed in the future but it's an adaptation of interleaved reads."

    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -1 <forward_reads.fq> -2 <reverse_reads.fq> -n <name> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-1","--forward_reads", type=str, help = "path/to/forward_reads.fq")
    parser_io.add_argument("-2","--reverse_reads", type=str, help = "path/to/reverse_reads.fq")
    parser_io.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    parser_io.add_argument("-o","--project_directory", type=str, default="veba_output/assembly", help = "path/to/project_directory [Default: veba_output/assembly]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
    parser_utility.add_argument("--tmpdir", type=str, help="Set temporary directory")  #site-packges in future

    # Assembler
    parser_assembler = parser.add_argument_group('SPAdes arguments')
    parser_assembler.add_argument("-P", "--spades_program", type=str, default="metaspades.py", help="SPAdes | More options (e.g. --arg 1 ) [Default: 'metaspades.py']\nhttp://cab.spbu.ru/files/release3.11.1/manual.html")
    parser_assembler.add_argument("-m", "--memory", type=int, default=250, help="SPAdes | RAM limit in Gb (terminates if exceeded). [Default: 250]")
    parser_assembler.add_argument("--spades_options", type=str, default="", help="SPAdes | More options (e.g. --arg 1 ) [Default: '']\nhttp://cab.spbu.ru/files/release3.11.1/manual.html")

    # Aligner
    parser_aligner = parser.add_argument_group('Bowtie2 arguments')
    parser_aligner.add_argument("--bowtie2_index_options", type=str, default="", help="bowtie2-build | More options (e.g. --arg 1 ) [Default: '']")
    parser_aligner.add_argument("--bowtie2_options", type=str, default="", help="bowtie2 | More options (e.g. --arg 1 ) [Default: '']")

    # featureCounts
    parser_featurecounts = parser.add_argument_group('featureCounts arguments')
    parser_featurecounts.add_argument("--featurecounts_options", type=str, default="", help="featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Directories
    directories = dict()
    directories["project"] = create_directory(opts.project_directory)
    directories["sample"] = create_directory(os.path.join(directories["project"], opts.name))
    directories["output"] = create_directory(os.path.join(directories["sample"], "output"))
    directories["log"] = create_directory(os.path.join(directories["sample"], "log"))
    if not opts.tmpdir:
        opts.tmpdir = os.path.join(directories["sample"], "tmp")
    directories["tmp"] = create_directory(opts.tmpdir)
    directories["checkpoints"] = create_directory(os.path.join(directories["sample"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["sample"], "intermediate"))
    os.environ["TMPDIR"] = directories["tmp"]

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
