#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse
from collections import OrderedDict

import pandas as pd

from soothsayer_utils import *
from genopype import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.5.8"

# STAR
def get_star_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]

    # Command
    cmd = [

    os.environ["STAR"],
   "--genomeDir {}".format(opts.reference_index),
   "--readFilesIn {} {}".format(input_filepaths[0], input_filepaths[1]),
   "--outFileNamePrefix {}/".format(output_directory),
   "--runThreadN {}".format(opts.n_jobs),
   # "--outFilterScoreMinOverLread 0",
   # "--outFilterMatchNminOverLread 0",
   # "--outFilterMatchNmin {}".format(CONFIG_PARAMETERS["outFilterMatchNmin"]),
   "--outReadsUnmapped Fastx",
    opts.star_options,
    ]

    if opts.forward_reads.endswith(".gz"):
        cmd += [ 
        "--readFilesCommand zcat",
        ]
    
    cmd += [

        "&&",

    # Sort mapped BAM
    os.environ["samtools"],
    "sort",
    "-@ {}".format(opts.n_jobs),
    "-T {}".format(os.path.join(directories["tmp"], "samtools_sort")),
    os.path.join(output_directory, "Aligned.out.sam"),
    ">",
    os.path.join(output_directory, "mapped.sorted.bam"),

    # Index BAM
        "&&",

    os.environ["samtools"],
    "index",
    "-@",
    opts.n_jobs,
    output_filepaths[0],

    # Calculate coverage
        "&&",

    os.environ["samtools"],
    "coverage",
    output_filepaths[0],
    "|",
    "gzip",
    ">",
    "{}.coverage.tsv.gz".format(output_filepaths[0]),

    # Mapped R1 reads
        "&&",

    os.environ["samtools"],
    "view",
    os.path.join(output_directory, "mapped.sorted.bam"),
    "|",
    "cut -f1",
    "|",
    "sort -u",
    "|",
    # os.environ["pigz"], 
    # "-c -f -p {}".format(opts.n_jobs),
    "gzip",
    ">",
    os.path.join(output_directory, "mapped.reads.list.gz"),
    ]

    # Do something with unpaired reads 
    if opts.retain_unmapped_reads:
        cmd += [
        # Gzip unmapped R1 reads
            "&&",

        "gzip -c",
        os.path.join(output_directory, "Unmapped.out.mate1"),
        ">",
        os.path.join(output_directory, "unmapped_1.fastq.gz"),

        # Gzip unmapped R2 reads
            "&&",

        "gzip -c",
        os.path.join(output_directory, "Unmapped.out.mate2"),
        ">",
        os.path.join(output_directory, "unmapped_2.fastq.gz"),
        ]


    cmd += [

    # Compile stats
        "&&",
    os.environ["compile_star_statistics.py"],
    "-i {}".format(os.path.join(output_directory, "Log.final.out")),
    "-o {}".format(os.path.join(output_directory, "alignment_statistics.tsv.gz")),
    # "-n {}".format(opts.name),

    # Gzip logs
        "&&",

    "gzip",
    os.path.join(output_directory, "*.out"),
    os.path.join(output_directory, "SJ.out.tab"),



    # Remove sam file, temporary directory, and original unmapped fastqs
        "&&",

    "rm -rf {} {} {}".format(
        os.path.join(output_directory, "Aligned.out.sam"),
        os.path.join(output_directory, "_STARtmp"),
        os.path.join(output_directory, "*.mate*"),
        ),
    ]
    return cmd



# featureCounts
def get_featurecounts_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command

    # Gene-Level Counts
    cmd = [
    "mkdir -p {}".format(os.path.join(directories["tmp"], "featurecounts")),

        "&&",

    os.environ["featureCounts"],
    "-a {}".format(opts.reference_gtf),
    "-o {}".format(os.path.join(output_directory, "featurecounts.transcripts.tsv")),
    "-F GTF",
    "-g transcript_id",
    "--tmpDir {}".format(os.path.join(directories["tmp"], "featurecounts")),
    "-T {}".format(opts.n_jobs),
    opts.featurecounts_options,
    input_filepaths[0],
      
        "&&",

    "tail -n +3",
    os.path.join(output_directory, "featurecounts.transcripts.tsv"),
    "|",
    "cut -f1,7",
    "|",
    "gzip",
    ">",
    os.path.join(output_directory, "counts.transcripts.tsv.gz"),

        "&&",

    os.environ["featureCounts"],
    "-a {}".format(opts.reference_gtf),
    "-o {}".format(os.path.join(output_directory, "featurecounts.genes.tsv")),
    "-F GTF",
    "-g gene_id",
    "--tmpDir {}".format(os.path.join(directories["tmp"], "featurecounts")),
    "-T {}".format(opts.n_jobs),
    opts.featurecounts_options,
    input_filepaths[0],
      
        "&&",

    "tail -n +3",
    os.path.join(output_directory, "featurecounts.genes.tsv"),
    "|",
    "cut -f1,7",
    "|",
    "gzip",
    ">",
    os.path.join(output_directory, "counts.genes.tsv.gz"),

    ]

    if opts.retain_featurecounts:
        cmd += [
            "&&",
        "gzip {}".format(os.path.join(output_directory, "featurecounts.*.tsv")),
        ]
    else:
        cmd += [
            "&&",
        "rm {}".format(os.path.join(output_directory, "featurecounts.*.tsv")),
        ]

    return cmd

# Symlink
def get_symlink_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    cmd = ["("]
    for filepath in input_filepaths:
        cmd.append("ln -f -s -r {} {}".format(filepath, output_directory))
        cmd.append("&&")
    cmd[-1] = ")"
    return cmd


# ============
# Run Pipeline
# ============
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = {
        "compile_star_statistics.py",
    }

    required_executables={
                # "pigz",
                # 1
                "STAR",
                "samtools",
                # 2
                "featureCounts",
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
    # os.environ[]
    # forward_reads != reverse_reads
    assert opts.forward_reads != opts.reverse_reads, "You probably mislabeled the input files because `-1` should not be the same as `-2`: {}".format(opts.forward_reads)
    # assert opts.forward_reads.endswith(".gz") & opts.reverse_reads.endswith(".gz"), "--forward_reads and --reverse_reads must be gzipped"
    if opts.forward_reads.endswith(".gz"):
        assert opts.reverse_reads.endswith(".gz"), "If --forward_reads are gzipped then --reverse_reads must be gzipped"
    if opts.reverse_reads.endswith(".gz"):
        assert opts.forward_reads.endswith(".gz"), "If --reverse_reads are gzipped then --forward_reads must be gzipped"
    # Set environment variables
    add_executables_to_environment(opts=opts)




def create_pipeline(opts, directories, f_cmds):

    # .................................................................
    # Primordial
    # .................................................................
    # Commands file
    pipeline = ExecutablePipeline(name="STAR Mapping Pipeline", description=opts.name, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])

    # ==========
    # STAR
    # ==========
    step = 1
    program = "star"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    step = 1
    description = "Aligning reads to reference"

    # i/o
    input_filepaths = [opts.forward_reads, opts.reverse_reads]

    output_filenames = ["mapped.sorted.bam", "mapped.sorted.bam.bai", "mapped.sorted.bam.coverage.tsv.gz", "alignment_statistics.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_star_cmd(**params)
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
    input_filenames = ["mapped.sorted.bam"]
    input_filepaths = list(map(lambda filename: os.path.join(directories[("intermediate",  "1__star")], filename), input_filenames))

    output_filenames = ["counts.transcripts.tsv.gz", "counts.genes.tsv.gz"]
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
    program = "symlink"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories["output"]

    # Info
    step = 3
    description = "Symlinking relevant output files"


    # i/o
    input_filepaths = [
        os.path.join(directories[("intermediate", "1__star")], "mapped.sorted.bam"),
        os.path.join(directories[("intermediate", "1__star")], "mapped.sorted.bam.bai"),
        os.path.join(directories[("intermediate", "1__star")], "mapped.sorted.bam.coverage.tsv.gz"),
        os.path.join(directories[("intermediate", "1__star")], "mapped.reads.list.gz"),
        os.path.join(directories[("intermediate", "1__star")], "alignment_statistics.tsv.gz"),
        os.path.join(directories[("intermediate", "2__featurecounts")], "counts.transcripts.tsv.gz"),
        os.path.join(directories[("intermediate", "2__featurecounts")], "counts.genes.tsv.gz"),
    ]
    if opts.retain_unmapped_reads:
        input_filepaths += [
        os.path.join(directories[("intermediate", "1__star")], "unmapped_1.fastq.gz"),
        os.path.join(directories[("intermediate", "1__star")], "unmapped_2.fastq.gz"),
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
            validate_inputs=False,
            validate_outputs=False,
    )


    return pipeline



def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -1 <reads_1.fq[.gz]> -2 <reads_2.fq[.gz]> -n <name> -o <output_directory> -x <reference_directory>"    .format(__program__)
    epilog = "Copyright 2020 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_required = parser.add_argument_group('Required arguments')
    parser_required.add_argument("-1","--forward_reads", type=str, help = "path/to/forward_reads.fq[.gz]", required=True)
    parser_required.add_argument("-2","--reverse_reads", type=str, help = "path/to/reverse_reads.fq[.gz]", required=True)
    parser_required.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    parser_required.add_argument("-o","--project_directory", type=str, default="star_output", help = "path/to/project_directory [Default: star_output]")


    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Index
    parser_reference = parser.add_argument_group('Reference arguments')
    parser_reference.add_argument("-x", "--reference_index",type=str, required=True, help="path/to/star_index")
    # parser_reference.add_argument("-r", "--reference_fasta", type=str, required=True, help = "path/to/reference.fasta" )
    parser_reference.add_argument("-g", "--reference_gtf",type=str, required=True, help="path/to/reference.saf. If not provided then it is created from --reference_gtf")


    # STAR
    parser_star = parser.add_argument_group('STAR arguments')
    parser_star.add_argument("--retain_unmapped_reads", default=1, type=int, help = "Retain reads that do not map to reference. 0=No, 1=yes [Default: 1]") 
    parser_star.add_argument("--star_options", type=str, default="", help="STAR | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/alexdobin/STAR")

    # featureCounts
    parser_featurecounts = parser.add_argument_group('featureCounts arguments')
    parser_featurecounts.add_argument("--retain_featurecounts", default=0, type=int, help = "Retain feature counts output table (a slimmer version is output regardless). 0=No, 1=yes [Default: 0]") 
    parser_featurecounts.add_argument("--featurecounts_options", type=str, default="", help="featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Potential names

    # Directories
    directories = dict()
    directories["project"] = create_directory(opts.project_directory)
    directories["sample"] = create_directory(os.path.join(directories["project"], opts.name))
    directories["output"] = create_directory(os.path.join(directories["sample"], "output"))
    directories["log"] = create_directory(os.path.join(directories["sample"], "log"))
    directories["tmp"] = create_directory(os.path.join(directories["sample"], "tmp"))
    directories["checkpoints"] = create_directory(os.path.join(directories["sample"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["sample"], "intermediate"))

    # Configure parameters


    # Info
    print(format_header("STAR Pipeline", "="), file=sys.stdout)
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
