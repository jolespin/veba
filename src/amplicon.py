#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
from collections import OrderedDict, defaultdict
import numpy as np
import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.5.8"

# Reads archive
def get_reads_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    cmd = [
        os.environ["qiime"],
        "tools",
        "import",
        "--type 'SampleData[PairedEndSequencesWithQuality]'",
        "--input-path {}".format(opts.reads_table),
        "--output-path {}".format(os.path.join(output_directory, "demultiplexed_reads.qza")),
        "--input-format {}".format(opts.input_reads_format),
        "&&",
        os.environ["qiime"],
        "demux",
        "summarize", 
        "--i-data {}".format(os.path.join(output_directory, "demultiplexed_reads.qza")),
        "--o-visualization {}".format(os.path.join(output_directory, "demultiplexed_reads.qzv")),
        "&&",
        "FORWARD_READS=$(tail -n +2 {} | cut -f2 )".format(opts.reads_table),
        "&&",
        os.environ["fastq_position_statistics.py"],
        "-f q=0.25",
        "-o {}".format(os.path.join(output_directory, "forward_reads.position-specific_statistics.tsv")),
        "${FORWARD_READS}",
        "&&",
        "REVERSE_READS=$(tail -n +2 {} | cut -f3 )".format(opts.reads_table),
        "&&",
        os.environ["fastq_position_statistics.py"],
        "-f q=0.25",
        "-o {}".format(os.path.join(output_directory, "reverse_reads.position-specific_statistics.tsv")),
        "${REVERSE_READS}",
        
        ]
    return cmd

# Trim detection
def get_trim_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    cmd = [
    ]
    if opts.forward_trim is not None:
        cmd += [
        "mkdir -p {}".format(os.path.join(output_directory, "forward_reads_trimming_suggestion")),
        "&&",
        "echo {} > {}".format(opts.forward_trim, os.path.join(output_directory, "forward_reads_trimming_suggestion", "position_to_trim.txt")),
        ]
    else:
        cmd += [
        os.environ["determine_trim_position.py"],
        "-i {}".format(input_filepaths[0]),
        "-o {}".format(os.path.join(output_directory, "forward_reads_trimming_suggestion")),
        "--minimum_quality {}".format(opts.minimum_quality),
        "--minimum_length {}".format(opts.minimum_length),
        "--window_size {}".format(opts.window_size),
        "--plot_title {}".format(opts.reads_table),
        ]
        if opts.maximum_average_loss is not None:
            cmd += ["--maximum_average_loss {}".format(opts.maximum_average_loss)]

    if opts.reverse_trim is not None:
        cmd += [
        "mkdir -p {}".format(os.path.join(output_directory, "reverse_reads_trimming_suggestion")),
        "&&",
        "echo {} > {}".format(opts.reverse_trim, os.path.join(output_directory, "reverse_reads_trimming_suggestion", "position_to_trim.txt")),
        ]
    else:
        cmd += [
        "&&",
        os.environ["determine_trim_position.py"],
        "-i {}".format(input_filepaths[1]),
        "-o {}".format(os.path.join(output_directory, "reverse_reads_trimming_suggestion")),
        "--minimum_quality {}".format(opts.minimum_quality),
        "--minimum_length {}".format(opts.minimum_length),
        "--window_size {}".format(opts.window_size),
        "--plot_title {}".format(opts.reads_table),
        ]
        if opts.maximum_average_loss is not None:
            cmd += ["--maximum_average_loss {}".format(opts.maximum_average_loss)]
    return cmd

# DADA2
def get_dada2_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
 
    # Command
    cmd = [
        "FORWARD_TRIM_POSITION=$(cat {})".format(input_filepaths[0]),
            "&&",
        "REVERSE_TRIM_POSITION=$(cat {})".format(input_filepaths[1]),
            "&&",
        os.environ["qiime"],
        "dada2",
        "denoise-paired",
        "--i-demultiplexed-seqs {}".format(input_filepaths[2]),
        "--p-trunc-len-f ${FORWARD_TRIM_POSITION}",
        "--p-trunc-len-r ${REVERSE_TRIM_POSITION}",
        "--o-table {}".format(output_filepaths[0]),
        "--o-representative-sequences {}".format(output_filepaths[1]),
       "--o-denoising-stats {}".format(output_filepaths[2]),
       "--p-n-threads {}".format(opts.n_jobs),
       "--p-min-overlap {}".format(opts.minimum_overlap),
       opts.dada2_options,
    ]
    return cmd

# Classify
def get_classify_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    cmd = [
        os.environ["qiime"],
        "feature-classifier",
        "classify-sklearn",
        "--i-reads {}".format(input_filepaths[0]),
        "--i-classifier {}".format(input_filepaths[1]),
        "--o-classification {}".format(output_filepaths[0]),
        "--p-n-jobs {}".format(opts.n_jobs),
    ]
    return cmd

# Phylogeny
def get_phylogeny_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # qiime phylogeny align-to-tree-mafft-fasttree \
    # --i-sequences rep-seqs.qza \
    # --o-alignment aligned-rep-seqs.qza \
    # --o-masked-alignment masked-aligned-rep-seqs.qza \
    # --o-tree unrooted-tree.qza \
    # --o-rooted-tree rooted-tree.qza

    # Command
    cmd = [
        os.environ["qiime"],
        "phylogeny",
        opts.phylogeny_mode,
        "--p-n-threads {}".format(opts.n_jobs),
        "--i-sequences {}".format(input_filepaths[0]),
        "--o-tree {}".format(output_filepaths[0]),
        "--o-rooted-tree {}".format(output_filepaths[1]),
        "--o-alignment {}".format(output_filepaths[2]),
        "--o-masked-alignment {}".format(output_filepaths[3]),
    ]
    return cmd

# Conversion
def get_conversion_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [
"""
for FP in {}; do
    {} tools export --input-path $FP --output-path {};
    done

""".format(" ".join(input_filepaths), os.environ["qiime"], output_directory),
    os.environ["biom"],
    "convert",
    "-i {}".format(os.path.join(output_directory, "feature-table.biom")),
    "-o {}".format(os.path.join(output_directory, "feature-table.tsv")),
    "--to-tsv",
        "&&",
    os.environ["replace_fasta_descriptions.py"],
    "-f {}".format(os.path.join(output_directory, "dna-sequences.fasta")),
    "-d {}".format(os.path.join(output_directory, "taxonomy.tsv")),
    "-o {}".format(os.path.join(output_directory, "dna-sequences.with_taxonomy.fasta")),

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
                "determine_trim_position.py",
                "fastq_position_statistics.py",
                "replace_fasta_descriptions.py",
                }

    required_executables={
                "qiime",
                "biom",
     } | accessory_scripts

    executables = dict()
    for name in required_executables:
        executables[name] = os.path.join(os.environ["CONDA_PREFIX"], "bin", name)

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
    pipeline = ExecutablePipeline(name=__program__, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])
    
    # -------------------------

    step = 1

    # Info
    program = "reads_archive"
    program_label = "{}__{}".format(step, program)
    description = "Create reads archive"

    # Add to directories
    output_directory = directories["intermediate"] 

    # i/o
    input_filepaths = [
        opts.reads_table,
    ]

   

    output_filenames = ["demultiplexed_reads.qza", "demultiplexed_reads.qzv", "forward_reads.position-specific_statistics.tsv","reverse_reads.position-specific_statistics.tsv"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_reads_cmd(**params)
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

    # -------------------------
    step = 2

    # Info
    program = "trimming"
    program_label = "{}__{}".format(step, program)
    description = "Determine trim position"

    # Add to directories
    output_directory = directories["intermediate"] 


    # i/o
    input_filepaths = [
        os.path.join(directories["intermediate"], "forward_reads.position-specific_statistics.tsv"),
        os.path.join(directories["intermediate"], "reverse_reads.position-specific_statistics.tsv"),
    ]

    output_filenames = ["forward_reads_trimming_suggestion/position_to_trim.txt", "reverse_reads_trimming_suggestion/position_to_trim.txt"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_trim_cmd(**params)
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

    continue_pipeline = True
    if ((opts.forward_trim is None) or (opts.reverse_trim is None)):
        if opts.inspect_trim_regions:
            continue_pipeline = False
            print("--inspect_trim_regions was provided. Please inspect trim regions for reads.", file=sys.stdout)

    if continue_pipeline:
        # -------------------------
        step = 3

        # Info
        program = "dada2"
        program_label = "{}__{}".format(step, program)
        description = "Identiy ASV reference sequences"

        # Add to directories
        output_directory = directories["intermediate"] 

        # qiime dada2 denoise-paired \
            # --i-demultiplexed-seqs paired-end-demux.qza \
            # --p-trunc-len-f 251 \
            # --p-trunc-len-r 231 \
            # --o-table table.qza \
            # --o-representative-sequences rep-seqs.qza \
            # --o-denoising-stats denoising-stats.qza \


        # i/o
        input_filepaths = [
            os.path.join(directories["intermediate"], "forward_reads_trimming_suggestion", "position_to_trim.txt"),
            os.path.join(directories["intermediate"], "reverse_reads_trimming_suggestion", "position_to_trim.txt"),
            os.path.join(directories["intermediate"], "demultiplexed_reads.qza"),
        ]

        output_filenames = ["asv_table.qza", "asv_sequences.qza", "asv_denoising-stats.qza"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_dada2_cmd(**params)
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

        # -------------------------
        step = 4

        # Info
        program = "classify"
        program_label = "{}__{}".format(step, program)
        description = "Classify ASVs"

        # Add to directories
        output_directory = directories["intermediate"] 



        # i/o
        input_filepaths = [
            os.path.join(directories["intermediate"], "asv_sequences.qza"),
            opts.classifier,
        ]

        output_filenames = ["asv_taxonomy.qza"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_classify_cmd(**params)
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

        # -------------------------
        step = 5

        # Info
        program = "phylogeny"
        program_label = "{}__{}".format(step, program)
        description = "Phylogenetic inference on ASVs"

        # Add to directories
        output_directory = directories["intermediate"] 

        # i/o
        input_filepaths = [
            os.path.join(directories["intermediate"], "asv_sequences.qza"),
        ]

        output_filenames = ["asv_sequences.tree.unrooted.qza", "asv_sequences.tree.rooted.qza", "asv_sequences.aligned.qza", "asv_sequences.aligned.masked.qza"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))



        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_phylogeny_cmd(**params)
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


        # -------------------------
        step = 6

        # Info
        program = "conversion"
        program_label = "{}__{}".format(step, program)
        description = "Convert QIIME2 artifact objects"

        # Add to directories
        output_directory = directories["output"] 

        # i/o
        input_filenames = [ 
            # ASV
            "asv_table.qza", 
            "asv_sequences.qza", 
            "asv_denoising-stats.qza",
            "asv_taxonomy.qza",

            # Phylogeny
            "asv_sequences.tree.unrooted.qza", 
            "asv_sequences.tree.rooted.qza", 
            "asv_sequences.aligned.qza", 
            "asv_sequences.aligned.masked.qza",
        ]

        input_filepaths = list(map(lambda filename: os.path.join(directories["intermediate"], filename), input_filenames))


        output_filenames = [
            "feature-table.tsv",
            "taxonomy.tsv",
            "dna-sequences.with_taxonomy.fasta",
            "stats.tsv",
            "tree.nwk",
        ]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))



        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_conversion_cmd(**params)
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
    df = pd.read_csv(opts.reads_table, sep="\t")
    assert all(df.columns == pd.Index(["sample-id", "forward-absolute-filepath","reverse-absolute-filepath"])), "Header of --reads_table must be [sample-id]<tab>[forward-absolute-filepath]<tab>[reverse-absolute-filepath]"
    assert df.notnull().all(axis=None), "--reads_table cannot include missing values"
    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <reads_table.tsv> -c <classifier.qza> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i","--reads_table", type=str, required=True, help = "path/to/reads_table.tsv. 3 columns separated by tabs with the following header: [sample-id <tab> forward-absolute-filepath <tab> reverse-absolute-filepath]\nA utility script is provided: compile_reads_table.py")
    parser_io.add_argument("-c","--classifier", type=str, required=True, help = "path/to/feature_classifier. Data Resources: https://docs.qiime2.org/2022.8/data-resources/")
    parser_io.add_argument("-o","--project_directory", type=str, default="veba_output/amplicon", help = "path/to/project_directory [Default: veba_output/amplicon]")
    # parser_io.add_argument("--input_reads_format", type=str, default="PairedEndFastqManifestPhred33V2", help = "qiime tools import --input-format [Default: PairedEndFastqManifestPhred33V2]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    # parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
    parser_utility.add_argument("--tmpdir", type=str, help="Set temporary directory")  

    # Trim detection
    parser_trimming = parser.add_argument_group('Trim detection arguments')
    parser_trimming.add_argument("--inspect_trim_regions", action="store_true", help = "Manually inspect trim regions then rerun [PLEASE USE THIS TO CHECK TRIMMING SUGGESTIONS AS THEY ARE CURRENTLY EXPERIMENTAL]")
    parser_trimming.add_argument("-f","--forward_trim", type=int, help = "Specify forward trim position")
    parser_trimming.add_argument("-r","--reverse_trim", type=int, help = "Specify reverse trim position")
    parser_trimming.add_argument("-q","--minimum_quality",  default=30.0,type=float, help = "Minimum quality value")
    parser_trimming.add_argument("-m","--minimum_length",  default=100,type=int, help = "Minimum length.  If minimum quality value makes length shorter than this then an error will yield with which samples are responsible [Default: 100]")
    parser_trimming.add_argument("-w", "--window_size",  default=4,type=int, help = "Window size [Default: 4]")
    parser_trimming.add_argument("-l", "--maximum_average_loss",  type=float, help = "Maximum average loss for window size [Default: --window_size]")

   # DADA2
    parser_dada2 = parser.add_argument_group('DADA2 arguments')
    parser_dada2.add_argument("--minimum_overlap", type=int, default=12, help = "DADA2 | The minimum length of the overlap required for merging the forward and reverse reads. [Default: 12]")
    parser_dada2.add_argument("--dada2_options", type=str, default="", help = "Additional DADA2 options. '--arg value'")

   # Phylogeny
    parser_dada2 = parser.add_argument_group('Phylogeny arguments')
    parser_dada2.add_argument("--phylogeny_mode", type=str, default="align-to-tree-mafft-fasttree", help = "QIIME2 phylogeny submodule [Default: align-to-tree-mafft-fasttree]")
    parser_dada2.add_argument("--phylogeny_options", type=str, default="", help = "Additional options. '--arg value'")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    opts.path_config = "CONDA_PREFIX"
    opts.input_reads_format = "PairedEndFastqManifestPhred33V2"

    if opts.n_jobs == -1:
        from multiprocessing import cpu_count 
        opts.n_jobs = cpu_count()
    assert opts.n_jobs >= 1, "--n_jobs must be ≥ 1.  To select all available threads, use -1."

    if opts.project_directory.endswith("/"):
        opts.project_directory = opts.project_directory[:-1]

    # Directories
    directories = dict()
    directories["project"] = create_directory(opts.project_directory)
    # directories["sample"] = create_directory(os.path.join(directories["project"], opts.name))
    directories["output"] = create_directory(os.path.join(directories["project"], "output"))
    directories["log"] = create_directory(os.path.join(directories["project"], "log"))
    if not opts.tmpdir:
        opts.tmpdir = os.path.join(directories["project"], "tmp")
    directories["tmp"] = create_directory(opts.tmpdir)
    directories["checkpoints"] = create_directory(os.path.join(directories["project"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["project"], "intermediate"))
    os.environ["TMPDIR"] = directories["tmp"]

    # Info
    print(format_header("Warning: Automatic detection of forward and reverse trim regions is experimental.", "+"), file=sys.stderr)
    print(format_header("         Please use --inspect_trim_regions before proceeding with ASV detection to ensure that regions make sense.", "+"), file=sys.stderr)
    print(format_header("         This warning message will remain until thoroughly benchmarked.", "+"), file=sys.stderr)

    print(format_header(__program__, "="), file=sys.stdout)
    print(format_header("Configuration:", "-"), file=sys.stdout)
    # print(format_header("Name: {}".format(opts.name), "."), file=sys.stdout)
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
        pipeline = create_pipeline(
                     opts=opts,
                     directories=directories,
                     f_cmds=f_cmds,
        )
        pipeline.compile()
        pipeline.execute(restart_from_checkpoint=opts.restart_from_checkpoint)

if __name__ == "__main__":
    main()
