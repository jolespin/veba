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


# TransDecoder
def get_transdecoder_longorfs_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

# #  Transdecoder.LongOrfs|http://transdecoder.github.io> - Transcriptome Protein Prediction
# #
# #
# #  Required:
# #
# #    -t <string>                            transcripts.fasta
# #
# #  Optional:
# #
# #   --gene_trans_map <string>              gene-to-transcript identifier mapping file (tab-delimited, gene_id<tab>trans_id<return> )
# #
# #   -m <int>                               minimum protein length (default: 100)
# #
# #   -G <string>                            genetic code (default: universal; see PerlDoc; options: Euplotes, Tetrahymena, Candida, Acetabularia)
# #
# #   -S                                     strand-specific (only analyzes top strand)
# #
# #   --output_dir | -O  <string>            path to intended output directory (default:  basename( -t val ) + ".transdecoder_dir")
# #
# #   --genetic_code <string>                Universal (default)
# #
# #        Genetic Codes (derived from: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

    cmd = [
        os.environ["TransDecoder.LongOrfs"],
        "-t {}".format(opts.fasta),
        "-m {}".format(opts.minimum_protein_length),
        "-G {}".format(opts.genetic_code),
        "--output_dir {}".format(os.path.join(output_directory)),
        opts.transdecoder_longorfs_options,
        ]
    if opts.genes_to_transcripts:
        cmd += [ 
            "--gene_trans_map {}".format(opts.genes_to_transcripts),
        ]

    return cmd

def get_transdecoder_predict_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

#  Transdecoder.LongOrfs|http://transdecoder.github.io> - Transcriptome Protein Prediction
#
#
#  Required:
#
#   -t <string>                            transcripts.fasta
#
#  Common options:
#
#
#   --retain_long_orfs_mode <string>        'dynamic' or 'strict' (default: dynamic)
#                                        In dynamic mode, sets range according to 1%FDR in random sequence of same GC content.
#
#
#   --retain_long_orfs_length <int>         under 'strict' mode, retain all ORFs found that are equal or longer than these many nucleotides even if no other evidence
#                                         marks it as coding (default: 1000000) so essentially turned off by default.)
#
#   --retain_pfam_hits <string>            domain table output file from running hmmscan to search Pfam (see transdecoder.github.io for info)
#                                        Any ORF with a pfam domain hit will be retained in the final output.
#
#   --retain_blastp_hits <string>          blastp output in '-outfmt 6' format.
#                                        Any ORF with a blast match will be retained in the final output.
#
#   --single_best_only                     Retain only the single best orf per transcript (prioritized by homology then orf length)
#
#   --output_dir | -O  <string>            output directory from the TransDecoder.LongOrfs step (default: basename( -t val ) + ".transdecoder_dir")
#
#   -G <string>                            genetic code (default: universal; see PerlDoc; options: Euplotes, Tetrahymena, Candida, Acetabularia, ...)
#
#   --no_refine_starts                     start refinement identifies potential start codons for 5' partial ORFs using a PWM, process on by default.
#
##  Advanced options
#
#    -T <int>                            Top longest ORFs to train Markov Model (hexamer stats) (default: 500)
#                                        Note, 10x this value are first selected for removing redundancies,
#                                        and then this -T value of longest ORFs are selected from the non-redundant set.
#  Genetic Codes
#
#
#   --genetic_code <string>                Universal (default)
#
#        Genetic Codes (derived from: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

    cmd = [
        os.environ["TransDecoder.Predict"],
        "-t {}".format(opts.fasta),
        "--retain_long_orfs_mode {}".format(opts.mode),
        "--output_dir {}".format(output_directory),
        # "-G {}".format(oTrapts.genetic_code),
        opts.transdecoder_predict_options,
    ]
    if not opts.allow_multiple_orfs_per_transcript:
        cmd += [ 
            "--single_best_only",
        ]

    if opts.transdecoder_method != "default":
        if opts.transdecoder_method == "diamond,hmmer":
                cmd += [
                    "--retain_blastp_hits {}".format(input_filepaths[1]),
                    "--retain_pfam_hits {}".format(input_filepaths[2]),
                    ]
        else:
            if opts.transdecoder_method == "diamond":
                cmd += [
                    "--retain_blastp_hits {}".format(input_filepaths[1]),
                    ]

            if opts.transdecoder_method == "hmmer":
                cmd += [
                    "--retain_pfam_hits {}".format(input_filepaths[1]),
                    ]
    cmd += [ 
            "&&",
        "NAME=$(basename {})".format(opts.fasta),
            "&&",
        "cat $NAME.transdecoder.gff3 | {} > {}".format(
            os.environ["append_geneid_to_transdecoder_gff.py"],
            os.path.join(output_directory, "gene_models.gff"),
            ),
            "&&",
        "mv $NAME.transdecoder.gff3 {}".format(os.path.join(output_directory, "transdecoder.original.gff")),
            "&&",
        "mv $NAME.transdecoder.pep {}".format(os.path.join(output_directory, "gene_models.faa")),
            "&&",
        "mv $NAME.transdecoder.cds {}".format(os.path.join(output_directory, "gene_models.ffn")),
            "&&",
        "mv $NAME.transdecoder.bed {}".format(os.path.join(output_directory, "gene_models.bed")),

            "&&",

        "rm -f pipeliner.*.cmds",
    ]
    return cmd

# Diamond
def get_diamond_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [

        os.environ["diamond"],
        "blastp",
        "--db {}".format(opts.diamond_database),
        "--query {}".format(input_filepaths[0]),
        "--threads {}".format(opts.n_jobs),
        "-f 6",
        "--evalue {}".format(opts.diamond_evalue),
        "-o {}".format(os.path.join(output_directory, "output.tsv")),
        "--max-target-seqs 1",
        "--tmpdir {}".format(directories["tmp"]),
    ]
    if bool(opts.diamond_sensitivity):
        cmd += [ 
            "--{}".format(opts.diamond_sensitivity),
        ]
    cmd += [ 
        opts.diamond_options,
    ]
 

    return cmd

# HMMER
def get_hmmsearch_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [
        os.environ["hmmsearch"],
        "--tblout {}".format(os.path.join(output_directory, "output.tblout.tsv")),
        "--domtblout {}".format(os.path.join(output_directory, "output.domtblout.tsv")),
        "--cpu {}".format(opts.n_jobs),
        "--seed {}".format(opts.random_state + 1),
    ]
    if bool(opts.hmmsearch_threshold):
        cmd += [ "--{}".format(opts.hmmsearch_threshold)]
    else:
        cmd += ["-E {}".format(opts.hmmsearch_evalue)]


    cmd += [
        opts.hmm_database,
        input_filepaths[0],
        ">",
        "/dev/null",
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

    program = "transdecoder_longorfs"

    program_label = "{}__{}".format(step, program)

    # Add to directories


    # Info
    description = "Identifying long ORFs with TransDecoder"
    
    # i/o
    input_filepaths = [
        opts.fasta, 
        ]

    output_directory = directories[("intermediate",  "1__transdecoder")] = create_directory(os.path.join(directories["intermediate"], "1__transdecoder"))
    output_filenames = [ 
        "longest_orfs.*",
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))


    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_transdecoder_longorfs_cmd(**params)


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

    transdecoder_predict_input_filepaths = [opts.fasta]

    # ==========
    # Diamond
    # ==========
    if opts.diamond_database:
        step += 1

        program = "diamond"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

        # Info
        description = "Diamond"

        # i/o
        input_filepaths = [
            os.path.join(directories[("intermediate",  "1__transdecoder")], "longest_orfs.pep"),
            opts.diamond_database,
            ]
        output_filenames = ["output.tsv"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_diamond_cmd(**params)

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

        transdecoder_predict_input_filepaths += output_filepaths


    # ==========
    # HMMSearch
    # ==========
    if opts.hmm_database:

        step += 1

        program = "hmmsearch"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

        # Info
        description = "HMMSearch"

        # i/o
        input_filepaths = [
            os.path.join(directories[("intermediate",  "1__transdecoder")], "longest_orfs.pep"),
            opts.hmm_database,
            ]
        output_filenames = ["output.domtblout.tsv"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_hmmsearch_cmd(**params)

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

        transdecoder_predict_input_filepaths += output_filepaths

    # ==========
    # Predict
    # ==========
    step += 1

    program = "transdecoder_predict"

    program_label = "{}__{}".format(step, program)

    # Add to directories


    # Info
    description = "Predicting proteins with TransDecoder"
    
    # i/o
    input_filepaths = transdecoder_predict_input_filepaths

    output_directory = directories[("intermediate",  "1__transdecoder")] 
    output_filenames = [ 
        "gene_models.*",
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))


    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_transdecoder_predict_cmd(**params)


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
    accessory_scripts = set([
        "append_geneid_to_transdecoder_gff.py",
    ])

    required_executables={
        "TransDecoder.LongOrfs",
        "TransDecoder.Predict",
     } 

    if opts.diamond_database:
        required_executables.add("diamond")
    if opts.hmm_database:
        required_executables.add("hmmsearch")

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

    opts.transdecoder_method = "default"
    if any([opts.diamond_database, opts.hmm_database]):
        if opts.diamond_database and opts.hmm_database:
            opts.transdecoder_method = "diamond,hmmer"
        else:
            if opts.diamond_database and not opts.hmm_database:
                opts.transdecoder_method = "diamond"
            if opts.hmm_database and not opts.diamond_database:
                opts.transdecoder_method = "hmmer"



def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = """
____________________________________________________________________________________________________________________________________
*************
   WARNING
*************

Since TransDecoder writes output files to current working directory you 
can only run this for one sample at a time and not multiple concurrently.

Feature request submitted here: https://github.com/TransDecoder/TransDecoder/issues/169

Due to dependency conflicts, TransDecoder does not come in any VEBA environments.  Please create a new environment:

`conda create -n transdecoder_env -c bioconda -c jolespin genopype transdecoder diamond hmmer`

____________________________________________________________________________________________________________________________________

{} -f <transcripta.fasta> -o <output_directory> |Optional: -g <genes_to_transcripts.tsv> -H <database.hmm> -D <database.dmnd>

    """.format(__program__)

    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-f","--fasta", type=str, required=True, help = "path/to/transcripts.fasta")
    parser_io.add_argument("-o","--output_directory", type=str, default="transdecoder_output", help = "path/to/project_directory [Default: transdecoder_output]")
    parser_io.add_argument("-g","--genes_to_transcripts", type=str, required=False,  help = "path/to/genes_to_transcripts.tsv, [Optional] Format: [id_gene]<tab>[id_transcript], No header")
    parser_io.add_argument("-D","--diamond_database", type=str, required=False,  help = "path/to/diamond_database.dmnd")
    parser_io.add_argument("-H","--hmm_database", type=str, required=False,  help = "path/to/domains.hmm")


    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # TransDecoder
    parser_transdecoder = parser.add_argument_group('TransDecoder arguments')
    parser_transdecoder.add_argument("-m", "--minimum_protein_length", type=int, default=100, help="TransDecoder.LongOrfs | Minimum protein length  [Default: 100]")
    parser_transdecoder.add_argument("--mode", type=str, default="dynamic", help="TransDecoder.Predict | Mode for retaining long ORFs.  In dynamic mode, sets range according to 1 percent of FDR in random sequence of same GC content. {dynamic,strict} [Default: dynamic]")
    parser_transdecoder.add_argument("--genetic_code", type=str, default="universal", help="TransDecoder | Genetic code [Default: universal]")
    parser_transdecoder.add_argument("--allow_multiple_orfs_per_transcript", action="store_true", help="TransDecoder | Allow multiple ORFs per transcript [Default: Don't]")
    parser_transdecoder.add_argument("--transdecoder_longorfs_options", type=str, default="", help="TransDecoder.LongOrfs | More options (e.g. --arg 1 ) [Default: ''] https://github.com/TransDecoder/TransDecoder/wiki")
    parser_transdecoder.add_argument("--transdecoder_predict_options", type=str, default="", help="TransDecoder.Predict | More options (e.g. --arg 1 ) [Default: ''] https://github.com/TransDecoder/TransDecoder/wiki")

    # Diamond
    parser_diamond = parser.add_argument_group('Diamond arguments')
    parser_diamond.add_argument("--diamond_sensitivity", type=str, default="", help="Diamond | Sensitivity [Default:  '']")
    parser_diamond.add_argument("--diamond_evalue", type=float, default=0.001, help="Diamond | E-Value [Default: 0.001]")
    parser_diamond.add_argument("--diamond_options", type=str, default="", help="Diamond | More options (e.g. --arg 1 ) [Default: '']")

    # HMMER
    parser_hmmer = parser.add_argument_group('HMMSearch arguments')
    parser_hmmer.add_argument("--hmmsearch_threshold", type=str, default="cut_ga", help="HMMSearch | Threshold {cut_ga, cut_nc, gut_tc} [Default:  cut_ga]")
    parser_hmmer.add_argument("--hmmsearch_evalue", type=float, default=10.0, help="Diamond | E-Value [Default: 10.0]")
    parser_hmmer.add_argument("--hmmsearch_options", type=str, default="", help="Diamond | More options (e.g. --arg 1 ) [Default: '']")

    # Options
    opts = parser.parse_args()

    opts.script_directory  = script_directory
    opts.script_filename = script_filename

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

    
   


if __name__ == "__main__":
    main(sys.argv[1:])
