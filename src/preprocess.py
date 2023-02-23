#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
from collections import OrderedDict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *
import fastq_preprocessor


__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.2.23"

# ============
# Run Pipeline
# ============
# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = set([])

    required_executables={
                "repair.sh",
                "bbduk.sh",
                "bowtie2",
                "fastp",
                "seqkit",
                "fastq_preprocessor",
    } | accessory_scripts

    if opts.path_config == "CONDA_PREFIX":
        executables = dict()
        for name in required_executables:
            executables[name] = os.path.join(os.environ["CONDA_PREFIX"], "bin", name)
    else:
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

    assert opts.forward_reads != opts.reverse_reads, "You probably mislabeled the input files because `r1` should not be the same as `r2`: {}".format(opts.forward_reads)
    assert_acceptable_arguments(opts.retain_trimmed_reads, {0,1})
    assert_acceptable_arguments(opts.retain_decontaminated_reads, {0,1})

    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Wrapper around github.com/jolespin/fastq_preprocessor
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -1 <reads_1.fq> -2 <reads_2.fq> -n <name> -o <output_directory> |Optional| -x <reference_index> -k <kmer_database>".format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-1","--forward_reads", type=str, help = "path/to/reads_1.fastq")
    parser_io.add_argument("-2","--reverse_reads", type=str, help = "path/to/reads_2.fastq")
    parser_io.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    parser_io.add_argument("-o","--project_directory", type=str, default="veba_output/preprocess", help = "path/to/project_directory [Default: veba_output/preprocess]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv. Must have at least 2 columns [name, executable] [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=int, help = "Restart from a particular checkpoint")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Fastp
    parser_fastp = parser.add_argument_group('Fastp arguments')
    parser_fastp.add_argument("-m", "--minimum_read_length", type=int, default=75, help="Fastp | Minimum read length [Default: 75]")
    parser_fastp.add_argument("-a", "--adapters", type=str, default="detect", help="Fastp | path/to/adapters.fasta [Default: detect]")
    parser_fastp.add_argument("--fastp_options", type=str, default="", help="Fastp | More options (e.g. --arg 1 ) [Default: '']")

    # Bowtie
    parser_bowtie2 = parser.add_argument_group('Bowtie2 arguments')
    parser_bowtie2.add_argument("-x", "--contamination_index", type=str, help="Bowtie2 | path/to/contamination_index\n(e.g., Human GRCh38 from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids//GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz)")
    parser_bowtie2.add_argument("--retain_trimmed_reads", default=0, type=int, help = "Retain fastp trimmed fastq after decontamination. 0=No, 1=yes [Default: 0]") 
    parser_bowtie2.add_argument("--retain_contaminated_reads", default=0, type=int, help = "Retain contaminated fastq after decontamination. 0=No, 1=yes [Default: 0]")
    parser_bowtie2.add_argument("--bowtie2_options", type=str, default="", help="Bowtie2 | More options (e.g. --arg 1 ) [Default: '']\nhttp://bowtie-bio.sourceforge.net/bowtie2/manual.shtml")

    # BBDuk
    parser_bbduk = parser.add_argument_group('BBDuk arguments')
    parser_bbduk.add_argument("-k","--kmer_database", type=str,  help="BBDuk | path/to/kmer_database\n(e.g., ribokmers.fa.gz from https://drive.google.com/file/d/0B3llHR93L14wS2NqRXpXakhFaEk/view?usp=sharing)")
    parser_bbduk.add_argument("--kmer_size", type=int, default=31, help="BBDuk | k-mer size [Default: 31]")
    parser_bbduk.add_argument("--retain_kmer_hits", default=0, type=int, help = "Retain reads that map to k-mer database. 0=No, 1=yes [Default: 0]")
    parser_bbduk.add_argument("--retain_non_kmer_hits", default=0, type=int, help = "Retain reads that do not map to k-mer database. 0=No, 1=yes [Default: 0]")
    parser_bbduk.add_argument("--bbduk_options", type=str, default="", help="BBDuk | More options (e.g., --arg 1) [Default: '']")

    # Options
    opts = parser.parse_args()
    # opts.script_directory  = script_directory
    # opts.script_filename = script_filename

    # Threads
    if opts.n_jobs == -1:
        from multiprocessing import cpu_count 
        opts.n_jobs = cpu_count()
    assert opts.n_jobs >= 1, "--n_jobs must be ≥ 1.  To select all available threads, use -1."

    #Get arguments
    args = list() 
    for k,v in opts.__dict__.items():
        if v is not None:
            args += ["--{}".format(k), str(v)]
    # args = flatten(map(lambda item: ("--{}".format(item[0]), item[1]), opts.__dict__.items()))
    sys.argv = [sys.argv[0]] + args

    # Wrapper
    fastq_preprocessor.main(args)



if __name__ == "__main__":
    main()
