#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, random
from collections import OrderedDict, defaultdict

import pandas as pd
import numpy as np

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.04.11"



def get_maxbin2_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Create dummy scaffolds_to_bins.tsv to overwrite later. This makes DAS_Tool easier to run
    cmd = [
         "echo '' > {}".format(os.path.join(output_directory, "scaffolds_to_bins.tsv")),
         "&&",
         "mkdir -p {}".format(os.path.join(output_directory, "bins")),

    ]

    if opts.bam:
        # Coverage for MaxBin2 
        cmd += [
            "&&",
            os.environ["coverm"],
            "contig",
            "--threads {}".format(opts.n_jobs),
            "--methods metabat",
            "--bam-files",
            " ".join(opts.bam),
            "|",
            # Coverage for MaxBin2
            "cut -f1,4",
            "|",
            "tail -n +2",
            ">",
            os.path.join(output_directory,"intermediate", "coverage.tsv"),

            "&&",

            os.environ["run_MaxBin.pl"],
            "-contig {}".format(opts.fasta), # scaffolds.fasta
            "-out {}".format(os.path.join(output_directory, "intermediate", "bin")),
            "-abund {}".format(os.path.join(output_directory, "intermediate", "coverage.tsv")),
            "-min_contig_length {}".format(opts.minimum_contig_length),
            "-markerset {}".format(opts.maxbin2_markerset),
            "-thread {}".format(opts.n_jobs),
            "-verbose",
            opts.maxbin2_options,
        ]
    else:
        cmd += [ 
            "&&",

            os.environ["run_MaxBin.pl"],
            "-contig {}".format(opts.fasta), # scaffolds.fasta
            "-out {}".format(os.path.join(output_directory, "intermediate", "bin")),
            "-abund {}".format(opts.coverage),
            "-min_contig_length {}".format(opts.minimum_contig_length),
            "-markerset {}".format(opts.maxbin2_markerset),
            "-thread {}".format(opts.n_jobs),
            "-verbose",
            opts.maxbin2_options,
        ]

    # Move bins and change extension
    cmd += [ 
"""

for FP in %s;
    do ID_GENOME=$(basename ${FP} .fasta);
    mv $FP %s/${ID_GENOME}.fa
    done

"""%( 
    os.path.join(output_directory, "intermediate", "*.fasta"),
    os.path.join(output_directory, "bins"),
    )
    ]

    # Remove small bins
    cmd += [ 
r"""
for FP in %s;
    do GENOME_SIZE=$(cat $FP | grep -v "^>" | tr -d "\n" | wc -m);
    echo "[GENOME SIZE] ${FP} = ${GENOME_SIZE}"

    if (( %d > ${GENOME_SIZE} )); then
        echo "[REMOVING] ${FP}"
        rm -f ${FP}
    fi

    done


"""%( 
    os.path.join(output_directory, "bins", "*.fa"),
    opts.minimum_genome_length,
    )
    
    ]


    return cmd

# Metabat2
def get_metabat2_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Create dummy scaffolds_to_bins.tsv to overwrite later. This makes DAS_Tool easier to run
    cmd = [
         "echo '' > {}".format(os.path.join(output_directory, "scaffolds_to_bins.tsv")),
         "&&",
         "mkdir -p {}".format(os.path.join(output_directory, "bins")),

    ]

    # Coverage for MetaBat2 
    if opts.bam:
        cmd += [
        "&&",
        os.environ["coverm"],
        "contig",
        "--threads {}".format(opts.n_jobs),
         "--methods metabat",
         "--bam-files",
        " ".join(opts.bam),
        ">",
        os.path.join(output_directory, "intermediate", "coverage.tsv"),
        
        "&&",
       
        os.environ["metabat2"],
        "-i {}".format(opts.fasta), # scaffolds.fasta
        "-o {}".format(os.path.join(output_directory, "bins", "bin")),
        "-a {}".format(os.path.join(output_directory, "intermediate", "coverage.tsv")),
        "-m {}".format(max(opts.minimum_contig_length, 1500)), # mininimum contig length must be >= 1500 nts for MetaBat2
        "-t {}".format(opts.n_jobs),
        "--minClsSize {}".format(opts.minimum_genome_length),
        # "--unbinned", # Don't forget about these...
        "--seed {}".format(opts.random_state),
        "--verbose",
        opts.metabat2_options,

        ]
    # Coeverage provided
    else:
        cmd += [
         "&&",
        os.environ["metabat2"],
        "-i {}".format(opts.fasta), # scaffolds.fasta
        "-o {}".format(os.path.join(output_directory, "bins", "bin")),
        "-a {}".format(opts.coverage),
        "-m {}".format(max(opts.minimum_contig_length, 1500)), # mininimum contig length must be >= 1500 nts for MetaBat2
        "-t {}".format(opts.n_jobs),
        "--minClsSize {}".format(opts.minimum_genome_length),
        # "--unbinned", # Don't forget about these...
        "--seed {}".format(opts.random_state),
        "--verbose",
        opts.metabat2_options,
        ]

#     # Remove small bins
#     cmd += [ 
# """

# for FP in %s;
#     do GENOME_SIZE=$(cat $FP | grep -v "^>" | tr -d "\n" | wc -m);
#     echo "[GENOME SIZE] ${FP} = ${GENOME_SIZE}"

#     if (( %d > ${GENOME_SIZE} )); then
#         echo "[REMOVING] ${FP}"
#         rm -f ${FP}
#     fi

#     done


# """%( 
#     os.path.join(output_directory, "bins", "*.fa"),
#     opts.minimum_genome_size,
#     )
#     ]
    
    return cmd


def get_concoct_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [
        # Create dummy scaffolds_to_bins.tsv to overwrite later. This makes DAS_Tool easier to run
        "echo '' > {}".format(os.path.join(output_directory, "scaffolds_to_bins.tsv")),
        "&&",
        "mkdir -p {}".format(os.path.join(output_directory, "bins")),
        "&&",
        # cut_up_fasta.py original_contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
        os.environ["cut_up_fasta.py"],
        opts.fasta, 
        "-c {}".format(opts.concoct_fragment_length),
        "-o {}".format(opts.concoct_overlap_length),
        "--merge_last",
        "-b {}".format(os.path.join(output_directory, "intermediate", "scaffolds_fragmented.bed")),
        ">",
        os.path.join(output_directory, "intermediate", "scaffolds_fragmented.fasta"),

        "&&",

        # concoct_coverage_table.py contigs_10K.bed mapping/Sample*.sorted.bam > coverage_table.tsv
        os.environ["concoct_coverage_table.py"],
        os.path.join(output_directory, "intermediate", "scaffolds_fragmented.bed"), 
        " ".join(opts.bam), 
        ">",
        os.path.join(output_directory, "intermediate", "coverage_table.tsv"), 

        "&&",

        # concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output/
        os.environ["concoct"],
        "--composition_file {}".format(os.path.join(output_directory, "intermediate", "scaffolds_fragmented.fasta")),
        "--coverage_file {}".format(os.path.join(output_directory, "intermediate", "coverage_table.tsv")),
        "--length_threshold {}".format(opts.minimum_contig_length),
        "--seed {}".format(opts.random_state),
        "--threads {}".format(opts.n_jobs),
        "-b {}".format(os.path.join(output_directory, "intermediate")),
        opts.concoct_options,

        "&&",
        # merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv

        os.environ["merge_cutup_clustering.py"],
        os.path.join(output_directory, "intermediate", "clustering_gt{}.csv".format(opts.minimum_contig_length)),
        ">",
        os.path.join(output_directory, "intermediate", "clustering_gt{}.merged.csv".format(opts.minimum_contig_length)),


        # extract_fasta_bins.py original_contigs.fa concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins

        "&&",

        os.environ["extract_fasta_bins.py"],
        opts.fasta,
        os.path.join(output_directory, "intermediate", "clustering_gt{}.merged.csv".format(opts.minimum_contig_length)),
        "--output_path {}".format(os.path.join(output_directory, "bins/")),

    ]

        # Remove small bins
    cmd += [ 
r"""

for FP in %s;
    do GENOME_SIZE=$(cat $FP | grep -v "^>" | tr -d "\n" | wc -m);
    echo "[GENOME SIZE] ${FP} = ${GENOME_SIZE}"

    if (( %d > ${GENOME_SIZE} )); then
        echo "[REMOVING] ${FP}"
        rm -f ${FP}
    fi

    done


"""%( 
    os.path.join(output_directory, "bins", "*.fa"),
    opts.minimum_genome_length,
    )
    ]

    return cmd



def get_compile_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # scaffolds_to_bins
    cmd = [
        os.environ["scaffolds_to_bins.py"],
        "-x fa",
        "-i {}".format(os.path.join(output_directory, "bins")),
    ]
    if opts.bin_prefix:
        cmd += ["--bin_prefix {}".format(opts.bin_prefix)]

    cmd += [
        ">",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),

        # "binned.list"
       "&&",

        "cut",
        "-f1",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),
        ">",
        os.path.join(output_directory, "binned.list"),

        # bins.list
        "&&",

        "cut -f2",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),
        "|",
        "sort",
        "-u",
        ">",
        os.path.join(output_directory, "bins.list"),

        # unbinned.list
        "&&",
        "cat",
        opts.fasta,
        "|",
        "grep",
        '"^>"',
        "|",
        "cut -c2-",
        ">",
        os.path.join(output_directory, "unbinned.list"),
    ]

    # unbinned.fasta
    if opts.retain_unbinned_fasta:
        cmd += [ 
        "&&",
        "cat",
        opts.fasta,
        "|",
        os.environ["seqkit"],
        "grep",
        "--pattern-file {}".format(os.path.join(output_directory, "binned.list")),
        "-v",
        ">",
        os.path.join(output_directory, "unbinned.fasta"),
        ]

    # relabel with prefix
    if opts.bin_prefix:
        cmd += [ 
    """

    for FP in %s;
        do ID_GENOME=$(basename ${FP} .fa);
        mv $FP %s/%s${ID_GENOME}.fa
        done

    """%( 
        os.path.join(output_directory, "bins", "*.fa"),
        os.path.join(output_directory, "bins"),
        opts.bin_prefix,
        )
        ]
    else:
        cmd += ["&&"]

    cmd += [ 
        # "genome_statistics.tsv"
        # "&&",
        os.environ["seqkit"],
        "stats",
        "--basename",
        "--all",
        "-T",
        "-j {}".format(opts.n_jobs),
        os.path.join(output_directory, "bins", "*.fa"),
        ">",
        os.path.join(output_directory,  "genome_statistics.tsv"),
    ]


    if opts.remove_bins:
        cmd += [ 
            "&&",
            "rm -rf {}".format(os.path.join(output_directory,"bins")),
        ]

    # if opts.remove_intermediate_files:
    #     cmd += [ 
    #         "&&",
    #         "rm -rf {}".format(os.path.join(output_directory,"intermediate")),
    #     ]

    return cmd


# ============
# Run Pipeline
# ============



def create_pipeline(opts, directories, f_cmds):

    # .................................................................
    # Primordial
    # .................................................................
    # Commands file
    pipeline = ExecutablePipeline(name=__program__, description=opts.algorithm, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])


    if opts.random_state == -1:
        opts.random_state = random.choice(range(1,10000))


    # ==========
    # Binning
    # ==========
    step = 1

    program = "binning_{}".format(opts.algorithm)

    program_label = "{}__{}".format(step, program)

    # Add to directories


    # Info
    description = "Binning via {}".format(opts.algorithm)
    
    # i/o
    output_directory = directories["output"]
    output_filepaths = [
        os.path.join(output_directory, "bins", "*.fa"),
    ]

    if opts.algorithm == "concoct":
        # Output directory


        input_filepaths = [opts.fasta, opts.bam]

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_concoct_cmd(**params)


    if opts.algorithm == "metabat2":
        # Output directory


        if opts.bam:
            input_filepaths = [
                opts.fasta, 
                *opts.bam,
                ]
        else:
            input_filepaths = [
                opts.fasta, 
                opts.coverage,
                ]            


        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_metabat2_cmd(**params)

    if opts.algorithm == "maxbin2":
        # Output directory
        output_filepaths = [
            os.path.join(output_directory, "bins", "*.fa"),
        ]
        if opts.bam:
            input_filepaths = [
                opts.fasta, 
                *opts.bam,
                ]
        else:
            input_filepaths = [
                opts.fasta, 
                opts.coverage,
                ]            


        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_maxbin2_cmd(**params)


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


   # =============
    # Compile
    # =============
    step = 2
    program = "compile"

    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories["output"]

    # Info

    description = "Compiling results for output"

    # i/o
    input_filepaths = [
        os.path.join(directories["output"],  "bins","*.fa"),
    ]


    output_filenames =  [
        "scaffolds_to_bins.tsv", 
        "binned.list",
        "unbinned.list",
        "bins.list",
        "genome_statistics.tsv",
    ]
    if not opts.remove_bins:
        output_filenames += ["bins/*.fa"]
    if opts.retain_unbinned_fasta:
        output_filenames += ["unbinned.fasta"]

    output_filepaths = list(map(lambda fn:os.path.join(directories["output"], fn), output_filenames))

    
    params = {
    "input_filepaths":input_filepaths,
    "output_filepaths":output_filepaths,
    "output_directory":output_directory,
    "opts":opts,
    "directories":directories,
    }

    cmd = get_compile_cmd(**params)
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
    accessory_scripts = {
        "scaffolds_to_bins.py",
    }

    required_executables={
                "seqkit",
                # "tiara",
     } 
    if opts.algorithm in {"maxbin2", "metabat2"}:
        if opts.bam:
            required_executables |= {"coverm"}
        if opts.algorithm == "maxbin2":
                required_executables |= {
                    "run_MaxBin.pl",
                    }
        if opts.algorithm == "metabat2":
                required_executables |= {
                    "metabat2",
                    }

    if opts.algorithm == "concoct":
         required_executables  |= {
                "concoct",
                "cut_up_fasta.py", 
                "concoct_coverage_table.py",
                "merge_cutup_clustering.py", 
                "extract_fasta_bins.py", 
            } 
    # if opts.algorithm == "metacoag":
    #     required_executables |= {"MetaCoAG"}
        

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
        if name.endswith(".py"):
            executables[name] = "python " + os.path.join(opts.script_directory, name)
        else: 
            executables[name] = os.path.join(opts.script_directory, name)


    print(format_header( "Adding executables to path from the following source: {}".format(opts.path_config), "-"), file=sys.stdout)
    for name, executable in executables.items():
        if name in required_executables:
            print(name, executable, sep = " --> ", file=sys.stdout)
            os.environ[name] = executable.strip()
    print("", file=sys.stdout)


# Configure parameters
def configure_parameters(opts, directories):
    # assert opts.reference_assembly is not None, "Must include --reference_assembly"
    assert not all([bool(opts.bam), bool(opts.coverage)]), "Cannot have both --bam and --coverage"
    if opts.algorithm in {"metabat2", "maxbin2"}:
        assert  any([bool(opts.bam), bool(opts.coverage)]), "Must have either --bam or --coverage for --algorithm in metabat2 or maxbin2"
    if opts.algorithm == "concoct":
        assert opts.bam is not None, "Must provide --bam if --algorithm = concoct"

    if opts.bin_prefix == "DEFAULT":
        if opts.algorithm == "maxbin2":
            opts.bin_prefix = "{}-{}__".format(opts.algorithm.upper(), opts.maxbin2_markerset)

        else:
            opts.bin_prefix = "{}__".format(opts.algorithm.upper())

    if opts.bin_prefix == "NONE":
        opts.bin_prefix = None
    




    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(argv=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -f <scaffolds.fasta> -b <mapped.sorted.bam>|-c <coverage.tsv> -a <algorithm> -o <output_directory>".format(__program__)

    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-f","--fasta", type=str, required=True, help = "path/to/scaffolds.fasta")
    parser_io.add_argument("-b","--bam", type=str, required=False, nargs="+", help = "path/to/mapped.sorted.bam files separated by spaces. Can be used with any binning algorithm but cannot be used with --coverage")
    parser_io.add_argument("-c","--coverage", type=str, required=False, help = "path/to/coverage.tsv. Can be used with --algorithm metabat2 maxbin2 but not available for --algorithm concoct and cannot be used with --bam")
    parser_io.add_argument("-o","--output_directory", type=str,  help = "path/to/binning_output [Default: binning_output/[algorithm_name]]")
    # parser_io.add_argument("-M","--multisplit", type=str,  help = "path/to/scaffolds_to_samples.tsv. Use this to perform multisplit binning. Expected table input is the following format: [id_scaffold]<tab>[id_sample], No header. [Optional]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')

    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Binning
    parser_binning = parser.add_argument_group('Binning arguments')
    parser_binning.add_argument("-a", "--algorithm", type=str, default="metabat2", help="Binning algorithm: {concoct, metabat2, maxbin2} Future: {metacoag, vamb} [Default: metabat2] ")
    parser_binning.add_argument("-m", "--minimum_contig_length", type=int, default=1500, help="Minimum contig length.  [Default: 1500] ")
    parser_binning.add_argument("-s", "--minimum_genome_length", type=int, default=150000, help="Minimum genome length.  [Default: 150000] ")
    parser_binning.add_argument("-P","--bin_prefix", type=str,  default="DEFAULT", help = "Prefix for bin names. Special strings include: 1) --bin_prefix NONE which does not include a bin prefix; and 2) --bin_prefix DEFAULT then prefix is [ALGORITHM_UPPERCASE]__")
    parser_binning.add_argument("-B", "--remove_bins", action="store_true", help="Remove fasta files for bins")
    parser_binning.add_argument("-U", "--retain_unbinned_fasta", action="store_true", help="Retain unbinned fasta")
    # parser_binning.add_argument("-R", "--remove_intermediate_files", action="store_true", help="Remove intermediate files")


    parser_maxbin2 = parser.add_argument_group('MaxBin2 arguments')
    parser_maxbin2.add_argument("--maxbin2_markerset", type=int, default=107, help="MaxBin2 marker set: {107, 40} [Default: 107] ")
    parser_maxbin2.add_argument("--maxbin2_options", type=str, default="", help="MaxBin2 | More options (e.g. --arg 1 ) [Default: ''] | https://bitbucket.org/berkeleylab/metabat/src/master/")


    parser_metabat2 = parser.add_argument_group('Metabat2 arguments')
    parser_metabat2.add_argument("--metabat2_options", type=str, default="", help="MetaBat2 | More options (e.g. --arg 1 ) [Default: ''] | https://bitbucket.org/berkeleylab/metabat/src/master/")

    parser_conccoct = parser.add_argument_group('CONCOCT arguments')
    parser_conccoct.add_argument("--concoct_fragment_length", type=int, default=10000, help="CONCOCT | Fragment length [Default: 10000] ")
    parser_conccoct.add_argument("--concoct_overlap_length", type=int, default=0, help="CONCOCT | Fragment overlap length [Default: 0] ")
    parser_conccoct.add_argument("--concoct_options", type=str, default="", help="CONCOCT | More options (e.g. --arg 1 ) [Default: '']")

    # parser_metacoag = parser.add_argument_group('MetaCoAG arguments')
    # parser_metacoag.add_argument("--metacoag_options", type=str, default="", help="MetaCoAG | More options (e.g. --arg 1 ) [Default: '']")


    # # Tiara
    # parser_domain = parser.add_argument_group('Domain classification arguments')
    # parser_domain.add_argument("--logit_transform", type=str, default="softmax", help = " Transformation for consensus_domain_classification: {softmax, tss} [Default: softmax")
    # parser_domain.add_argument("--tiara_minimum_length", type=int, default=3000, help="Tiara | Minimum contig length. Anything lower than 3000 is not recommended. [Default: 3000]")
    # parser_domain.add_argument("--tiara_options", type=str, default="", help="Tiara | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/ibe-uw/tiara")


    # Options
    opts = parser.parse_args(argv)

    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    

    # Directories

    assert_acceptable_arguments(opts.algorithm, {"concoct", "metabat2", "maxbin2"})
    # if opts.multisplit:
    #     assert opts.algorithm != "maxbin2", "MaxBin2 is marker-based and cannot be used for multisplit binning"



    if not opts.output_directory:
        if opts.algorithm == "maxbin2":
            opts.output_directory = "binning_output/{}-{}".format(opts.algorithm, opts.maxbin2_markerset)
        else:
            opts.output_directory = "binning_output/{}".format(opts.algorithm)

    directories = dict()
    directories["output"] = create_directory(opts.output_directory)
    directories["intermediate"] = create_directory(os.path.join(directories["output"], "intermediate"))
    directories["log"] = create_directory(os.path.join(directories["output"],"intermediate", "log"))
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
