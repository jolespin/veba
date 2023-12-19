#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, gzip
from collections import OrderedDict, defaultdict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from genopype import __version__ as genopype_version
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.12.19"

# Preprocess reads
def get_sylph_sketch_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [
    os.environ["sylph"],
    "sketch",
    "-t {}".format(opts.n_jobs),
    "-c {}".format(opts.sylph_sketch_subsampling_rate),
    "-k {}".format(opts.sylph_sketch_k),
    "--min-spacing {}".format(opts.sylph_sketch_minimum_spacing),
    "-1 {}".format(opts.forward_reads),
    "-2 {}".format(opts.reverse_reads),
    "-d {}".format(output_directory),

        "&&",

    "mv",
    "-v",
    os.path.join(output_directory, "{}.paired.sylsp".format(os.path.split(opts.forward_reads)[1])),
    os.path.join(output_directory, "reads.sylsp"),
    ]

    return cmd

def get_sylph_profile_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    cmd = [
        os.environ["sylph"],
        "profile",
        "-t {}".format(opts.n_jobs),
        "--minimum-ani {}".format(opts.sylph_profile_minimum_ani),
        "--min-number-kmers {}".format(opts.sylph_profile_minimum_number_kmers),
        "--min-count-correct {}".format(opts.sylph_profile_minimum_count_correct),
        opts.sylph_profile_options,
        " ".join(input_filepaths),
        "|",
        "gzip",
        ">",
        os.path.join(output_directory, "sylph_profile.tsv.gz"),

            "&&",

        os.environ["reformat_sylph_profile_single_sample_output.py"],
        "-i {}".format(os.path.join(output_directory, "sylph_profile.tsv.gz")),
        "-o {}".format(output_directory),
        "-c {}".format(opts.genome_clusters) if opts.genome_clusters else "",
        "-f Taxonomic_abundance",
        "-x {}".format(opts.extension),
        "--header" if opts.header else "",
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
    accessory_scripts = set([
            "reformat_sylph_profile_single_sample_output.py",
    ]
    )

    required_executables={ 
                "sylph",
                # "seqkit",
                
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
    pipeline = ExecutablePipeline(name=__program__, description=opts.name, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])

    # ==========
    # Preprocess reads
    # ==========

    if opts.input_reads_format == "paired":

        step = 0

        # Info
        program = "sylph_sketch"
        program_label = "{}__{}".format(step, program)
        description = "Sketch input reads"
        
        # Add to directories
        output_directory = directories["output"]
        # i/o
        input_filepaths = [opts.forward_reads, opts.reverse_reads]
        output_filepaths = [
            os.path.join(output_directory, "reads.sylsp"),
            ]

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_sylph_sketch_cmd(**params)
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
    else:
        output_filepaths = [opts.reads_sketch]
    

    # ==========
    # Profile
    # ==========
    
    step = 1

    # Info
    program = "sylph_profile"
    program_label = "{}__{}".format(step, program)
    description = "Profile genome databases"
    
    # Add to directories
    output_directory = directories["output"] 

    # i/o
    input_filepaths = output_filepaths + opts.sylph_databases


    output_filepaths = [
            os.path.join(output_directory,  "sylph_profile.tsv.gz"),
            os.path.join(output_directory,  "taxonomic_abundance.tsv.gz"),
        ]
    if opts.genome_clusters:
        input_filepaths += [
            opts.genome_clusters,
        ]
        output_filepaths += [
            os.path.join(output_directory,  "taxonomic_abundance.clusters.tsv.gz"),
        ]

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_sylph_profile_cmd(**params)
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

    for db in opts.sylph_databases:
        assert db.endswith(".syldb"), "{} must have .syldb file extension".format(db)

 # --input_reads_format
    assert_acceptable_arguments(opts.input_reads_format, {"paired",  "sketch", "auto"})
    if opts.input_reads_format == "auto":
        if any([opts.forward_reads, opts.reverse_reads]):
            assert opts.forward_reads != opts.reverse_reads, "You probably mislabeled the input files because `forward_reads` should not be the same as `reverse_reads`: {}".format(opts.forward_reads)
            assert opts.forward_reads is not None, "If running in --input_reads_format paired mode, --forward_reads and --reverse_reads are needed."
            assert opts.reverse_reads is not None, "If running in --input_reads_format paired mode, --forward_reads and --reverse_reads are needed."
            opts.input_reads_format = "paired"
        if opts.reads_sketch is not None:
            assert opts.forward_reads is None, "If running in --input_reads_format sketch mode, you cannot provide --forward_reads, --reverse_reads"
            assert opts.reverse_reads is None, "If running in --input_reads_format sketch mode, you cannot provide --forward_reads, --reverse_reads"
            opts.input_reads_format = "sketch"

        print("Auto detecting reads format: {}".format(opts.input_reads_format), file=sys.stdout)
    assert_acceptable_arguments(opts.input_reads_format, {"paired", "sketch"})

    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -1 <forward_reads.fq> -2 <reverse_reads.fq>|-s <sketch> -n <name> -o <output_directory> -d <db_1.syldb db_2.syldb ... db_n.syldb>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    
    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-1","--forward_reads", type=str,  help = "path/to/forward_reads.fq[.gz]")
    parser_io.add_argument("-2","--reverse_reads", type=str,  help = "path/to/reverse_reads.fq[.gz]]")
    parser_io.add_argument("-s","--reads_sketch", type=str, help = "path/to/reads_sketch.sylsp (e.g., sylph sketch output) (Cannot be used with --forward_reads and --reverse_reads)")
    parser_io.add_argument("-n", "--name", type=str, required=True, help="Name of sample")
    parser_io.add_argument("-d","--sylph_databases", type=str, nargs="+", required=True, help = "Sylph database(s) with all genomes.  Can be multiple databases delimited by spaces.  Use compile_custom_sylph_sketch_database_from_genomes.py to build database.") 
    parser_io.add_argument("-o","--project_directory", type=str, default="veba_output/profiling/taxonomy", help = "path/to/project_directory [Default: veba_output/profiling/taxonomy]")
    parser_io.add_argument("-c","--genome_clusters", type=str, help = "path/to/mags_to_slcs.tsv. [id_genome]<tab>[id_genome-cluster], No header. Aggregates counts for genome clusters.")
    parser_io.add_argument("-F", "--input_reads_format", choices={"paired", "sketch"}, type=str, default="auto", help = "Input reads format {paired, sketch} [Default: auto]")
    parser_io.add_argument("-x","--extension", type=str, default="fa", help = "Fasta file extension for bins. Assumes all genomes have the same file extension. [Default: fa]")


    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
    parser_utility.add_argument("--tmpdir", type=str, help="Set temporary directory")  #site-packges in future

    # Sylph
    parser_sylph_sketch = parser.add_argument_group('Sylph sketch arguments (Fastq)')
    parser_sylph_sketch.add_argument("--sylph_sketch_k", type=int, choices={21,31}, default=31,  help="Sylph sketch [Fastq] |  Value of k. Only k = 21, 31 are currently supported. [Default: 31]")
    parser_sylph_sketch.add_argument("--sylph_sketch_minimum_spacing", type=int,  default=30,  help="Sylph sketch [Fastq] |  Minimum spacing between selected k-mers on the genomes [Default: 30]")
    parser_sylph_sketch.add_argument("--sylph_sketch_subsampling_rate", type=int, default=100,  help="Sylph sketch [Fastq] |  Subsampling rate.	 sylph runs without issues if the -c for all genomes is ≥ the -c for reads.  [Default: 100]")
    parser_sylph_sketch.add_argument("--sylph_sketch_options", type=str, default="", help="Sylph sketch [Fastq] | More options for `sylph sketch` (e.g. --arg 1 ) [Default: '']")

    parser_sylph_profile = parser.add_argument_group('Sylph profile arguments')
    parser_sylph_profile.add_argument("--sylph_profile_minimum_ani", type=float, default=95, help="Sylph profile | Minimum adjusted ANI to consider (0-100). [Default: 95]")
    parser_sylph_profile.add_argument("--sylph_profile_minimum_number_kmers", type=int, default=20, help="Sylph profile | Exclude genomes with less than this number of sampled k-mers.  Default is 50 in Sylph but lowering to 20 accounts for viruses and small CPR genomes. [Default: 20]")
    parser_sylph_profile.add_argument("--sylph_profile_minimum_count_correct", type=int, default=3, help="Sylph profile | Minimum k-mer multiplicity needed for coverage correction. Higher values gives more precision but lower sensitivity [Default: 3]")
    parser_sylph_profile.add_argument("--sylph_profile_options", type=str, default="", help="Sylph profile | More options for `sylph profile` (e.g. --arg 1 ) [Default: '']")
    parser_sylph_profile.add_argument("--header", action="store_true",  help = "Include header in taxonomic abundance tables")

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
    print("GenoPype version:", genopype_version, file=sys.stdout) #sys.path[2]
    print("Script version:", __version__, file=sys.stdout)
    print("Moment:", get_timestamp(), file=sys.stdout)
    print("Directory:", os.getcwd(), file=sys.stdout)
    if "TMPDIR" in os.environ: print(os.environ["TMPDIR"], file=sys.stdout)
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
