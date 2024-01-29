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
__version__ = "2023.11.30"

DIAMOND_DATABASE_SUFFIX = "_v201901b.dmnd"

# Preprocess reads
def get_preprocess_reads_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    if opts.input_reads_format == "paired":
        cmd = [
            os.environ["repair.sh"],
            "in1={}".format(opts.forward_reads),
            "in2={}".format(opts.reverse_reads),
            "out=stdout.fastq",
            "|",
            os.environ["bbmerge.sh"],
            "minoverlap={}".format(opts.minimum_merge_overlap),
            "int=t",
            "in=stdin.fastq",
            "out={}".format(os.path.join(output_directory, "joined.fastq.gz")),
            opts.bbmerge_options,

                "&&",

            os.environ["seqkit"],
            "stats",
            "-T",
            "-j {}".format(opts.n_jobs),
            opts.forward_reads,
            opts.reverse_reads,
            os.path.join(output_directory, "joined.fastq.gz"),
            ">",
            os.path.join(output_directory, "reads.seqkit_stats.tsv"),

        ]

    if opts.input_reads_format == "bam":
        cmd = [
            os.environ["samtools"],
            "fastq",
            "--threads {}".format(opts.n_jobs),
            opts.bam,
            "|",
            os.environ["repair.sh"],
            "in=stdin.fastq",
            "int=t",
            "out=stdout.fastq",
            "|",
            os.environ["bbmerge.sh"],
            "minoverlap={}".format(opts.minimum_merge_overlap),
            "int=t",
            "in=stdin.fastq",
            "out={}".format(os.path.join(output_directory, "joined.fastq.gz")),
            opts.bbmerge_options,

                "&&",

            os.environ["seqkit"],
            "stats",
            "-T",
            "-j {}".format(opts.n_jobs),
            os.path.join(output_directory, "joined.fastq.gz"),
            ">",
            os.path.join(output_directory, "reads.seqkit_stats.tsv"),
        ]

    if opts.input_reads_format == "joined":
        cmd = [
            "echo",
            '"[Skipping] Reads are already joined: {}"'.format(opts.joined_reads),

                "&&",

            os.environ["seqkit"],
            "stats",
            "-T",
            "-j {}".format(opts.n_jobs),
            opts.joined_reads,
            ">",
            os.path.join(output_directory, "reads.seqkit_stats.tsv"),
        ]
   
    return cmd

def get_preprocess_database_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command

    if opts.database_format == "fasta":
        cmd = [
            os.environ["diamond"],
            "makedb",
            "--threads {}".format(opts.n_jobs),
            "--in {}".format(opts.fasta),
            "--db {}".format(os.path.join(output_directory, "database{}".format(DIAMOND_DATABASE_SUFFIX)))
        ]

    if opts.database_format == "diamond_database":
        cmd = [
            "DST={}; SRC={}; SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST/{}".format(
                output_directory,
                opts.diamond_database,
                "database{}".format(DIAMOND_DATABASE_SUFFIX),
            )
        ]
    
    cmd += [ 

            "&&",

        os.environ["diamond"],
        "dbinfo",
        "--threads {}".format(opts.n_jobs),
        "--db {}".format(os.path.join(output_directory,  "database{}".format(DIAMOND_DATABASE_SUFFIX))),
        ">",
        os.path.join(output_directory, "database{}.dbinfo.txt".format(DIAMOND_DATABASE_SUFFIX)),
    ]
   
    return cmd

def get_humann_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    cmd = [
        os.environ["humann"],
        "-i {}".format(input_filepaths[0]),
        "-o {}".format(output_directory),
        "--threads {}".format(opts.n_jobs),
        "-v",
        "--search-mode {}".format(opts.search_mode),
        "--memory-use {}".format(opts.humann_memory),
        "--input-format {}".format("fastq.gz" if input_filepaths[0].endswith(".gz") else "fastq"),
        "--bypass-nucleotide-search",
        "--diamond {}".format(os.path.split(os.environ["diamond"])[0]),
        "--evalue {}".format(opts.evalue),
        "--protein-database {}".format(directories[("intermediate", "2__diamond_database")]),
        "--translated-identity-threshold {}".format(opts.translated_identity_threshold) if opts.translated_identity_threshold else "",
        "--translated-query-coverage-threshold {}".format(opts.translated_query_coverage_threshold),
        "--translated-subject-coverage-threshold {}".format(opts.translated_subject_coverage_threshold),
        "--output-basename humann",
        "--pathways {}".format(opts.pathways),
        "--id-mapping {}".format(opts.identifier_mapping),
        "--o-log {}".format(os.path.join(output_directory, "humann.log")),
        "--remove-column-description-output",
        opts.humann_options,

            "&&",

        # humann_diamond_unaligned.fa
        "mv",
        os.path.join(output_directory, "humann_humann_temp", "humann_diamond_unaligned.fa"),
        output_directory,

            "&&",

        # humann_diamond_aligned.tsv
        "echo", # Write headers
        "-e",
        '"qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"',
        ">",
        os.path.join(output_directory, "humann_diamond_aligned.tsv"),

            "&&",

        "cat",
        os.path.join(output_directory, "humann_humann_temp", "humann_diamond_aligned.tsv"),
        ">>",
        os.path.join(output_directory, "humann_diamond_aligned.tsv"), #blast6 format

            "&&",

        "rm -rf",
        os.path.join(output_directory, "humann_humann_temp"),

            "&&",

        "gzip",
        os.path.join(output_directory, "humann.log"),
        os.path.join(output_directory, "humann_diamond_aligned.tsv"),
        os.path.join(output_directory, "humann_diamond_unaligned.fa"),
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
    accessory_scripts = set([
                ]
    )

    required_executables={ 
                "bbmerge.sh",
                "repair.sh",
                "samtools",
                "humann",
                # "humann_rename_table", # humann_rename_table --input demo_fastq/rxn-cpm.tsv --output demo_fastq/rxn-cpm-named.tsv --names metacyc-rxn
                "seqkit",
                "diamond",
                
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
    
    step = 1

    # Info
    program = "preprocess_reads"
    program_label = "{}__{}".format(step, program)
    description = "Preprocessing input reads"
    
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # i/o
    if opts.input_reads_format == "paired":
        input_filepaths = [opts.forward_reads, opts.reverse_reads]
        output_filepaths = [
            os.path.join(output_directory, "joined.fastq.gz"),
            ]

    if opts.input_reads_format == "bam":
        input_filepaths = [opts.bam]
        output_filepaths = [
            os.path.join(output_directory, "joined.fastq.gz"),
        ]

    if opts.input_reads_format == "joined":
        input_filepaths = [opts.joined_reads]
        output_filepaths = [opts.joined_reads]

    output_filepaths += [os.path.join(output_directory, "reads.seqkit_stats.tsv")]


    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_preprocess_reads_cmd(**params)
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
    # Preprocess database
    # ==========
    
    step = 2

    # Info
    program = "diamond_database"
    program_label = "{}__{}".format(step, program)
    description = "Preprocessing Diamond database"
    
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # i/o
    if opts.database_format == "fasta":
        input_filepaths = [opts.fasta]
    if opts.database_format == "diamond_database":
        input_filepaths = [opts.diamond_database]

    output_filepaths = [
            os.path.join(output_directory,  "database{}".format(DIAMOND_DATABASE_SUFFIX)),
            os.path.join(output_directory,  "database{}.dbinfo.txt".format(DIAMOND_DATABASE_SUFFIX)),
        ]

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_preprocess_database_cmd(**params)
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
    # HUMAnN
    # ==========
    
    step = 3

    # Info
    program = "humann"
    program_label = "{}__{}".format(step, program)
    description = "HUMAnN pathway profiling"
    
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # i/o
    if opts.input_reads_format == "joined":
        input_filepaths = [opts.joined_reads]
    else:
        input_filepaths = [os.path.join(directories[("intermediate", "1__preprocess_reads")], "joined.fastq.gz")]
    input_filepaths += [
        os.path.join(directories[("intermediate", "2__diamond_database")],  "database{}".format(DIAMOND_DATABASE_SUFFIX)),
        opts.identifier_mapping,
    ]


    output_filenames =  [
        "humann_pathabundance.tsv",
        "humann_pathcoverage.tsv",
        "humann_genefamilies.tsv",
        "humann.log.gz",
        "humann_diamond_unaligned.fa.gz",
        "humann_diamond_aligned.tsv.gz",
    ]

    output_filepaths = list(map(lambda fn:os.path.join(output_directory, fn), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_humann_cmd(**params)
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
    step = 4

    # Info
    program = "symlink"
    program_label = "{}__{}".format(step, program)
    description = "Symlinking relevant output files"

    # Add to directories
    output_directory = directories["output"]

    # i/o
    input_filepaths = [
            os.path.join(directories[("intermediate", "1__preprocess_reads")], "reads.seqkit_stats.tsv"),
            os.path.join(directories[("intermediate", "3__humann")], "humann_pathabundance.tsv"),
            os.path.join(directories[("intermediate", "3__humann")], "humann_pathcoverage.tsv"),
            os.path.join(directories[("intermediate", "3__humann")], "humann_genefamilies.tsv"),
            os.path.join(directories[("intermediate", "3__humann")], "humann_diamond_unaligned.fa.gz"),
            os.path.join(directories[("intermediate", "3__humann")], "humann_diamond_aligned.tsv.gz"),
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

    # --input_reads_format
    assert_acceptable_arguments(opts.input_reads_format, {"paired", "joined", "bam", "auto"})
    if opts.input_reads_format == "auto":
        if any([opts.forward_reads, opts.reverse_reads]):
            assert opts.forward_reads != opts.reverse_reads, "You probably mislabeled the input files because `forward_reads` should not be the same as `reverse_reads`: {}".format(opts.forward_reads)
            assert opts.forward_reads is not None, "If running in --input_reads_format paired mode, --forward_reads and --reverse_reads are needed."
            assert opts.reverse_reads is not None, "If running in --input_reads_format paired mode, --forward_reads and --reverse_reads are needed."
            opts.input_reads_format = "paired"
        if opts.joined_reads is not None:
            assert opts.forward_reads is None, "If running in --input_reads_format joined mode, you cannot provide --forward_reads, --reverse_reads, or --bam."
            assert opts.reverse_reads is None, "If running in --input_reads_format joined mode, you cannot provide --forward_reads, --reverse_reads, or --bam."
            assert opts.bam is None, "If running in --input_reads_format joined mode, you cannot provide --forward_reads, --reverse_reads, or --bam."
            opts.input_reads_format = "joined"
        if opts.bam is not None:
            assert opts.forward_reads is None, "If running in --input_reads_format joined mode, you cannot provide --forward_reads, --reverse_reads, or --joined_reads."
            assert opts.reverse_reads is None, "If running in --input_reads_format joined mode, you cannot provide --forward_reads, --reverse_reads, or --joined_reads."
            assert opts.joined_reads is None, "If running in --input_reads_format joined mode, you cannot provide --forward_reads, --reverse_reads, or --joined_reads."
            opts.input_reads_format = "bam"
        print("Auto detecting reads format: {}".format(opts.input_reads_format), file=sys.stdout)
    assert_acceptable_arguments(opts.input_reads_format, {"paired", "joined", "bam"})

    # --search_mode
    assert_acceptable_arguments(opts.search_mode, {"uniref50", "uniref90", "auto"})
    if opts.search_mode == "auto":
        if opts.identifier_mapping.endswith(".gz"):
            f = gzip.open(opts.identifier_mapping, "rt")
        else:
            f = open(opts.identifier_mapping, "r")
        line = f.readline().strip('\n')
        f.close()

        fields = line.split("\t")
        id_uniref = fields[1]
        assert id_uniref.startswith("UniRef"), "--identifier_mapping is not formatted correctly.  The 2nd column should have UniRef identifiers: {}".format(id_uniref)
        opts.search_mode = id_uniref.split("_")[0].lower()
    assert_acceptable_arguments(opts.search_mode, {"uniref50", "uniref90"})

    # --database and --fasta
    assert any([opts.fasta, opts.diamond_database]), "Either --diamond_database (preferred) or --fasta must be provided but not both"
    assert not all([opts.fasta, opts.diamond_database]), "Either --diamond_database (preferred) or --fasta must be provided but not both"
    if opts.fasta:
        opts.database_format = "fasta"
    if opts.diamond_database:
        opts.database_format = "diamond_database"

    # --pathways
    assert_acceptable_arguments(opts.pathways, {"metacyc", "unipathway"})

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
    parser_reads = parser.add_argument_group('Required reads arguments')
    parser_reads.add_argument("-1","--forward_reads", type=str, help = "path/to/forward_reads.fq (Requires --reverse_reads, cannot be used with --joined_reads or --bam)")
    parser_reads.add_argument("-2","--reverse_reads", type=str, help = "path/to/reverse_reads.fq (Requires --forward_reads, cannot be used with --joined_reads or --bam)")
    parser_reads.add_argument("-j","--joined_reads", type=str, help = "path/to/joined_reads.fq (e.g., bbmerge.sh output) (Cannot be used with --forward_reads, --reverse_reads, or --bam)")
    parser_reads.add_argument("-b","--bam", type=str, help = "path/to/mapped.sorted.bam file aligned to genomes (Cannot be used with --forward_reads, --reverse_reads, or --joined_reads)")
    parser_reads.add_argument("-F", "--input_reads_format", type=str, default="auto", help = "Input reads format {paired, joined, bam} [Default: auto]")

    parser_database = parser.add_argument_group('Required database arguments')
    parser_database.add_argument("-i", "--identifier_mapping", type=str, required=True, help = "Identifier mapping which includes [id_protein]<tab>[id_uniref]<tab>[length]<tab>[lineage].  In VEBA, you can use `compile_custom_humann_database_from_annotations.py`. \nhttps://github.com/biobakery/humann#custom-reference-database-annotations ")
    parser_database.add_argument("-f", "--fasta", type=str, help = "Protein fasta to build database")
    parser_database.add_argument("-d","--diamond_database", type=str, help = "Diamond database with all proteins from --identifier_mapping") #! Future versions allow multiple databases


    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    parser_io.add_argument("-o","--project_directory", type=str, default="veba_output/profiling/pathways", help = "path/to/project_directory [Default: veba_output/profiling/pathways]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
    parser_utility.add_argument("--tmpdir", type=str, help="Set temporary directory")  #site-packges in future

    # BBMerge
    parser_bbmerge = parser.add_argument_group('bbmerge.sh arguments')
    parser_bbmerge.add_argument("--minimum_merge_overlap", type=int, default=12, help="bbmerge.sh | Minimum number of overlapping bases to allow merging. [Default: 12]")
    parser_bbmerge.add_argument("--bbmerge_options", type=str, default="", help="bbmerge.sh options (e.g. --arg 1 ) [Default: '']")

    # HUMAnN
    parser_humann = parser.add_argument_group('HUMAnN arguments')
    parser_humann.add_argument("--search_mode", type=str, default="auto", help="HUMAnN | Search for uniref50 or uniref90 gene families {uniref50, uniref90, auto} [Default: 'auto']")
    parser_humann.add_argument("--pathways", type=str, default="metacyc", help="HUMAnN | The database to use for pathway computations {metacyc, unipathway} [Default: 'metacyc']")
    parser_humann.add_argument("-e", "--evalue", type=float, default=1.0, help="HUMAnN | The evalue threshold to use with the translated search [Default: 1.0]")
    parser_humann.add_argument("-m", "--translated_identity_threshold", type=float,  help="HUMAnN | Identity threshold for translated alignments [Default: Tuned automatically (based on uniref mode) unless a custom value is specified]")
    parser_humann.add_argument("-q", "--translated_query_coverage_threshold", type=float, default=90.0,  help="HUMAnN | Query coverage threshold for translated alignments [Default: 90.0]")
    parser_humann.add_argument("-s", "--translated_subject_coverage_threshold", type=float, default=50.0,  help="HUMAnN | Subject coverage threshold for translated alignments [Default: 50.0]")
    parser_humann.add_argument("--humann_memory", type=str, default="minimum", help="HUMAnN | Memory use mode {minimum, maximum} [Default: 'minimum']")
    parser_humann.add_argument("--humann_options", type=str, default="", help="HUMAnN options (e.g. --arg 1 ) [Default: '']")

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
