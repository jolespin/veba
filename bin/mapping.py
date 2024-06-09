#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
from collections import OrderedDict, defaultdict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from genopype import __version__ as genopype_version
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.4.29"


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

    if opts.scaffolds_to_bins:
        cmd += [ 
            "&&",
            "(",
            os.environ["genome_spatial_coverage.py"],
            "-i {}".format(opts.scaffolds_to_bins),
            "-f {}".format(opts.reference_fasta),
            "-o {}".format(os.path.join(output_directory, "genome_spatial_coverage.tsv.gz")),
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
    "(",
        os.environ["featureCounts"],
        # "-G {}".format(opts.reference_fasta),
        "-a {}".format(opts.reference_gff),
        "-o {}".format(os.path.join(output_directory, "featurecounts.orfs.tsv")),
        "-F GTF",
        "--tmpDir {}".format(os.path.join(directories["tmp"], "featurecounts")),
        "-T {}".format(min(64, opts.n_jobs)), # The maximum number of threads featureCounts can use is 64 so any more will throw this error: "Value for argumant -T is out of range: 1 to 64"
        "-g {}".format(opts.attribute_type),
        "-t {}".format(opts.feature_type),
        "-p --countReadPairs",
        opts.featurecounts_options,
        input_filepaths[0],
    ")",


    # Scaffold-Level Counts
        "&&",

    "(",
        os.environ["featureCounts"],
        # "-G {}".format(opts.reference_fasta),
        "-a {}".format(opts.reference_saf),
        "-o {}".format(os.path.join(output_directory, "featurecounts.scaffolds.tsv")),
        "-F SAF",
        "--tmpDir {}".format(os.path.join(directories["tmp"], "featurecounts")),
        "-T {}".format(opts.n_jobs),
        "-p --countReadPairs",
        opts.featurecounts_options,
        input_filepaths[0],
    ")",
        "&&",
        "tail -n +3",
        os.path.join(output_directory, "featurecounts.orfs.tsv"),
        "|",
        "cut -f1,7",
        "|",
        "gzip",
        ">",
        os.path.join(output_directory, "counts.orfs.tsv.gz"),
        "&&",
        "tail -n +3",
        os.path.join(output_directory, "featurecounts.scaffolds.tsv"),
        "|",
        "cut -f1,7",
        "|",
        "gzip",
        ">",
        os.path.join(output_directory, "counts.scaffolds.tsv.gz"),
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

    if opts.proteins_to_orthogroups:
        cmd += [ 
            "&&",
            os.environ["groupby_table.py"],
            "-m {}".format(opts.proteins_to_orthogroups),
            "-t {}".format(os.path.join(output_directory, "counts.orfs.tsv.gz")),
            "-o {}".format(os.path.join(output_directory, "counts.orthogroups.tsv.gz")),
        ]

    if opts.scaffolds_to_bins:
        cmd += [ 
            "&&",
            os.environ["groupby_table.py"],
            "-m {}".format(opts.scaffolds_to_bins),
            "-t {}".format(os.path.join(output_directory, "counts.scaffolds.tsv.gz")),
            "-o {}".format(os.path.join(output_directory, "counts.mags.tsv.gz")),
        ]

    if opts.scaffolds_to_clusters:
        cmd += [ 
            "&&",
            os.environ["groupby_table.py"],
            "-m {}".format(opts.scaffolds_to_clusters),
            "-t {}".format(os.path.join(output_directory, "counts.scaffolds.tsv.gz")),
            "-o {}".format(os.path.join(output_directory, "counts.clusters.tsv.gz")),
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

    accessory_scripts = {"groupby_table.py", "genome_spatial_coverage.py"}

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
    if opts.scaffolds_to_bins:
        output_filenames += ["genome_spatial_coverage.tsv.gz"]
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

    output_filenames = ["counts.orfs.tsv.gz", "counts.scaffolds.tsv.gz"]
    if opts.proteins_to_orthogroups:
        output_filenames.append("counts.orthogroups.tsv.gz")
    if opts.scaffolds_to_bins:
        output_filenames.append("counts.mags.tsv.gz")
    if opts.scaffolds_to_clusters:
        output_filenames.append("counts.clusters.tsv.gz")
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
        os.path.join(directories[("intermediate", "2__featurecounts")], "counts.orfs.tsv.gz"),
        os.path.join(directories[("intermediate", "2__featurecounts")], "counts.scaffolds.tsv.gz"),
    ]

    if opts.proteins_to_orthogroups:
        input_filepaths += [ 
            os.path.join(directories[("intermediate", "2__featurecounts")], "counts.orthogroups.tsv.gz"),
        ]

    if opts.scaffolds_to_bins:
        input_filepaths += [ 
             
            os.path.join(directories[("intermediate", "1__bowtie2")], "genome_spatial_coverage.tsv.gz"),
            os.path.join(directories[("intermediate", "2__featurecounts")], "counts.mags.tsv.gz"),
        ]   

    if opts.scaffolds_to_clusters:
        input_filepaths += [ 
            os.path.join(directories[("intermediate", "2__featurecounts")], "counts.clusters.tsv.gz"),
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
    # Assert that if --reference_gff or --reference_saf are not provided, then --reference_index is a directory
    if opts.reference_gff is None:
        assert os.path.isdir(opts.reference_index), "If --reference_gff is not provided, then --reference_index must be provided as a directory containing a file 'reference.gff'"
        opts.reference_gff = os.path.join(opts.reference_index, "reference.gff")

    if opts.reference_saf is None:
        assert os.path.isdir(opts.reference_index), "If --reference_saf is not provided, then --reference_index must be provided as a directory containing a file 'reference.saf'"
        opts.reference_saf = os.path.join(opts.reference_index, "reference.saf")

    # Check if --reference_index is a directory, if it is then set reference.fa as the directory
    if os.path.isdir(opts.reference_index):
        if opts.reference_gzipped:
            opts.reference_index = os.path.join(opts.reference_index, "reference.fa.gz")
        else:
            opts.reference_index = os.path.join(opts.reference_index, "reference.fa")

    # If --reference_fasta isn't provided then set it to the --reference_index
    if opts.reference_fasta is None:
        opts.reference_fasta =  opts.reference_index

    if opts.proteins_to_orthogroups is not None:
        assert os.path.exists(opts.proteins_to_orthogroups)

    if opts.scaffolds_to_bins is not None:
        assert os.path.exists(opts.scaffolds_to_bins)

    if opts.scaffolds_to_clusters is not None:
        assert os.path.exists(opts.scaffolds_to_clusters)

    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -1 <reads_1.fq> -2 <reads_2.fq> -n <name> -o <output_directory> -x <reference_directory>"    .format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-1","--forward_reads", type=str, help = "path/to/reads_1.fastq", required=True)
    parser_io.add_argument("-2","--reverse_reads", type=str, help = "path/to/reads_2.fastq", required=True)
    parser_io.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    parser_io.add_argument("-o","--project_directory", type=str, default="veba_output/mapping", help = "path/to/project_directory [Default: veba_output/mapping]")

    parser_reference = parser.add_argument_group('Reference arguments')
    parser_reference.add_argument("-x", "--reference_index",type=str, required=True, help="path/to/bowtie2_index. Either a file or directory. If directory, then it assumes the index is named `reference.fa`")
    parser_reference.add_argument("-r", "--reference_fasta", type=str, required=False, help = "path/to/reference.fasta. If not provided then it is set to the --reference_index" ) # ; or (2) a directory of fasta files [Must all have the same extension.  Use `query_ext` argument]
    parser_reference.add_argument("-a", "--reference_gff",type=str, required=False, help="path/to/reference.gff. If not provided then --reference_index must be a directory that contains the file: 'reference.gff'")
    parser_reference.add_argument("-s", "--reference_saf",type=str, required=False, help="path/to/reference.saf. If not provided then --reference_index must be a directory that contains the file: 'reference.saf'")
    parser_reference.add_argument("-z", "--reference_gzipped",action="store_true", help="If --reference_index directory, then it assumes the index is named `reference.fa.gz` instead of `reference.fa`")

    # parser_io.add_argument("-S","--scaffold_identifier_mapping", type=str, required=False,  help = "path/to/scaffold_identifiers.tsv, Format: [id_scaffold]<tab>[id_mag]<tab>[id_cluster], No header")
    # parser_io.add_argument("-O","--orf_identifier_mapping", type=str, required=False,  help = "path/to/scaffold_identifiers.tsv, Format: [id_scaffold]<tab>[id_mag]<tab>[id_cluster], No header")

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
    #! NOT SURE HOW THIS WILL WORK WITH PRODIGAL AND METAEUK. WILL PROBABLY NEED TO POST PROCESS METAEUK
    parser_featurecounts = parser.add_argument_group('featureCounts arguments')
    parser_featurecounts.add_argument("-g", "--attribute_type", type=str, default="gene_id", help = "Attribute type in GTF/GFF file. [Default: gene_id]")
    parser_featurecounts.add_argument("-t", "--feature_type", type=str, default="CDS", help = "Feature type in GTF/GFF file. [Default: CDS]")
    parser_featurecounts.add_argument("--retain_featurecounts", default=0, type=int, help = "Retain feature counts output table (a slimmer version is output regardless). 0=No, 1=yes [Default: 0]") 
    # parser_featurecounts.add_argument("--long_reads", action="store_true", help="featureCounts | Use this if long reads are being used")
    parser_featurecounts.add_argument("--featurecounts_options", type=str, default="", help="featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/")

    parser_identifiers = parser.add_argument_group('Identifier arguments')
    parser_identifiers.add_argument("--proteins_to_orthogroups", type=str, help = "path/to/protein_to_orthogroup.tsv, [id_orf]<tab>[id_orthogroup], No header")
    parser_identifiers.add_argument("--scaffolds_to_bins", type=str, help = "path/to/scaffold_to_bins.tsv, [id_scaffold]<tab>[id_bin], No header")
    parser_identifiers.add_argument("--scaffolds_to_clusters", type=str, help = "path/to/scaffold_to_cluster.tsv, [id_scaffold]<tab>[id_cluster], No header")

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
    print("GenoPype version:", genopype_version, file=sys.stdout) #sys.path[2]
    print("Script version:", __version__, file=sys.stdout)
    print("Moment:", get_timestamp(), file=sys.stdout)
    print("Directory:", os.getcwd(), file=sys.stdout)
    if "TMPDIR" in os.environ: print(os.environ["TMPDIR"], file=sys.stdout)
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
