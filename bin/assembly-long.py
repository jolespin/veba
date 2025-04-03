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
__version__ = "2024.12.11"

# Assembly
def get_assembly_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    cmd = [
        os.environ["flye"],
        "--{} {}".format(opts.reads_type, input_filepaths[0]),
        "-g {}".format(opts.estimated_assembly_size) if opts.estimated_assembly_size else "",
        "-o {}".format(output_directory),
        "-t {}".format(opts.n_jobs),
        "--deterministic" if opts.deterministic else "",
        "--meta" if opts.program == "metaflye" else "",
        opts.assembler_options,

            # Get failed length cutoff fasta
                "&&",

            "mv",
            os.path.join(output_directory, "assembly.fasta"),
            os.path.join(output_directory, "assembly_original.fasta"),

                "&&",

            "cat",
            os.path.join(output_directory, "assembly_original.fasta"),
            "|",
            os.environ["seqkit"],
            "seq",
            "-M {}".format(max(opts.minimum_contig_length - 1, 1)),
            "|",
            "gzip",
            ">",
            os.path.join(output_directory, "assembly_failed_length_cutoff.fasta.gz"),

            # Filter out small scaffolds and add prefix if applicable
                "&&",

            "cat",
            os.path.join(output_directory, "assembly_original.fasta"),
            "|",
            os.environ["seqkit"],
            "seq",
            "-m {}".format(opts.minimum_contig_length),
            "|",
            os.environ["seqkit"],
            "replace",
            "-r {}".format(opts.scaffold_prefix),
            "-p '^'",
            ">",
            os.path.join(output_directory, "assembly.fasta"),

                "&&",

            "rm -rf",
            os.path.join(output_directory, "assembly_original.fasta"),

                "&&",
                
            os.environ["fasta_to_saf.py"],
            "-i",
            os.path.join(output_directory, "assembly.fasta"),
            ">",
            os.path.join(output_directory, "assembly.fasta.saf"),
            
                "&&",
                
            os.environ["prepend_de-bruijn_path.py"],
            "-i",
            os.path.join(output_directory, "assembly_graph.gfa"),
            "--prefix",
            opts.scaffold_prefix,
            "-o",
            os.path.join(output_directory, "assembly_graph.prefixed.gfa"),
            "--program",
            "flye",
            
                "&&",
                
            "mv",
            os.path.join(output_directory, "assembly_graph.prefixed.gfa"),
            os.path.join(output_directory, "assembly_graph.gfa"),
            
                "&&",

            os.environ["prepend_de-bruijn_path.py"],
            "-i",
            os.path.join(output_directory, "assembly_info.txt"),
            "--prefix",
            opts.scaffold_prefix,
            "-o",
            os.path.join(output_directory, "assembly_info.prefixed.txt"),
            "--program",
            "flye",
            
                "&&",
                
            "mv",
            os.path.join(output_directory, "assembly_info.prefixed.txt"),
            os.path.join(output_directory, "assembly_info.txt"),

        ]



    # files_to_remove = [ 
    # ]

    # for fn in files_to_remove:
    #     cmd += [ 
    #         "&&",
    #         "rm -rf {}".format(os.path.join(output_directory, fn)),
    #     ]
    return cmd

# Bowtie2
def get_alignment_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [
    # Clear temporary directory just in case
    "rm -rf {}".format(os.path.join(directories["tmp"], "*")),
    "&&",

    # MiniMap2 Index
    "(",
    os.environ["minimap2"],
    "-t {}".format(opts.n_jobs),
    "-d {}".format(output_filepaths[0]), # Index
    opts.minimap2_index_options,
    input_filepaths[1], # Reference
    ")",

    "&&",

    # MiniMap2
    "(",
    os.environ["minimap2"],
    "-a",
    "-t {}".format(opts.n_jobs),
    "-x {}".format(opts.minimap2_preset),
    opts.minimap2_options,
    output_filepaths[0],
    input_filepaths[0],



    # Convert to sorted BAM
    "|",

    os.environ["samtools"],
    "view",
    "-b",
    "-h",
    "-F 4",

    "|",

    os.environ["samtools"],
    "sort",
    "--threads {}".format(opts.n_jobs),
    "--reference {}".format(input_filepaths[1]),
    "-T {}".format(os.path.join(directories["tmp"], "samtools_sort")),
    ">",
    output_filepaths[1],
    ")",

    "&&",

    "(",
    os.environ["samtools"],
    "index",
    "-@ {}".format(opts.n_jobs),
    output_filepaths[1],
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
        "-L",
        "--tmpDir {}".format(os.path.join(directories["tmp"], "featurecounts")),
        "-T {}".format(min(64, opts.n_jobs)), # The maximum number of threads featureCounts can use is 64 so any more will throw this error: "Value for argumant -T is out of range: 1 to 64"
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
    accessory_scripts = {
                "prepend_de-bruijn_path.py",
                "fasta_to_saf.py",
                }

    required_executables={
                "flye",
                "minimap2",
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
    # Assembly
    # ==========
    
    step = 1

    # Info
    program = "assembly"
    program_label = "{}__{}".format(step, program)
    description = "Assembling long reads via {}".format(opts.program.capitalize())
    
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # i/o
    input_filepaths = [opts.reads]
    output_filenames = ["assembly.fasta", "assembly.fasta.saf"]
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
        opts.reads,
        os.path.join(directories[("intermediate", "1__assembly")], "assembly.fasta"),
        ] 

    output_filepaths = [
        os.path.join(directories[("intermediate", "1__assembly")], "assembly.fasta.mmi"),
        os.path.join(output_directory, "mapped.sorted.bam"),
    ]

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

    input_filepaths = [ 
        os.path.join(directories[("intermediate", "1__assembly")], "assembly.fasta"),
        os.path.join(directories[("intermediate", "1__assembly")], "assembly.fasta.saf"),
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

    input_filepaths = [ 
        os.path.join(directories[("intermediate", "1__assembly")], "assembly.fasta"),
        os.path.join(directories[("intermediate", "1__assembly")], "assembly.fasta.mmi"),
        os.path.join(directories[("intermediate", "1__assembly")], "assembly_graph.gfa"),
        os.path.join(directories[("intermediate", "1__assembly")], "assembly_info.txt"),
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

    # Scaffold prefix
    if opts.scaffold_prefix == "NONE":
        opts.scaffold_prefix = ""
    else:
        if "NAME" in opts.scaffold_prefix:
            opts.scaffold_prefix = opts.scaffold_prefix.replace("NAME", opts.name)
        print("Using the following prefix for all {} scaffolds: {}".format(opts.program, opts.scaffold_prefix), file=sys.stdout)
    
    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <reads.fq[.gz]> -n <name> -g <estimated_genome_size> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i","--reads", type=str, required=True, help = "path/to/reads.fq[.gz]")
    parser_io.add_argument("-n", "--name", type=str, required=True, help="Name of sample")
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
    parser_assembler = parser.add_argument_group('Assembler arguments')
    parser_assembler.add_argument("-P", "--program", type=str, default="metaflye", choices={"flye", "metaflye"}, help="Assembler |  {flye, metaflye}} [Default: 'metaflye']")
    parser_assembler.add_argument("-s", "--scaffold_prefix", type=str, default="NAME__", help="Assembler |  Special options:  Use NAME to use --name.  Use NONE to not include a prefix. [Default: 'NAME__']")
    parser_assembler.add_argument("-m", "--minimum_contig_length", type=int, default=1, help="Minimum contig length.  Should be lenient here because longer thresholds can be used for binning downstream. Recommended for metagenomes to use 1000 here. [Default: 1] ")
    parser_assembler.add_argument("-t", "--reads_type", type=str, default="nano-hq", choices={"nano-hq", "nano-corr", "nano-raw", "pacbio-hifi", "pacbio-corr", "pacbio-raw"}, help="Reads type for (meta)flye.  {nano-hq, nano-corr, nano-raw, pacbio-hifi, pacbio-corr, pacbio-raw} [Default: nano-hq] ")
    parser_assembler.add_argument("-g", "--estimated_assembly_size", type=str,  help="Estimated assembly size (e.g., 5m, 2.6g)")
    parser_assembler.add_argument("--deterministic", action="store_true", help="Use deterministic mode.  This will result in a slower assembly and will not be threaded but should produce the same assembly each run.")
    parser_assembler.add_argument("--assembler_options", type=str, default="", help="Assembler options for Flye-based programs (e.g. --arg 1 ) [Default: '']")

    # Aligner
    parser_aligner = parser.add_argument_group('MiniMap2 arguments')
    parser_aligner.add_argument("--minimap2_preset", type=str, default="map-ont", help="MiniMap2 | MiniMap2 preset {map-pb, map-ont, map-hifi} [Default: map-ont]")
    # parser_aligner.add_argument("--no_create_index", action="store_true", help="Do not create a MiniMap2 index")
    parser_aligner.add_argument("--minimap2_index_options", type=str, default="", help="MiniMap2 | More options (e.g. --arg 1 ) [Default: '']\nhttps://github.com/lh3/minimap2")
    parser_aligner.add_argument("--minimap2_options", type=str, default="", help="MiniMap2 | More options (e.g. --arg 1 ) [Default: '']\nhttps://github.com/lh3/minimap2")

    # featureCounts
    parser_featurecounts = parser.add_argument_group('featureCounts arguments')
    parser_featurecounts.add_argument("--featurecounts_options", type=str, default="", help="featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/")


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
    # os.environ["TMPDIR"] = directories["tmp"]

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
