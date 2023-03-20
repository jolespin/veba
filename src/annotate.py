#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob,  shutil
from collections import OrderedDict, defaultdict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *
from soothsayer_utils.soothsayer_utils import assert_acceptable_arguments

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.3.14"


# Diamond
def get_diamond_nr_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    tmp = os.path.join(directories["tmp"], "diamond-nr")
    # Command
    cmd = [

        "mkdir -p {}".format(tmp),

            "&&",

        os.environ["diamond"],
        "blastp",
        "--db {}".format(input_filepaths[1]),
        "--query {}".format(input_filepaths[0]),
        "--threads {}".format(opts.n_jobs),
        "-f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle",
        "--evalue {}".format(opts.diamond_evalue),
        "-o {}".format(os.path.join(output_directory, "output.tsv")),
        "--max-target-seqs 1",
        "--tmpdir {}".format(tmp),
    ]
    if bool(opts.diamond_sensitivity):
        cmd += [ 
            "--{}".format(opts.diamond_sensitivity),
        ]
    cmd += [ 
        opts.diamond_options,
    ]

    cmd += [ 
            "&&",
            
        "pigz",
        "-f",
        "-p {}".format(opts.n_jobs),
        os.path.join(output_directory, "output.tsv"),
    ]           

    return cmd

def get_diamond_mibig_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    tmp = os.path.join(directories["tmp"], "diamond-mibig")
    # Command
    cmd = [

        "mkdir -p {}".format(tmp),

            "&&",

        os.environ["diamond"],
        "blastp",
        "--db {}".format(input_filepaths[1]),
        "--query {}".format(input_filepaths[0]),
        "--threads {}".format(opts.n_jobs),
        "-f 6 qseqid sseqid pident length mismatch qlen qstart qend slen sstart send evalue bitscore qcovhsp scovhsp",
        "--evalue {}".format(opts.diamond_evalue),
        "-o {}".format(os.path.join(output_directory, "output.tsv")),
        "--max-target-seqs 1",
        "--tmpdir {}".format(tmp),
    ]
    if bool(opts.diamond_sensitivity):
        cmd += [ 
            "--{}".format(opts.diamond_sensitivity),
        ]
    cmd += [ 
        opts.diamond_options,
    ]

    cmd += [ 
            "&&",
            
        "pigz",
        "-f",
        "-p {}".format(opts.n_jobs),
        os.path.join(output_directory, "output.tsv"),
    ]           

    return cmd


# HMMER
def get_hmmsearch_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [
        "(",
        os.environ["hmmsearch"],
        "--tblout {}".format(os.path.join(output_directory, "output.tsv")),
        "--cut_ga",
        "--cpu {}".format(opts.n_jobs),
        "--seed {}".format(opts.random_state + 1),
        input_filepaths[1],
        input_filepaths[0],
        ">",
        "/dev/null",
        ")",

            "&&",

        "pigz",
        "-f",
        "-p {}".format(opts.n_jobs),
        os.path.join(output_directory, "output.tsv"),
    ]    
    return cmd

# KOFAMSCAN
def get_kofamscan_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [
        "(",
        "mkdir -p {}".format(os.path.join(directories["tmp"], "kofamscan")),

            "&&",

        os.environ["exec_annotation"],
        "--cpu {}".format(opts.n_jobs),
        "-f detail-tsv",
        "-p {}".format(input_filepaths[1]),
        "-k {}".format(input_filepaths[2]),
        "--tmp-dir {}".format(os.path.join(directories["tmp"], "kofamscan")),
        opts.kofamscan_options,
        input_filepaths[0],
        "|",
        'grep "*"',
        ">",
        os.path.join(output_directory, "output.tsv"),
        ")",

            "&&",

        "pigz",
        "-f",
        "-p {}".format(opts.n_jobs),
        os.path.join(output_directory, "output.tsv"),
    ]    
    return cmd

def get_merge_and_score_taxonomy_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [
        os.environ["merge_annotations_and_score_taxonomy.py"],
        "--diamond_nr {}".format(input_filepaths[0]),
        "--diamond_mibig {}".format(input_filepaths[1]),
        "--hmmsearch_pfam {}".format(input_filepaths[2]),
        "--hmmsearch_amr {}".format(input_filepaths[3]),
        "--hmmsearch_antifam {}".format(input_filepaths[4]),
        "--kofam {}".format(input_filepaths[5]),
        "--veba_database {}".format(opts.veba_database),
        "--fasta {}".format(opts.proteins),
        "-o {}".format(output_directory)
    ]
    if opts.identifier_mapping:
        cmd += [ 
            "-i {}".format(opts.identifier_mapping),
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
                "merge_annotations_and_score_taxonomy.py",
                }

    required_executables={
                # 1
                "diamond",
                # 2 
                "hmmsearch",
                # 3
                "exec_annotation",

     } | accessory_scripts

    if opts.path_config == "CONDA_PREFIX":
        executables = dict()
        for name in sorted(required_executables):
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
    pipeline = ExecutablePipeline(name=__program__,  f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])

    # ==========
    # Diamond [nr]
    # ==========
    step = 1

    program = "diamond-nr"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Diamond [nr]"

    # i/o
    input_filepaths = [
        opts.proteins, 
         os.path.join(opts.veba_database, "Annotate", "nr", "nr.dmnd"),
        ]
    output_filenames = ["output.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_diamond_nr_cmd(**params)

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
    # Diamond
    # ==========
    step = 2

    program = "diamond-mibig"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Diamond [MIBiG]"

    # i/o
    input_filepaths = [
        opts.proteins, 
         os.path.join(opts.veba_database, "Annotate", "MIBiG", "mibig_v3.1.dmnd"),
        ]
    output_filenames = ["output.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_diamond_mibig_cmd(**params)

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
    # HMMSearch-Pfam
    # ==========
    step = 3

    program = "hmmsearch-pfam"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "HMMSearch [PFAM]"

    # i/o
    input_filepaths = [
        opts.proteins, 
        os.path.join(opts.veba_database, "Annotate", "Pfam", "Pfam-A.hmm.gz"),
        ]
    output_filenames = ["output.tsv.gz"]
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

    # =============
    # HMMSearch-AMR
    # =============
    step = 4

    program = "hmmsearch-amr"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "HMMSearch [NCBIfam-AMR]"

    # i/o
    input_filepaths = [
        opts.proteins, 
        os.path.join(opts.veba_database, "Annotate", "NCBIfam-AMRFinder", "NCBIfam-AMRFinder.hmm.gz"),
        ]
    output_filenames = ["output.tsv.gz"]
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

    # =================
    # HMMSearch-AntiFam
    # =================
    step = 5

    program = "hmmsearch-antifam"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "HMMSearch [AntiFam]"

    # i/o
    input_filepaths = [
        opts.proteins, 
        os.path.join(opts.veba_database, "Contamination", "AntiFam", "AntiFam.hmm.gz"),
        ]
    output_filenames = ["output.tsv.gz"]
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
    # ==========
    # KOFAMSCAN
    # ==========
    step = 6

    program = "kofamscan"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "KOFAMSCAN [KOFAM]"
    # i/o
    input_filepaths = [
        opts.proteins, 
        os.path.join(opts.veba_database, "Annotate", "KOFAM", "profiles"),
        os.path.join(opts.veba_database, "Annotate", "KOFAM", "ko_list"),
        ]
    output_filenames = ["output.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_kofamscan_cmd(**params)

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
    # Output
    # ==========
    step = 7

    program = "merge_and_score_taxonomy"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories["output"]

    # Info
    description = "Merging annotation results and scoring taxonomy [NCBI Taxonomy]"

    # i/o
    input_filepaths = [
        os.path.join(directories[("intermediate",  "1__diamond-nr")], "output.tsv.gz"),
        os.path.join(directories[("intermediate",  "2__diamond-mibig")], "output.tsv.gz"),
        os.path.join(directories[("intermediate",  "3__hmmsearch-pfam")], "output.tsv.gz"),
        os.path.join(directories[("intermediate",  "4__hmmsearch-amr")], "output.tsv.gz"),
        os.path.join(directories[("intermediate",  "5__hmmsearch-antifam")], "output.tsv.gz"),
        os.path.join(directories[("intermediate",  "6__kofamscan")], "output.tsv.gz"),
        os.path.join(opts.veba_database, "Classify", "NCBITaxonomy", "nodes.dmp"),
        os.path.join(opts.veba_database, "Classify", "NCBITaxonomy", "names.dmp"),
        os.path.join(opts.veba_database, "Classify", "NCBITaxonomy", "merged.dmp"),


    ]
    output_filenames = ["protein_annotations.tsv.gz"]
    if opts.identifier_mapping:
        output_filenames += ["lineage.weighted_majority_vote.contigs.tsv.gz", "lineage.weighted_majority_vote.genomes.tsv.gz"]
    
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_merge_and_score_taxonomy_cmd(**params)

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

    assert_acceptable_arguments(opts.diamond_sensitivity, {"", "fast", "mid-sensitive", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"})

    # Set environment variables
    add_executables_to_environment(opts=opts)

    # Concatenate proteins if a directory is provided
    if os.path.isdir(opts.proteins):
        fp_protein_cat = os.path.join(directories["tmp"], "proteins.faa")
        with open(fp_protein_cat, "w") as f_out:
            for fp in glob.glob(os.path.join(opts.proteins,"*.{}".format(opts.extension))):
                with open(fp, "r") as f_in:
                    shutil.copyfileobj(f_in, f_out)
        print("Concatenating files from --proteins {}".format(opts.proteins), file=sys.stdout)
        print("Changing --proteins to {}".format(fp_protein_cat), file=sys.stdout)
        opts.proteins = fp_protein_cat


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -a <proteins> -o <output_directory> |Optional: -i <identifier_mapping> ".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-a","--proteins", type=str, required=True, help = "Either path/to/proteins.faa or a directory of fasta files using [-x]")
    parser_io.add_argument("-o","--output_directory", type=str, default="veba_output/annotation", help = "path/to/project_directory [Default: veba_output/annotation]")
    parser_io.add_argument("-i","--identifier_mapping", type=str, required=False, help = "Tab-seperated value table of [id_orf]<tab>[id_contig]<tab>[id_mag]")
    parser_io.add_argument("-x", "--extension", type=str, default="faa", help = "Fasta file extension for proteins if a directory is provided for --proteins [Default: faa]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("--keep_temporary_directory", action="store_true",  help = "Keep temporary directory [Default is to remove]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Databases
    parser_databases = parser.add_argument_group('Database arguments')
    parser_databases.add_argument("--veba_database", type=str,  help=f"VEBA database location.  [Default: $VEBA_DATABASE environment variable]")
    
    # Diamond
    parser_diamond = parser.add_argument_group('Diamond arguments')
    parser_diamond.add_argument("--diamond_sensitivity", type=str, default="", help="Diamond | Sensitivity [Default:  '']")
    parser_diamond.add_argument("--diamond_evalue", type=float, default=0.001, help="Diamond | E-Value [Default: 0.001]")
    parser_diamond.add_argument("--diamond_options", type=str, default="", help="Diamond | More options (e.g. --arg 1 ) [Default: '']")

    # HMMER
    parser_hmmer = parser.add_argument_group('HMMSearch arguments')
    # parser_hmmer.add_argument("--hmmsearch_threshold", type=str, default="cut_ga", help="HMMSearch | Threshold {cut_ga, cut_nc, gut_tc} [Default:  cut_ga]")
    # parser_hmmer.add_argument("--hmmsearch_evalue", type=float, default=10.0, help="Diamond | E-Value [Default: 10.0]")
    parser_hmmer.add_argument("--hmmsearch_options", type=str, default="", help="Diamond | More options (e.g. --arg 1 ) [Default: '']")

    # KOFAMSCAN
    parser_kofamscan = parser.add_argument_group('KOFAMSCAN arguments')
    parser_kofamscan.add_argument("--kofamscan_options", type=str, default="", help="Diamond | More options (e.g. --arg 1 ) [Default: '']")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Threads
    if opts.n_jobs == -1:
        from multiprocessing import cpu_count 
        opts.n_jobs = cpu_count()
    assert opts.n_jobs >= 1, "--n_jobs must be ≥ 1.  To select all available threads, use -1."

    # Database
    if opts.veba_database is None:
        assert "VEBA_DATABASE" in os.environ, "Please set the following environment variable 'export VEBA_DATABASE=/path/to/veba_database' or provide path to --veba_database"
        opts.veba_database = os.environ["VEBA_DATABASE"]

    # Directories
    directories = dict()
    directories["project"] = create_directory(opts.output_directory)
    directories["output"] = create_directory(os.path.join(directories["project"], "output"))
    directories["log"] = create_directory(os.path.join(directories["project"], "log"))
    directories["tmp"] = create_directory(os.path.join(directories["project"], "tmp"))
    directories["checkpoints"] = create_directory(os.path.join(directories["project"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["project"], "intermediate"))
    os.environ["TMPDIR"] = directories["tmp"]
    
    # Info
    print(format_header(__program__, "="), file=sys.stdout)
    print(format_header("Configuration:", "-"), file=sys.stdout)
    print("Python version:", sys.version.replace("\n"," "), file=sys.stdout)
    print("Python path:", sys.executable, file=sys.stdout) #sys.path[2]
    print("Script version:", __version__, file=sys.stdout)
    print("VEBA Database:", opts.veba_database, file=sys.stdout)
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

    if not opts.keep_temporary_directory: # Add this to genopype
        for fp in glob.glob(os.path.join(directories["tmp"],"*")):
            if os.path.isdir(fp):
                shutil.rmtree(fp, ignore_errors=True)
            else:
                os.remove(fp)

if __name__ == "__main__":
    main()
