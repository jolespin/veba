#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
from collections import OrderedDict, defaultdict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.2.1"

# Global clustering
def get_global_clustering_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command

    cmd = [
        os.environ["global_clustering.py"],
        "-i {}".format(input_filepaths[0]),
        "-o {}".format(output_directory),
        "--no_singletons" if bool(opts.no_singletons) else "",
        "-p {}".format(opts.n_jobs),

        "--ani_threshold {}".format(opts.ani_threshold),
        "--genome_cluster_prefix {}".format(opts.genome_cluster_prefix) if bool(opts.genome_cluster_prefix) else "",
        "--genome_cluster_suffix {}".format(opts.genome_cluster_suffix) if bool(opts.genome_cluster_suffix) else "",
        "--genome_cluster_prefix_zfill {}".format(opts.genome_cluster_prefix_zfill) if bool(opts.genome_cluster_prefix_zfill) else "",
        "--fastani_options {}".format(opts.fastani_options) if bool(opts.fastani_options) else "",

        "--minimum_identity_threshold {}".format(opts.minimum_identity_threshold),
        "--minimum_coverage_threshold {}".format(opts.minimum_coverage_threshold),
        "--protein_cluster_prefix {}".format(opts.protein_cluster_prefix) if bool(opts.protein_cluster_prefix) else "",
        "--protein_cluster_suffix {}".format(opts.protein_cluster_suffix) if bool(opts.protein_cluster_suffix) else "",
        "--protein_cluster_prefix_zfill {}".format(opts.protein_cluster_prefix_zfill) if bool(opts.protein_cluster_prefix_zfill) else "",
        "--mmseqs2_options {}".format(opts.mmseqs2_options) if bool(opts.mmseqs2_options) else "",

            "&&",

        "ln -sf {} {}".format(os.path.realpath(os.path.join(output_directory, "output")), os.path.join(directories["output"], "global")),

    ]

    return cmd

# Local clustering
def get_local_clustering_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command

    cmd = [
        os.environ["local_clustering.py"],
        "-i {}".format(input_filepaths[0]),
        "-o {}".format(output_directory),
        "--no_singletons" if bool(opts.no_singletons) else "",
        "-p {}".format(opts.n_jobs),

        "--ani_threshold {}".format(opts.ani_threshold),
        "--genome_cluster_prefix {}".format(opts.genome_cluster_prefix) if bool(opts.genome_cluster_prefix) else "",
        "--genome_cluster_suffix {}".format(opts.genome_cluster_suffix) if bool(opts.genome_cluster_suffix) else "",
        "--genome_cluster_prefix_zfill {}".format(opts.genome_cluster_prefix_zfill) if bool(opts.genome_cluster_prefix_zfill) else "",
        "--fastani_options {}".format(opts.fastani_options) if bool(opts.fastani_options) else "",

        "--minimum_identity_threshold {}".format(opts.minimum_identity_threshold),
        "--minimum_coverage_threshold {}".format(opts.minimum_coverage_threshold),
        "--protein_cluster_prefix {}".format(opts.protein_cluster_prefix) if bool(opts.protein_cluster_prefix) else "",
        "--protein_cluster_suffix {}".format(opts.protein_cluster_suffix) if bool(opts.protein_cluster_suffix) else "",
        "--protein_cluster_prefix_zfill {}".format(opts.protein_cluster_prefix_zfill) if bool(opts.protein_cluster_prefix_zfill) else "",
        "--mmseqs2_options {}".format(opts.mmseqs2_options) if bool(opts.mmseqs2_options) else "",


            "&&",

        "ln -sf {} {}".format(os.path.realpath(os.path.join(output_directory, "output")), os.path.join(directories["output"], "local")),

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
        "global_clustering.py",
        "local_clustering.py",
        "edgelist_to_clusters.py",
        }

    required_executables={
                # 1
                "fastANI",
                "mmseqs",
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
    # Global clustering
    # ==========
    step = 1

    program = "global_clustering"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Global clustering of genomes (FastANI) and proteins (MMSEQS2)"

    # i/o
    input_filepaths = [opts.input]

    output_filenames = [
        "output/*.tsv",
        ]

    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_global_clustering_cmd(**params)

    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
                errors_ok=False,
    )

    # ==========
    # Local clustering
    # ==========
    if not opts.no_local_clustering:
        step = 2

        program = "local_clustering"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

        # Info
        description = "Local clustering of genomes (FastANI) and proteins (MMSEQS2)"

        # i/o
        input_filepaths = [opts.input]

        output_filenames = [
            "output/*.tsv",
            ]

        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_local_clustering_cmd(**params)

        pipeline.add_step(
                    id=program,
                    description = description,
                    step=step,
                    cmd=cmd,
                    input_filepaths = input_filepaths,
                    output_filepaths = output_filepaths,
                    validate_inputs=True,
                    validate_outputs=True,
                    errors_ok=False,
        )




    return pipeline

# Configure parameters
def configure_parameters(opts, directories):
    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -m <mags> -a <proteins> -o <output_directory> -t 95".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i", "--input", type=str, required=True,  help = "path/to/input.tsv, Format: Must include the follow columns (No header) [organism_type]<tab>[id_sample]<tab>[id_mag]<tab>[genome]<tab>[proteins] but can include additional columns to the right (e.g., [cds]<tab>[gene_models]).  Suggested input is from `compile_genomes_table.py` script.")
    parser_io.add_argument("-o","--output_directory", type=str, default="veba_output/cluster", help = "path/to/project_directory [Default: veba_output/cluster]")
    parser_io.add_argument("-e", "--no_singletons", action="store_true", help="Exclude singletons") #isPSLC-1_SSO-3345__SRR178126
    parser_io.add_argument("--no_local_clustering", action="store_true", help = "Only do global clustering")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # FastANI
    parser_fastani = parser.add_argument_group('FastANI arguments')
    parser_fastani.add_argument("-a", "--ani_threshold", type=float, default=95.0, help="FastANI | Species-level cluster (SLC) ANI threshold (Range (0.0, 100.0]) [Default: 95.0]")
    parser_fastani.add_argument("--genome_cluster_prefix", type=str, default="SLC-", help="Cluster prefix [Default: 'SLC-")
    parser_fastani.add_argument("--genome_cluster_suffix", type=str, default="", help="Cluster suffix [Default: '")
    parser_fastani.add_argument("--genome_cluster_prefix_zfill", type=int, default=0, help="Cluster prefix zfill. Use 7 to match identifiers from OrthoFinder.  Use 0 to add no zfill. [Default: 0]") #7
    parser_fastani.add_argument("--fastani_options", type=str, default="", help="FastANI | More options (e.g. --arg 1 ) [Default: '']")


    # MMSEQS2
    parser_mmseqs2 = parser.add_argument_group('MMSEQS2 arguments')
    parser_mmseqs2.add_argument("-t", "--minimum_identity_threshold", type=float, default=50.0, help="MMSEQS2 | SLC-Specific Protein Cluster (SSPC, previously referred to as SSO) percent identity threshold (Range (0.0, 100.0]) [Default: 50.0]")
    parser_mmseqs2.add_argument("-c", "--minimum_coverage_threshold", type=float, default=0.8, help="MMSEQS2 | SSPC coverage threshold (Range (0.0, 1.0]) [Default: 0.8]")
    parser_mmseqs2.add_argument("--protein_cluster_prefix", type=str, default="SSPC-", help="Cluster prefix [Default: 'SSPC-")
    parser_mmseqs2.add_argument("--protein_cluster_suffix", type=str, default="", help="Cluster suffix [Default: '")
    parser_mmseqs2.add_argument("--protein_cluster_prefix_zfill", type=int, default=0, help="Cluster prefix zfill. Use 7 to match identifiers from OrthoFinder.  Use 0 to add no zfill. [Default: 0]") #7
    parser_mmseqs2.add_argument("--mmseqs2_options", type=str, default="", help="MMSEQS2 | More options (e.g. --arg 1 ) [Default: '']")

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
    print("Moment:", get_timestamp(), file=sys.stdout)
    print("Directory:", os.getcwd(), file=sys.stdout)
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
