#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, warnings
from collections import OrderedDict, defaultdict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from genopype import __version__ as genopype_version
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.12.11"

# Global clustering
def get_global_clustering_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command

    cmd = [
        os.environ["global_clustering.py"],
        "-i {}".format(input_filepaths[0]),
        "-o {}".format(output_directory),
        # "--no_singletons" if bool(opts.no_singletons) else "",
        "-p {}".format(opts.n_jobs),

        "--genome_clustering_algorithm {}".format(opts.genome_clustering_algorithm),
        "--ani_threshold {}".format(opts.ani_threshold),
        "--genome_cluster_prefix {}".format(opts.genome_cluster_prefix) if bool(opts.genome_cluster_prefix) else "",
        "--genome_cluster_suffix {}".format(opts.genome_cluster_suffix) if bool(opts.genome_cluster_suffix) else "",
        "--genome_cluster_prefix_zfill {}".format(opts.genome_cluster_prefix_zfill) if bool(opts.genome_cluster_prefix_zfill) else "",
        "--skani_target_ani {}".format(opts.skani_target_ani),
        "--skani_minimum_af {}".format(opts.skani_minimum_af),
        "--skani_no_confidence_interval" if opts.skani_no_confidence_interval else "",

        "--skani_nonviral_preset {}".format(opts.skani_nonviral_preset),
        "--skani_nonviral_compression_factor {}".format(opts.skani_nonviral_compression_factor),
        "--skani_nonviral_marker_kmer_compression_factor {}".format(opts.skani_nonviral_marker_kmer_compression_factor),
        "--skani_nonviral_options {}".format(opts.skani_nonviral_options) if bool(opts.skani_nonviral_options) else "",

        "--skani_viral_preset {}".format(opts.skani_viral_preset),
        "--skani_viral_compression_factor {}".format(opts.skani_viral_compression_factor),
        "--skani_viral_marker_kmer_compression_factor {}".format(opts.skani_viral_marker_kmer_compression_factor),
        "--skani_viral_options {}".format(opts.skani_viral_options) if bool(opts.skani_viral_options) else "",

        "--fastani_options {}".format(opts.fastani_options) if bool(opts.fastani_options) else "",

        "--protein_clustering_algorithm {}".format(opts.protein_clustering_algorithm),
        "--minimum_identity_threshold {}".format(opts.minimum_identity_threshold),
        "--minimum_coverage_threshold {}".format(opts.minimum_coverage_threshold),
        "--protein_cluster_prefix {}".format(opts.protein_cluster_prefix) if bool(opts.protein_cluster_prefix) else "",
        "--protein_cluster_suffix {}".format(opts.protein_cluster_suffix) if bool(opts.protein_cluster_suffix) else "",
        "--protein_cluster_prefix_zfill {}".format(opts.protein_cluster_prefix_zfill) if bool(opts.protein_cluster_prefix_zfill) else "",
        "--mmseqs2_options {}".format(opts.mmseqs2_options) if bool(opts.mmseqs2_options) else "",
        "--diamond_options {}".format(opts.diamond_options) if bool(opts.diamond_options) else "",
        "--minimum_core_prevalence {}".format(opts.minimum_core_prevalence),

            "&&",

        "SRC={}; DST={}; SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST/global".format(
            os.path.join(output_directory, "output"), 
            directories["output"],
            ),
    ]

    return cmd

# Local clustering
def get_local_clustering_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command

    cmd = [
        os.environ["local_clustering.py"],
        "-i {}".format(input_filepaths[0]),
        "-o {}".format(output_directory),
        # "--no_singletons" if bool(opts.no_singletons) else "",
        "-p {}".format(opts.n_jobs),

        "--genome_clustering_algorithm {}".format(opts.genome_clustering_algorithm),
        "--ani_threshold {}".format(opts.ani_threshold),
        "--genome_cluster_prefix {}".format(opts.genome_cluster_prefix) if bool(opts.genome_cluster_prefix) else "",
        "--genome_cluster_suffix {}".format(opts.genome_cluster_suffix) if bool(opts.genome_cluster_suffix) else "",
        "--genome_cluster_prefix_zfill {}".format(opts.genome_cluster_prefix_zfill) if bool(opts.genome_cluster_prefix_zfill) else "",
        "--skani_target_ani {}".format(opts.skani_target_ani),
        "--skani_minimum_af {}".format(opts.skani_minimum_af),
        "--skani_no_confidence_interval" if opts.skani_no_confidence_interval else "",

        "--skani_nonviral_preset {}".format(opts.skani_nonviral_preset),
        "--skani_nonviral_compression_factor {}".format(opts.skani_nonviral_compression_factor),
        "--skani_nonviral_marker_kmer_compression_factor {}".format(opts.skani_nonviral_marker_kmer_compression_factor),
        "--skani_nonviral_options {}".format(opts.skani_nonviral_options) if bool(opts.skani_nonviral_options) else "",

        "--skani_viral_preset {}".format(opts.skani_viral_preset),
        "--skani_viral_compression_factor {}".format(opts.skani_viral_compression_factor),
        "--skani_viral_marker_kmer_compression_factor {}".format(opts.skani_viral_marker_kmer_compression_factor),
        "--skani_viral_options {}".format(opts.skani_viral_options) if bool(opts.skani_viral_options) else "",

        "--fastani_options {}".format(opts.fastani_options) if bool(opts.fastani_options) else "",

        "--protein_clustering_algorithm {}".format(opts.protein_clustering_algorithm),
        "--minimum_identity_threshold {}".format(opts.minimum_identity_threshold),
        "--minimum_coverage_threshold {}".format(opts.minimum_coverage_threshold),
        "--protein_cluster_prefix {}".format(opts.protein_cluster_prefix) if bool(opts.protein_cluster_prefix) else "",
        "--protein_cluster_suffix {}".format(opts.protein_cluster_suffix) if bool(opts.protein_cluster_suffix) else "",
        "--protein_cluster_prefix_zfill {}".format(opts.protein_cluster_prefix_zfill) if bool(opts.protein_cluster_prefix_zfill) else "",
        "--mmseqs2_options {}".format(opts.mmseqs2_options) if bool(opts.mmseqs2_options) else "",
        "--diamond_options {}".format(opts.diamond_options) if bool(opts.diamond_options) else "",
        "--minimum_core_prevalence {}".format(opts.minimum_core_prevalence),

            "&&",

        "SRC={}; DST={}; SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST/local".format(
            os.path.join(output_directory, "output"), 
            directories["output"],
            ),

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
                "skani",
                "fastANI",
                "mmseqs",
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


    # Genome clustering algorithm
    GENOME_CLUSTERING_ALGORITHM = opts.genome_clustering_algorithm.lower()
    if GENOME_CLUSTERING_ALGORITHM == "fastani":
        GENOME_CLUSTERING_ALGORITHM = "FastANI"
    if GENOME_CLUSTERING_ALGORITHM == "skani":
        GENOME_CLUSTERING_ALGORITHM = "skani"

    # Protein clustering algorithm
    PROTEIN_CLUSTERING_ALGORITHM = opts.protein_clustering_algorithm.split("-")[0].lower()
    if PROTEIN_CLUSTERING_ALGORITHM == "mmseqs":
        PROTEIN_CLUSTERING_ALGORITHM = PROTEIN_CLUSTERING_ALGORITHM.upper()
    if PROTEIN_CLUSTERING_ALGORITHM == "diamond":
        PROTEIN_CLUSTERING_ALGORITHM = PROTEIN_CLUSTERING_ALGORITHM.capitalize()

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
    description = "Global clustering of genomes ({}) and proteins ({})".format(GENOME_CLUSTERING_ALGORITHM, PROTEIN_CLUSTERING_ALGORITHM)

    # i/o
    input_filepaths = [opts.genomes_table]

    output_filenames = [
        "output/*.tsv.gz",
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
    if opts.local_clustering:
        step = 2

        program = "local_clustering"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

        # Info
        description = "Local clustering of genomes ({}) and proteins ({})".format(GENOME_CLUSTERING_ALGORITHM, PROTEIN_CLUSTERING_ALGORITHM)

        # i/o
        input_filepaths = [opts.genomes_table]

        output_filenames = [
            "output/*.tsv.gz",
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
    
    assert_acceptable_arguments(opts.protein_clustering_algorithm, {"easy-cluster", "easy-linclust", "mmseqs-cluster", "mmseqs-linclust", "diamond-cluster", "diamond-linclust"})
    if opts.protein_clustering_algorithm in {"easy-cluster", "easy-linclust"}:
        d = {"easy-cluster":"mmseqs-cluster", "easy-linclust":"mmseqs-linclust"}
        warnings.warn("\n\nPlease use `{}` instead of `{}` for MMSEQS2 clustering.".format(d[opts.protein_clustering_algorithm], opts.protein_clustering_algorithm))
        opts.protein_clustering_algorithm = d[opts.protein_clustering_algorithm]

    if opts.skani_nonviral_preset.lower() == "none":
        opts.skani_nonviral_preset = None

    if opts.skani_viral_preset.lower() == "none":
        opts.skani_viral_preset = None

    assert 0 < opts.minimum_core_prevalence <= 1.0, "--minimum_core_prevalence must be a float between (0.0,1.0])"
    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <genomes_table.tsv> -o <output_directory> -A 95 -a mmseqs-cluster".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i", "--genomes_table", type=str, default="stdin",  help = "path/to/genomes_table.tsv, Format: Must include the following columns (No header) [organism_type]<tab>[id_sample]<tab>[id_mag]<tab>[genome]<tab>[proteins]<tab>[cds] but can include additional columns to the right (e.g., [gene_models]).  Suggested input is from `compile_genomes_table.py` script. [Default: stdin]")
    parser_io.add_argument("-o","--output_directory", type=str, default="veba_output/cluster", help = "path/to/project_directory [Default: veba_output/cluster]")
    # parser_io.add_argument("-e", "--no_singletons", action="store_true", help="Exclude singletons") #isPSLC-1_SSO-3345__SRR178126
    parser_io.add_argument("-l", "--local_clustering", action="store_true", help = "Perform local clustering after global clustering")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # ANI
    parser_genome_clustering = parser.add_argument_group('Genome clustering arguments')
    parser_genome_clustering.add_argument("-G", "--genome_clustering_algorithm", type=str,  choices={"fastani", "skani"}, default="skani", help="Program to use for ANI calculations.  `skani` is faster and more memory efficient. For v1.0.0 - v1.3.x behavior, use `fastani`. [Default: skani]")
    parser_genome_clustering.add_argument("-A", "--ani_threshold", type=float, default=95.0, help="Species-level cluster (SLC) ANI threshold (Range (0.0, 100.0]) [Default: 95.0]")
    parser_genome_clustering.add_argument("--genome_cluster_prefix", type=str, default="SLC-", help="Cluster prefix [Default: 'SLC-")
    parser_genome_clustering.add_argument("--genome_cluster_suffix", type=str, default="", help="Cluster suffix [Default: '")
    parser_genome_clustering.add_argument("--genome_cluster_prefix_zfill", type=int, default=0, help="Cluster prefix zfill. Use 7 to match identifiers from OrthoFinder.  Use 0 to add no zfill. [Default: 0]") #7

    parser_skani = parser.add_argument_group('Skani triangle arguments')
    parser_skani.add_argument("--skani_target_ani",  type=float, default=80, help="skani | If you set --skani_target_ani to --ani_threshold, you may screen out genomes ANI ≥ --ani_threshold [Default: 80]")
    parser_skani.add_argument("--skani_minimum_af",  type=float, default=15, help="skani | Minimum aligned fraction greater than this value [Default: 15]")
    parser_skani.add_argument("--skani_no_confidence_interval",  action="store_true", help="skani | Output [5,95] ANI confidence intervals using percentile bootstrap on the putative ANI distribution")
    # parser_skani.add_argument("--skani_low_memory", action="store_true", help="Skani | More options (e.g. --arg 1 ) https://github.com/bluenote-1577/skani [Default: '']")

    parser_skani = parser.add_argument_group('[Prokaryotic & Eukaryotic] Skani triangle arguments')
    parser_skani.add_argument("--skani_nonviral_preset", type=str, default="medium", choices={"fast", "medium", "slow", "none"}, help="skani [Prokaryotic & Eukaryotic]| Use `none` if you are setting skani -c (compression factor) {fast, medium, slow, none} [Default: medium]")
    parser_skani.add_argument("--skani_nonviral_compression_factor", type=int, default=125,  help="skani [Prokaryotic & Eukaryotic]|  Compression factor (k-mer subsampling rate).	[Default: 125]")
    parser_skani.add_argument("--skani_nonviral_marker_kmer_compression_factor", type=int, default=1000,  help="skani [Prokaryotic & Eukaryotic] | Marker k-mer compression factor. Markers are used for filtering. [Default: 1000]")
    parser_skani.add_argument("--skani_nonviral_options", type=str, default="", help="skani [Prokaryotic & Eukaryotic] | More options for `skani triangle` (e.g. --arg 1 ) [Default: '']")

    parser_skani = parser.add_argument_group('[Viral] Skani triangle arguments')
    parser_skani.add_argument("--skani_viral_preset", type=str, default="slow", choices={"fast", "medium", "slow", "none"}, help="skani | Use `none` if you are setting skani -c (compression factor) {fast, medium, slow, none} [Default: slow]")
    parser_skani.add_argument("--skani_viral_compression_factor", type=int, default=30,  help="skani [Viral] | Compression factor (k-mer subsampling rate).	[Default: 30]")
    parser_skani.add_argument("--skani_viral_marker_kmer_compression_factor", type=int, default=200,  help="skani [Viral] | Marker k-mer compression factor. Markers are used for filtering. Consider decreasing to ~200-300 if working with small genomes (e.g. plasmids or viruses). [Default: 200]")
    parser_skani.add_argument("--skani_viral_options", type=str, default="", help="skani [Viral] | More options for `skani triangle` (e.g. --arg 1 ) [Default: '']")

    parser_fastani = parser.add_argument_group('FastANI arguments')
    parser_fastani.add_argument("--fastani_options", type=str, default="", help="FastANI | More options (e.g. --arg 1 ) [Default: '']")

    # Clustering
    parser_protein_clustering = parser.add_argument_group('Protein clustering arguments')
    parser_protein_clustering.add_argument("-P", "--protein_clustering_algorithm", type=str, choices={"mmseqs-cluster", "mmseqs-linclust", "diamond-cluster", "diamond-linclust"}, default="mmseqs-cluster", help="Clustering algorithm | Diamond can only be used for clustering proteins {mmseqs-cluster, mmseqs-linclust, diamond-cluster, mmseqs-linclust} [Default: mmseqs-cluster]")
    parser_protein_clustering.add_argument("-t", "--minimum_identity_threshold", type=float, default=50.0, help="Clustering | Percent identity threshold (Range (0.0, 100.0]) [Default: 50.0]")
    parser_protein_clustering.add_argument("-c", "--minimum_coverage_threshold", type=float, default=0.8, help="Clustering | Coverage threshold (Range (0.0, 1.0]) [Default: 0.8]")
    parser_protein_clustering.add_argument("--protein_cluster_prefix", type=str, default="SSPC-", help="Cluster prefix [Default: 'SSPC-")
    parser_protein_clustering.add_argument("--protein_cluster_suffix", type=str, default="", help="Cluster suffix [Default: '")
    parser_protein_clustering.add_argument("--protein_cluster_prefix_zfill", type=int, default=0, help="Cluster prefix zfill. Use 7 to match identifiers from OrthoFinder.  Use 0 to add no zfill. [Default: 0]") #7
    parser_protein_clustering.add_argument("--mmseqs2_options", type=str, default="", help="MMSEQS2 | More options (e.g. --arg 1 ) [Default: '']")
    parser_protein_clustering.add_argument("--diamond_options", type=str, default="", help="Diamond | More options (e.g. --arg 1 ) [Default: '']")

    # Pangenome
    parser_pangenome = parser.add_argument_group('Pangenome arguments')
    parser_pangenome.add_argument("--minimum_core_prevalence", type=float, default=1.0, help="Minimum ratio of genomes detected in a SLC for a SSPC to be considered core (Range (0.0, 1.0]) [Default: 1.0]")

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
