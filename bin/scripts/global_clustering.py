#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, shutil, time, gzip, warnings
from multiprocessing import cpu_count
from collections import OrderedDict, defaultdict

import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser 

# Soothsayer Ecosystem
from genopype import *
from genopype import __version__ as genopype_version
from soothsayer_utils import *

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.3.26"

def get_basename(x):
    _, fn = os.path.split(x)
    if fn.endswith(".gz"):
        fn = fn[:-3]
    return ".".join(fn.split(".")[:-1])

def get_protein_cluster_prevalence(df_input:pd.DataFrame):
    # Read Input
    genomes = sorted(df_input.iloc[:,0].unique())
    clusters = sorted(df_input.iloc[:,2].unique())

    # Create array
    A = np.zeros((len(genomes), len(clusters)), dtype=int)

    for _, (id_genome, id_protein, id_cluster) in df_input.iterrows():
        i = genomes.index(id_genome)
        j = clusters.index(id_cluster)
        A[i,j] += 1

    # Create output
    df_output = pd.DataFrame(A, index=genomes, columns=clusters)
    df_output.index.name = "id_genome"
    df_output.columns.name = "id_protein_cluster"

    return df_output

# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = {
        "edgelist_to_clusters.py",
        "clustering_wrapper.py",
        # "table_to_fasta.py",
    }

    required_executables={
        "skani",
        "fastANI",
        "mmseqs",
        "diamond",
     } 
  

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
    usage = "{} -i <genomes_table.tsv> -o <output_directory> -A 95 -a easy-cluster".format(__program__)

    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i", "--genomes_table", type=str, default="stdin",  help = "path/to/genomes_table.tsv, Format: Must include the following columns (No header) [organism_type]<tab>[id_sample]<tab>[id_mag]<tab>[genome]<tab>[proteins]<tab>[cds] but can include additional columns to the right (e.g., [gene_models]).  Suggested input is from `compile_genomes_table.py` script. [Default: stdin]")
    parser_io.add_argument("-o","--output_directory", type=str, default="global_clustering_output", help = "path/to/project_directory [Default: global_clustering_output]")
    parser_io.add_argument("-e", "--no_singletons", action="store_true", help="Exclude singletons")
    parser_io.add_argument("-R", "--no_representative_sequences", action="store_true", help="Do not write representative sequences to fasta") 
    parser_io.add_argument("-C", "--no_core_sequences", action="store_true", help="Do not write core pagenome sequences to fasta") 
    # parser_io.add_argument("-M", "--no_marker_sequences", action="store_true", help="Do not write core pagenome sequences to fasta") 

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
    # parser_utility.add_argument("--verbose", action='store_true')

    # ANI
    parser_genome_clustering = parser.add_argument_group('Genome clustering arguments')
    parser_genome_clustering.add_argument("-G", "--genome_clustering_algorithm", type=str,  choices={"fastani", "skani"}, default="skani", help="Program to use for ANI calculations.  `skani` is faster and more memory efficient. For v1.0.0 - v1.3.x behavior, use `fastani`. [Default: skani]")
    parser_genome_clustering.add_argument("-A", "--ani_threshold", type=float, default=95.0, help="Species-level cluster (SLC) ANI threshold (Range (0.0, 100.0]) [Default: 95.0]")
    parser_genome_clustering.add_argument("-F", "--af_threshold", type=float, default=30.0, help="Species-level cluster (SLC) alignment fraction threshold. Only available if `skani` is used as --genome_clustering_algorithm. (Range (0.0, 100.0]) [Default: 30.0]")

    parser_genome_clustering.add_argument("--genome_cluster_prefix", type=str, default="SLC-", help="Cluster prefix [Default: 'SLC-")
    parser_genome_clustering.add_argument("--genome_cluster_suffix", type=str, default="", help="Cluster suffix [Default: '")
    parser_genome_clustering.add_argument("--genome_cluster_prefix_zfill", type=int, default=0, help="Cluster prefix zfill. Use 7 to match identifiers from OrthoFinder.  Use 0 to add no zfill. [Default: 0]") #7

    parser_skani = parser.add_argument_group('Skani triangle arguments')
    parser_skani.add_argument("--skani_target_ani",  type=float, default=80, help="skani | If you set --skani_target_ani to --ani_threshold, you may screen out genomes ANI ≥ --ani_threshold [Default: 80]")
    parser_skani.add_argument("--skani_minimum_af",  type=float, default=15, help="skani | Minimum aligned fraction greater than this value. If you set --skani_minimum_af to --af_threshold, you may screen out genomes AF ≥ --af_threshold [Default: 15]")
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
        opts.n_jobs = cpu_count()
    assert opts.n_jobs >= 1, "--n_jobs must be ≥ 1 (or -1 to use all available threads)"

    # Directories
    directories = dict()
    directories["project"] = create_directory(opts.output_directory)
    directories["output"] = create_directory(os.path.join(directories["project"], "output"))
    directories["pangenome_tables"] = create_directory(os.path.join(directories["output"], "pangenome_tables"))
    if not opts.no_representative_sequences:
        directories["pangenome_core_sequences"] = create_directory(os.path.join(directories["output"], "pangenome_core_sequences"))
    directories["serialization"] = create_directory(os.path.join(directories["output"], "serialization"))
    directories["representatives"] = create_directory(os.path.join(directories["output"], "representatives"))
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
    print("Commands:", list(filter(bool,sys.argv)),  sep="\n", file=sys.stdout)
    configure_parameters(opts, directories)
    sys.stdout.flush()

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

    # Make directories
    t0 = time.time()
    print(format_header(" " .join(["* ({}) Creating directories:".format(format_duration(t0)), directories["intermediate"]])), file=sys.stdout)
    os.makedirs(opts.output_directory, exist_ok=True)
    
    # Load input
    if opts.genomes_table == "stdin":
        opts.genomes_table = sys.stdin
    df_genomes = pd.read_csv(opts.genomes_table, sep="\t", header=None)
    assert df_genomes.shape[1] >= 6, "Must include the follow columns (No header) [organism_type]<tab>[id_sample]<tab>[id_mag]<tab>[genome]<tab>[proteins]<tab>[cds] but can include additional columns to the right (e.g., [gene_models]).  Suggested input is from `compile_genomes_table.py` script.  You have provided a table with {} columns.".format(df_genomes.shape[1])
    df_genomes = df_genomes.iloc[:,:6]
    df_genomes.columns = ["organism_type", "id_sample", "id_mag", "genome", "proteins", "cds"]
    assert not np.any(df_genomes.isnull()), "--genomes_table has missing values.  Please correct this.\n{}".format(df_genomes.loc[df_genomes.isnull().sum(axis=1)[lambda x: x > 0].index].to_string())


    for organism_type, df_1 in df_genomes.groupby("organism_type"):
        # Organism directories
        os.makedirs(os.path.join(directories["intermediate"], organism_type), exist_ok=True)
        os.makedirs(os.path.join(directories["intermediate"], organism_type, "clusters"), exist_ok=True)

        # Write genomes to file
        with open(os.path.join(directories["intermediate"], organism_type, "genomes.list"), "w") as f:
            print(*sorted(df_1["genome"]), sep="\n", file=f)

        # Write genome identifiers to file
        with open(os.path.join(directories["intermediate"], organism_type, "genome_identifiers.list"), "w") as f:
            print(*sorted(df_1["genome"].map(get_basename)), sep="\n", file=f)


    # Get identifier mapping
    mag_to_numberofscaffolds = defaultdict(int)
    mag_to_numberofproteins = defaultdict(int)
    mag_to_sample = dict()
    scaffold_to_mag = dict()
    scaffold_to_sample = dict()
    protein_to_orthogroup = dict()
    protein_to_mag = dict()
    protein_to_sample = dict()
    protein_to_sequence = dict()
    protein_to_cds = dict()

    for i, row in pv(df_genomes.iterrows(), "Organizing identifiers"):
        id_mag = row["id_mag"]
        id_sample = row["id_sample"]
        mag_to_sample[id_mag] = id_sample
        with get_file_object(row["genome"], "read", verbose=False) as f:
            for header, seq in SimpleFastaParser(f):
                id_scaffold = header.split(" ")[0]
                assert id_scaffold not in scaffold_to_mag, "Duplicate scaffold/contig identifiers are not allowed: {} from {}".format(id_scaffold, row["genome"])
                scaffold_to_mag[id_scaffold] = id_mag 
                scaffold_to_sample[id_scaffold] = id_sample 
                mag_to_numberofscaffolds[id_mag] += 1

        with get_file_object(row["proteins"], "read", verbose=False) as f:
            for header, seq in SimpleFastaParser(f):
                id_protein = header.split(" ")[0]
                assert id_protein not in protein_to_orthogroup, "Duplicate ORF/protein identifiers are not allowed: {} from {}".format(id_protein, row["proteins"])
                protein_to_mag[id_protein] = id_mag
                protein_to_sample[id_protein] = id_sample 
                protein_to_sequence[id_protein] = seq
                mag_to_numberofproteins[id_mag] += 1

        with get_file_object(row["cds"], "read", verbose=False) as f:
            for header, seq in SimpleFastaParser(f):
                id_protein = header.split(" ")[0]
                assert id_protein in protein_to_sequence, "CDS sequence identifier must be in protein fasta: {} from {}".format(id_protein, row["cds"])
                protein_to_cds[id_protein] = seq

    mag_to_numberofscaffolds = pd.Series(mag_to_numberofscaffolds)
    mag_to_numberofproteins = pd.Series(mag_to_numberofproteins)
    mag_to_sample = pd.Series(mag_to_sample)
    mag_to_organismtype = pd.Series(dict(zip(df_genomes["id_mag"], df_genomes["organism_type"])))
    scaffold_to_mag = pd.Series(scaffold_to_mag)
    scaffold_to_sample = pd.Series(scaffold_to_sample)
    protein_to_mag = pd.Series(protein_to_mag)
    protein_to_sample = pd.Series(protein_to_sample)
    protein_to_sequence = pd.Series(protein_to_sequence)
    protein_to_cds = pd.Series(protein_to_cds)

    # Commands
    f_cmds = open(os.path.join(opts.output_directory, "commands.sh"), "w")

    # Pairwise ANI
    print(format_header("* ({}) Running {}:".format(format_duration(t0), GENOME_CLUSTERING_ALGORITHM)), file=sys.stdout)
    for fp in pv(glob.glob(os.path.join(directories["intermediate"], "*",  "genomes.list")), "Running pairwise ANI"):
        fields = fp.split("/")
        organism_type = fields[-2]
        output_directory = os.path.split(fp)[0]

        if opts.genome_clustering_algorithm == "skani":
            name = "skani__{}".format(organism_type)
            description = "[Program = skani] [Organism_Type = {}]".format(organism_type)

            arguments = list()

            if organism_type.lower() in {"viral", "virus", "virion"}:
                arguments += [
                os.environ["skani"],
                "triangle",
                "--sparse",
                "-t {}".format(opts.n_jobs),
                "-l {}".format(fp),
                "-o {}".format(os.path.join(output_directory, "skani_output.tsv")),
                "--ci" if not opts.skani_no_confidence_interval else "",
                "--min-af {}".format(opts.skani_minimum_af),
                "-s {}".format(opts.skani_target_ani),
                "-c {}".format(opts.skani_viral_compression_factor),
                "-m {}".format(opts.skani_viral_marker_kmer_compression_factor),
                "--{}".format(opts.skani_viral_preset) if opts.skani_viral_preset else "",
                opts.skani_viral_options,
                ]

            else:
                arguments += [
                os.environ["skani"],
                "triangle",
                "--sparse",
                "-t {}".format(opts.n_jobs),
                "-l {}".format(fp),
                "-o {}".format(os.path.join(output_directory, "skani_output.tsv")),
                "--ci" if not opts.skani_no_confidence_interval else "",
                "--min-af {}".format(opts.skani_minimum_af),
                "-s {}".format(opts.skani_target_ani),
                "-c {}".format(opts.skani_nonviral_compression_factor),
                "-m {}".format(opts.skani_nonviral_marker_kmer_compression_factor),
                "--{}".format(opts.skani_nonviral_preset) if opts.skani_nonviral_preset else "",
                opts.skani_nonviral_options,
                ]

            arguments += [
                    "&&",

                "cat",
                os.path.join(output_directory, "skani_output.tsv"),
                "|",
                "cut -f1-4",
                "|",
                "tail -n +2",
                "|",
                os.environ["edgelist_to_clusters.py"],
                "--basename",
                "-t {}".format(opts.ani_threshold),
                "-a {}".format(opts.af_threshold),

                "--no_singletons" if bool(opts.no_singletons) else "",
                "--cluster_prefix {}{}".format(organism_type[0].upper(), opts.genome_cluster_prefix),
                "--cluster_suffix {}".format(opts.genome_cluster_suffix) if bool(opts.genome_cluster_suffix) else "",
                "--cluster_prefix_zfill {}".format(opts.genome_cluster_prefix_zfill),
                "-o {}".format(os.path.join(output_directory, "genome_clusters.tsv")),
                "--identifiers {}".format(os.path.join(directories["intermediate"], organism_type, "genome_identifiers.list")),
                "--export_graph {}".format(os.path.join(directories["serialization"], f"{organism_type}.networkx_graph.pkl")),
                "--export_dict {}".format(os.path.join(directories["serialization"], f"{organism_type}.dict.pkl")),
                "--export_representatives {}".format(os.path.join(directories["representatives"], f"{organism_type}.representatives.tsv")),

                    "&&",

                "rm -rf {}".format(os.path.join(directories["tmp"], "*")),

                ] 
            
            cmd = Command(
                arguments,
                name=name, 
                f_cmds=f_cmds,
                )
        
        if opts.genome_clustering_algorithm == "fastani":
            name = "fastani__{}".format(organism_type)
            description = "[Program = FastANI] [Organism_Type = {}]".format(organism_type)
            cmd = Command([
                os.environ["fastANI"],
                "-t {}".format(opts.n_jobs),
                "--rl {}".format(fp),
                "--ql {}".format(fp),
                "-o {}".format(os.path.join(output_directory, "fastani_output.tsv")),
                opts.fastani_options,

                    "&&",

                "cat",
                os.path.join(output_directory, "fastani_output.tsv"),
                "|",
                "cut -f1-3",
                "|",
                os.environ["edgelist_to_clusters.py"],
                "--basename",
                "-t {}".format(opts.ani_threshold),
                "--no_singletons" if bool(opts.no_singletons) else "",
                "--cluster_prefix {}{}".format(organism_type[0].upper(), opts.genome_cluster_prefix),
                "--cluster_suffix {}".format(opts.genome_cluster_suffix) if bool(opts.genome_cluster_suffix) else "",
                "--cluster_prefix_zfill {}".format(opts.genome_cluster_prefix_zfill),
                "-o {}".format(os.path.join(output_directory, "genome_clusters.tsv")),
                "--identifiers {}".format(os.path.join(directories["intermediate"], organism_type, "genome_identifiers.list")),
                "--export_graph {}".format(os.path.join(directories["serialization"], f"{organism_type}.networkx_graph.pkl")),
                "--export_dict {}".format(os.path.join(directories["serialization"], f"{organism_type}.dict.pkl")),
                "--export_representatives {}".format(os.path.join(directories["representatives"], f"{organism_type}.representatives.tsv")),

                    "&&",

                "rm -rf {}".format(os.path.join(directories["tmp"], "*")),

                ], 
                name=name, 
                f_cmds=f_cmds,
                )

        # Run command
        cmd.run(
            checkpoint_message_notexists="[Running ({})] | {}".format(format_duration(t0), description),
            checkpoint_message_exists="[Loading Checkpoint ({})] | {}".format(format_duration(t0), description),
            write_stdout=os.path.join(directories["log"], "{}.o".format(name)),
            write_stderr=os.path.join(directories["log"], "{}.e".format(name)),
            write_returncode=os.path.join(directories["log"], "{}.returncode".format(name)),
            checkpoint=os.path.join(directories["checkpoints"], name),
            )
        if hasattr(cmd, "returncode_"):
            if cmd.returncode_ != 0:
                print("[Error] | {}".format(description), file=sys.stdout)
                print("Check the following files:\ncat {}".format(os.path.join(directories["log"], "{}.*".format(name))), file=sys.stdout)
                sys.exit(cmd.returncode_)

    # Protein Clustering
    print(format_header(" * ({}) Running {}:".format(format_duration(t0), PROTEIN_CLUSTERING_ALGORITHM)), file=sys.stdout)
    mag_to_genomecluster = dict()
    protein_to_proteincluster = dict()
    for fp in pv(glob.glob(os.path.join(directories["intermediate"], "*", "genome_clusters.tsv")), "Running {}".format(PROTEIN_CLUSTERING_ALGORITHM)):
        fields = fp.split("/")
        organism_type = fields[-2]

        mag_to_genomecluster_within_organism_type = pd.read_csv(fp, sep="\t", index_col=0, header=None).iloc[:,0]
        mag_to_genomecluster.update(mag_to_genomecluster_within_organism_type)

        for id_genomecluster, data in mag_to_genomecluster_within_organism_type.groupby(mag_to_genomecluster_within_organism_type):
            genomecluster_directory = os.path.join(directories["intermediate"], organism_type, "clusters", id_genomecluster)
            os.makedirs(genomecluster_directory, exist_ok=True)

            # Get MAGs and proteins
            mags = data.index
            proteins = protein_to_mag[protein_to_mag.map(lambda x: x in mags)].index 
            with open(os.path.join(genomecluster_directory, "genomes.list" ), "w") as f:
                print(*sorted(mags), sep="\n", file=f)
            with open(os.path.join(genomecluster_directory, "protein_identifiers.list"), "w") as f:
                print(*proteins, sep="\n", file=f)
            write_fasta(protein_to_sequence[proteins], os.path.join(genomecluster_directory, "proteins.faa" ))

            # Run Clustering
            name = "{}__{}__{}".format(PROTEIN_CLUSTERING_ALGORITHM.lower(), organism_type, id_genomecluster)
            description = "[Program = {}] [Organism_Type = {}] [Genome_Cluster = {}]".format(PROTEIN_CLUSTERING_ALGORITHM, organism_type, id_genomecluster)

            cmd = Command([
                os.environ["clustering_wrapper.py"],
                "--fasta {}".format(os.path.join(genomecluster_directory, "proteins.faa" )),
                "--output_directory {}".format(genomecluster_directory),
                "--no_singletons" if bool(opts.no_singletons) else "",
                "--algorithm {}".format(opts.protein_clustering_algorithm),
                "--n_jobs {}".format(opts.n_jobs),
                "--minimum_identity_threshold {}".format(opts.minimum_identity_threshold),
                "--minimum_coverage_threshold {}".format(opts.minimum_coverage_threshold),
                "--mmseqs2_options='{}'" if bool(opts.mmseqs2_options) else "",
                "--diamond_options='{}'" if bool(opts.diamond_options) else "",
                "--cluster_prefix {}_{}".format(id_genomecluster, opts.protein_cluster_prefix),
                "--cluster_suffix {}".format(opts.protein_cluster_suffix) if bool(opts.protein_cluster_suffix) else "",
                "--cluster_prefix_zfill {}".format(opts.protein_cluster_prefix_zfill),
                "--basename protein_clusters",
                "--identifiers {}".format(os.path.join(genomecluster_directory, "protein_identifiers.list")),
                "--no_sequences_and_header",
                "--representative_output_format table",

                    "&&",

                "rm -rf",
                os.path.join(genomecluster_directory, "proteins.faa" ),
                ], 
                name=name, 
                f_cmds=f_cmds,
                )

            # Run command
            cmd.run(
                checkpoint_message_notexists="[Running ({})] | {}".format(format_duration(t0), description),
                checkpoint_message_exists="[Loading Checkpoint ({})] | {}".format(format_duration(t0), description),
                write_stdout=os.path.join(directories["log"], "{}.o".format(name)),
                write_stderr=os.path.join(directories["log"], "{}.e".format(name)),
                write_returncode=os.path.join(directories["log"], "{}.returncode".format(name)),
                checkpoint=os.path.join(directories["checkpoints"], name),
                )
            if hasattr(cmd, "returncode_"):
                if cmd.returncode_ != 0:
                    print("[Error] | {}".format(description), file=sys.stdout)
                    print("Check the following files:\ncat {}".format(os.path.join(directories["log"], "{}.*".format(name))), file=sys.stdout)
                    sys.exit(cmd.returncode_)

            # Get the protein clusters
            protein_to_proteinluster_within_organism_type = pd.read_csv(os.path.join(genomecluster_directory, "output", "protein_clusters.tsv"), sep="\t", index_col=0, header=None).iloc[:,0]
            protein_to_proteincluster.update(protein_to_proteinluster_within_organism_type.to_dict())
    mag_to_genomecluster = pd.Series(mag_to_genomecluster)
    protein_to_proteincluster = pd.Series(protein_to_proteincluster)

    f_cmds.close()
    
    # Format output tables 
    print(format_header(" * ({}) Formatting Output Tables:".format(format_duration(t0))), file=sys.stdout)
    # MAGs
    print(" * ({}) Compiling genome identifier mappings".format(format_duration(t0)), file=sys.stdout)

    df_mags = pd.DataFrame(
        OrderedDict([
            ("organism_type", mag_to_organismtype),
            ("sample_of_origin", mag_to_sample),
            ("id_genome_cluster", mag_to_genomecluster),
            ("number_of_proteins", mag_to_numberofproteins),
        ])
    ).sort_values(["sample_of_origin", "id_genome_cluster"])
    df_mags.index.name = "id_genome"

    # Scaffolds
    print(" * ({}) Compiling scaffold identifier mappings".format(format_duration(t0)), file=sys.stdout)
    df_scaffolds = pd.DataFrame(
        OrderedDict([
            ("organism_type", scaffold_to_mag.map(lambda id_mag: mag_to_organismtype[id_mag])),
            ("id_genome", scaffold_to_mag),
            ("sample_of_origin", scaffold_to_sample),
            ("id_genome_cluster", scaffold_to_mag.map(lambda id_mag: mag_to_genomecluster[id_mag])),
        ])
    ).sort_values(["sample_of_origin", "id_genome", "id_genome_cluster"])
    df_scaffolds.index.name = "id_scaffold"

    # Proteins
    print(" * ({}) Compiling protein identifier mappings".format(format_duration(t0)), file=sys.stdout)
    df_proteins = pd.DataFrame(
        OrderedDict([
            ("organism_type", protein_to_mag.map(lambda id_mag: mag_to_organismtype[id_mag])),
            ("id_genome", protein_to_mag),
            ("sample_of_origin", protein_to_sample),
            ("id_genome_cluster", protein_to_mag.map(lambda id_mag: mag_to_genomecluster[id_mag])),
            ("id_protein_cluster", protein_to_proteincluster),
        ])
    ).sort_values(["sample_of_origin", "id_genome", "id_genome_cluster", "id_protein_cluster"])
    df_proteins.index.name = "id_protein"

    print(" * ({}) Compiling genome cluster identifier mappings".format(format_duration(t0)), file=sys.stdout)
    genomecluster_data = defaultdict(dict) 
    for id_genomecluster, df in df_mags.groupby("id_genome_cluster"):
        genomecluster_data[id_genomecluster]["number_of_components"] = df["id_genome_cluster"].dropna().size
        genomecluster_data[id_genomecluster]["components"] = set(df["id_genome_cluster"].dropna().index)
        genomecluster_data[id_genomecluster]["number_of_samples_of_origin"] = df["sample_of_origin"].dropna().nunique()
        genomecluster_data[id_genomecluster]["samples_of_origin"] = set(df["sample_of_origin"].dropna().unique())
    df_genomeclusters = pd.DataFrame(genomecluster_data).T.sort_index()[["number_of_components", "number_of_samples_of_origin", "components", "samples_of_origin"]]
    df_genomeclusters.index.name = "id_genome_cluster"

    print(" * ({}) Compiling protein cluster identifier mappings".format(format_duration(t0)), file=sys.stdout)
    proteincluster_data = defaultdict(dict) 
    for id_proteincluster, df in df_proteins.groupby("id_protein_cluster"):
        proteincluster_data[id_proteincluster]["number_of_components"] = df["id_protein_cluster"].dropna().size
        proteincluster_data[id_proteincluster]["components"] = set(df["id_protein_cluster"].dropna().index)
        proteincluster_data[id_proteincluster]["number_of_samples_of_origin"] = df["sample_of_origin"].dropna().nunique()
        proteincluster_data[id_proteincluster]["samples_of_origin"] = set(df["sample_of_origin"].dropna().unique())
    df_proteinclusters = pd.DataFrame(proteincluster_data).T.sort_index()[["number_of_components", "number_of_samples_of_origin", "components", "samples_of_origin"]]
    df_proteinclusters.index.name = "id_protein_cluster"

    print(" * ({}) Calculating feature compression ratios for each sample".format(format_duration(t0)), file=sys.stdout)
    fcr_data = defaultdict(dict)
    for organism_type, df in df_mags.groupby("organism_type"):
        feature_to_cluster = df["id_genome_cluster"].dropna()
        fcr_data[organism_type]["number_of_genomes"] = feature_to_cluster.size
        fcr_data[organism_type]["number_of_genome-clusters"] = feature_to_cluster.nunique()
        fcr_data[organism_type]["genomic_fcr"] = 1 - (feature_to_cluster.nunique()/feature_to_cluster.size)
    for organism_type, df in df_proteins.groupby("organism_type"):
        feature_to_cluster = df["id_protein_cluster"].dropna()
        fcr_data[organism_type]["number_of_proteins"] = feature_to_cluster.size
        fcr_data[organism_type]["number_of_protein-clusters"] = feature_to_cluster.nunique()
        fcr_data[organism_type]["functional_fcr"] = 1 - (feature_to_cluster.nunique()/feature_to_cluster.size)
    df_fcr = pd.DataFrame(fcr_data).T.sort_index()
    df_fcr = df_fcr.loc[:,["number_of_genomes", "number_of_genome-clusters", "genomic_fcr", "number_of_proteins", "number_of_protein-clusters", "functional_fcr"]]
    for id_field in ["number_of_genomes", "number_of_genome-clusters", "number_of_proteins", "number_of_protein-clusters"]:
        df_fcr[id_field] = df_fcr[id_field].astype(int)
    df_fcr.index.name = "organism_type"

    # Get representative sequences
    if not opts.no_representative_sequences:
        # print(format_header(" * ({}) Compiling representative sequences:".format(format_duration(t0))), file=sys.stdout)

        proteincluster_to_representative = list()
        for fp in glob.glob(os.path.join(directories["intermediate"], "*", "clusters", "*","output", "representative_sequences.tsv.gz")):
            a = pd.read_csv(fp, sep="\t", index_col=0, header=None).iloc[:,0]
            proteincluster_to_representative.append(a)
        proteincluster_to_representative = pd.concat(proteincluster_to_representative)
        representatives = set(proteincluster_to_representative.values)

        df_proteins["protein_cluster_representative"] = df_proteins.index.map(lambda x: x in representatives)

        with open(os.path.join(directories["output"], "representative_sequences.faa"), "w") as f_representatives:
            for id_proteincluster, id_representative in pv(proteincluster_to_representative.items(), description=" * ({}) Writing protein cluster representative sequences".format(format_duration(t0)), unit="sequence", total=proteincluster_to_representative.size):
                seq = protein_to_sequence[id_representative]
                header = f"{id_proteincluster} {id_representative}"
                print(f">{header}\n{seq}", file=f_representatives)
    
    # Write prevalence table
    # print(format_header(" * ({}) Compiling pangenome protein cluster prevalence tables:".format(format_duration(t0))), file=sys.stdout)
    genome_to_number_of_singletons = list()
    genome_to_ratio_of_singletons = list()
    protein_to_number_of_genomes_detected = dict()
    protein_to_ratio_of_genomes_detected = dict()
    genomecluster_to_corepangenome = dict()
    genomecluster_to_singletons = dict()

    for id_genomecluster, df in pv(df_proteins.groupby("id_genome_cluster"), description=" * ({}) Writing protein prevalence pangenome tables, counting singletons, and core pangenome sequences".format(format_duration(t0)), unit="genome cluster", total=df_proteins["id_genome_cluster"].nunique()):
        # Prevalence
        df = df.reset_index().loc[:,["id_genome", "id_protein", "id_protein_cluster"]]
        df_prevalence = get_protein_cluster_prevalence(df)
        df_prevalence.to_csv(os.path.join(directories["pangenome_tables"], f"{id_genomecluster}.tsv.gz"), sep="\t")

        # Detected/Not-detected
        df_prevalence = df_prevalence > 0
        number_of_genomes_detected = df_prevalence.sum(axis=0)
        ratio_of_genomes_detected = df_prevalence.mean(axis=0)
        protein_to_number_of_genomes_detected.update(number_of_genomes_detected.to_dict())
        protein_to_ratio_of_genomes_detected.update(ratio_of_genomes_detected.to_dict())

        # Core
        core_proteinclusters = ratio_of_genomes_detected[lambda x: x == 1].index
        genomecluster_to_corepangenome[id_genomecluster] = set(core_proteinclusters)
        if not opts.no_core_sequences:
            with open(os.path.join(directories["pangenome_core_sequences"], f"{id_genomecluster}.faa"), "w") as f_core:
                for id_proteincluster, id_representative in proteincluster_to_representative[core_proteinclusters].items():
                    seq = protein_to_sequence[id_representative]
                    header = f"{id_proteincluster} {id_representative}"
                    print(f">{header}\n{seq}", file=f_core)

            with open(os.path.join(directories["pangenome_core_sequences"], f"{id_genomecluster}.ffn"), "w") as f_core:
                for id_proteincluster, id_representative in proteincluster_to_representative[core_proteinclusters].items():
                    seq = protein_to_cds[id_representative]
                    header = f"{id_proteincluster} {id_representative}"
                    print(f">{header}\n{seq}", file=f_core)


        # Singletons
        genomecluster_to_singletons[id_genomecluster] = set(number_of_genomes_detected[lambda x: x == 1].index)

        if df_prevalence.shape[0] > 1:  
            number_of_singletons = df_prevalence.loc[:,df_prevalence.sum(axis=0)[lambda x: x == 1].index].sum(axis=1)
            ratio_of_singletons = number_of_singletons/df_prevalence.sum(axis=1)
            genome_to_number_of_singletons.append(number_of_singletons)
            genome_to_ratio_of_singletons.append(ratio_of_singletons)

    if len(genome_to_number_of_singletons):
        df_mags["number_of_singleton_protein_clusters"] = pd.concat(genome_to_number_of_singletons).reindex(df_mags.index).astype("Int64")
        df_mags["ratio_of_protein_cluster_are_singletons"] = pd.concat(genome_to_ratio_of_singletons)
    else:
        warnings.warn("No clusters with more than 1 genome so singleton analysis does not apply")

    genomecluster_to_corepangenome = pd.Series(genomecluster_to_corepangenome)
    genomecluster_to_singletons = pd.Series(genomecluster_to_singletons)


    # Add detection number and ratios
    df_proteinclusters["number_of_genomes_detected"] = pd.Series(protein_to_number_of_genomes_detected)
    df_proteinclusters["ratio_of_genomes_detected"] = pd.Series(protein_to_ratio_of_genomes_detected)
    df_proteinclusters["core_pangenome"] = df_proteinclusters["ratio_of_genomes_detected"] >= opts.minimum_core_prevalence
    df_proteinclusters["singleton"] = df_proteinclusters["number_of_genomes_detected"] == 1


    df_genomeclusters["number_of_proteins_in_core_pangenome"] = genomecluster_to_corepangenome.map(len)
    df_genomeclusters["core_pangenome"] = genomecluster_to_corepangenome
    df_genomeclusters["number_of_singleton_proteins"] = genomecluster_to_singletons.map(len)
    df_genomeclusters["singletons"] = genomecluster_to_singletons

    # Number of copies of SSPC per genome
    number_of_gene_copies_per_sspc_per_genome  = df_proteins.groupby(["id_genome", "id_protein_cluster"]).size()
    df_proteinclusters["average_number_of_copies_per_genome"] = number_of_gene_copies_per_sspc_per_genome.groupby(lambda x: x[1]).mean()


    # # Symlink pangenome graphs
    # print(format_header(" * ({}) Symlinking pangenome protein cluster NetworkX Graphs:".format(format_duration(t0))), file=sys.stdout)
    # for src in glob.glob(os.path.join(directories["intermediate"], "*", "clusters", "*","output", "protein_clusters.networkx_graph.pkl")):
    #     id_genomecluster = src.split("/")[-3]
    #     dst = os.path.join(directories[("misc", "pangenomes", "networkx_graphs")], f"{id_genomecluster}.pkl")
    #     if os.path.exists(dst):
    #         os.remove(dst)
    #     src_real_path = os.path.realpath(src)
    #     src_relative_path = os.path.relpath(src_real_path, os.path.split(dst)[0])
    #     os.symlink(src_relative_path, dst)
          
    # Writing output files
    print(format_header(" * ({}) Writing Output Tables:".format(format_duration(t0))), file=sys.stdout)
    df_mags.to_csv(os.path.join(directories["output"], "identifier_mapping.genomes.tsv.gz"), sep="\t")
    df_scaffolds.to_csv(os.path.join(directories["output"], "identifier_mapping.scaffolds.tsv.gz"), sep="\t")
    df_proteins.to_csv(os.path.join(directories["output"], "identifier_mapping.proteins.tsv.gz"), sep="\t")
    df_genomeclusters.to_csv(os.path.join(directories["output"], "genome_clusters.tsv.gz"), sep="\t")
    df_proteinclusters.to_csv(os.path.join(directories["output"], "protein_clusters.tsv.gz"), sep="\t")
    df_fcr.to_csv(os.path.join(directories["output"], "feature_compression_ratios.tsv.gz"), sep="\t")
    df_mags["id_genome_cluster"].to_frame().dropna(how="any", axis=0).to_csv(os.path.join(directories["output"], "mags_to_slcs.tsv"), sep="\t", header=None)
    df_scaffolds["id_genome"].to_frame().dropna(how="any", axis=0).to_csv(os.path.join(directories["output"], "scaffolds_to_mags.tsv"), sep="\t", header=None)
    df_scaffolds["id_genome_cluster"].to_frame().dropna(how="any", axis=0).to_csv(os.path.join(directories["output"], "scaffolds_to_slcs.tsv"), sep="\t", header=None)
    df_proteins["id_protein_cluster"].to_frame().dropna(how="any", axis=0).to_csv(os.path.join(directories["output"], "proteins_to_orthogroups.tsv"), sep="\t", header=None) # Change labels?
    print(*map(lambda fp: " * {}".format(fp), glob.glob(os.path.join(directories["output"],"*.tsv")) + glob.glob(os.path.join(directories["output"],"*.faa"))), sep="\n", file=sys.stdout )


if __name__ == "__main__":
    main(sys.argv[1:])