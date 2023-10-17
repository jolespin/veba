#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, shutil, time, gzip, warnings
from multiprocessing import cpu_count
from collections import OrderedDict, defaultdict
from tqdm import tqdm

import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser 

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.10.6"

PROTEIN_MINIMUM_IDENTITY_THRESHOLD = 50.0
NUCLEOTIDE_MINIMUM_IDENTITY_THRESHOLD = 75.0

def get_basename(x):
    _, fn = os.path.split(x)
    if fn.endswith(".gz"):
        fn = fn[:-3]
    return ".".join(fn.split(".")[:-1])

def get_marker_gene_cluster_prevalence(df_input:pd.DataFrame): # genome, protein_cluster, marker_gene_cluster
    # Read Input
    genome_clusters = sorted(df_input.iloc[:,0].unique())
    protein_clusters = sorted(df_input.iloc[:,1].unique())
    marker_gene_clusters = sorted(df_input.iloc[:,2].unique())

    # Create array
    A = np.zeros((len(genome_clusters), len(marker_gene_clusters)), dtype=int)

    for _, (id_genome_cluster, id_protein_cluster, id_marker_gene_cluster) in df_input.iterrows():
        i = genome_clusters.index(id_genome_cluster)
        j = marker_gene_clusters.index(id_marker_gene_cluster)
        A[i,j] += 1

    # Create output
    df_output = pd.DataFrame(A, index=genome_clusters, columns=marker_gene_clusters)
    df_output.index.name = "id_genome_cluster"
    df_output.columns.name = "id_marker_gene_cluster"

    return df_output

# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = {
        "edgelist_to_clusters.py",
        "mmseqs2_wrapper.py",
        # "table_to_fasta.py",
    }

    required_executables={
        # "fastANI",
        "mmseqs",
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
    
    assert_acceptable_arguments(opts.algorithm, {"easy-cluster", "easy-linclust"})
    assert_acceptable_arguments(opts.sequence_space, {"protein", "nucleotide"})
    if opts.minimum_identity_threshold is None:
        if opts.sequence_space == "protein":
            opts.minimum_identity_threshold = PROTEIN_MINIMUM_IDENTITY_THRESHOLD
        if opts.sequence_space == "nucleotide":
            opts.minimum_identity_threshold = NUCLEOTIDE_MINIMUM_IDENTITY_THRESHOLD
    
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
    parser_io.add_argument("-i", "--core_pangenomes_table", type=str, default="stdin",  help = "path/to/core_pangenomes_table.tsv, Format: Must include the follow columns (No header) [id_genome-cluster]<tab>[core_pangenome(protein.fasta)]<tab>[core_pangenome(cds.fasta)] but can include additional columns to the right (e.g., [metadata]).  Suggested input is from `compile_core_pangenome_table.py` script. [Default: stdin]")
    parser_io.add_argument("-n", "--protein_cluster_table", type=str, required=True,  help = "path/to/protein_clusters.tsv, Format: Must include the following column [average_number_of_copies_per_genome] but can include additional columns.  Suggested input is `protein_clusters.tsv` from the `cluster.py` module.")
    parser_io.add_argument("-o","--output_directory", type=str, default="veba_output/markers", help = "path/to/project_directory [Default: veba_output/markers]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
    # parser_utility.add_argument("--verbose", action='store_true')

    # MMSEQS2
    parser_mmseqs2 = parser.add_argument_group('MMSEQS2 arguments')
    parser_mmseqs2.add_argument("-a", "--algorithm", type=str, default="easy-cluster", help="MMSEQS2 | {easy-cluster, easy-linclust} [Default: easy-cluster]")
    parser_mmseqs2.add_argument("-s", "--sequence_space", type=str, default="nucleotide", help="MMSEQS2 | {protein, nucleotide} [Default: nucleotide]")
    parser_mmseqs2.add_argument("-t", "--minimum_identity_threshold", type=float, help="MMSEQS2 | SLC-Specific Marker Cluster percent identity threshold (Range (0.0, 100.0]) [Default: {} for protein-space and {} for nucleotide-space]".format(PROTEIN_MINIMUM_IDENTITY_THRESHOLD, NUCLEOTIDE_MINIMUM_IDENTITY_THRESHOLD))
    parser_mmseqs2.add_argument("-c", "--minimum_coverage_threshold", type=float, default=0.8, help="MMSEQS2 | SSPC coverage threshold (Range (0.0, 1.0]) [Default: 0.8]")
    parser_mmseqs2.add_argument("--marker_cluster_prefix", type=str, default="MGC-", help="Marker cluster prefix [Default: 'MGC-]")
    parser_mmseqs2.add_argument("--marker_cluster_suffix", type=str, default="", help="Marker cluster suffix [Default: '")
    parser_mmseqs2.add_argument("--marker_cluster_prefix_zfill", type=int, default=0, help="Cluster prefix zfill. Use 7 to match identifiers from OrthoFinder.  Use 0 to add no zfill. [Default: 0]") #7
    parser_mmseqs2.add_argument("--mmseqs2_options", type=str, default="", help="MMSEQS2 | More options (e.g. --arg 1 ) [Default: '']")

    # Prevalence
    parser_prevalence = parser.add_argument_group('Prevalence arguments')
    parser_prevalence.add_argument("-M", "--maximum_number_of_copies_per_genome", type=float, default=1.0, help="Maximum number of copies per genome for a marker gene to be considered [Default: 1.0]")
    parser_prevalence.add_argument("-m", "--minimum_number_of_copies_per_genome", type=float, default=1.0, help="Minimum number of copies per genome for a marker gene to be considered.  If coreness is relaxed (not recommended) then this will need to be adjusted [Default: 1.0]")

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
    # if not opts.no_representative_sequences:
    directories["marker_sequences"] = create_directory(os.path.join(directories["output"], "marker_sequences"))
    # directories["serialization"] = create_directory(os.path.join(directories["output"], "serialization"))
    # directories[("misc", "pangenomes", "networkx_graphs")] = create_directory(os.path.join(directories[("misc", "pangenomes")], "networkx_graphs"))
    directories["log"] = create_directory(os.path.join(directories["project"], "log"))
    directories["tmp"] = create_directory(os.path.join(directories["project"], "tmp"))
    directories["checkpoints"] = create_directory(os.path.join(directories["project"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["project"], "intermediate"))
    os.environ["TMPDIR"] = directories["tmp"]

    # Make directories
    t0 = time.time()
    print(format_header(" " .join(["* ({}) Creating directories:".format(format_duration(t0)), directories["intermediate"]])), file=sys.stdout)
    os.makedirs(opts.output_directory, exist_ok=True)

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

    # Load core pangenome table
    if opts.core_pangenomes_table == "stdin":
        opts.core_pangenomes_table = sys.stdin
    df_core_pangenomes = pd.read_csv(opts.core_pangenomes_table, sep="\t", header=None)
    assert df_core_pangenomes.shape[1] >= 3, "Must include the follow columns (No header) [id_genome_cluster]<tab>[proteins]<tab>[cds] but can include additional columns to the right (e.g., [metadata]).  Suggested input is from `compile_core_pangenome_table.py` script.  You have provided a table with {} columns.".format(df_core_pangenomes.shape[1])
    df_core_pangenomes = df_core_pangenomes.iloc[:,:3]
    df_core_pangenomes.columns = ["id_genome_cluster",  "proteins", "cds"]
    df_core_pangenomes = df_core_pangenomes.set_index("id_genome_cluster")
    assert not np.any(df_core_pangenomes.isnull()), "--core_pangenomes_table has missing values.  Please correct this.\n{}".format(df_core_pangenomes.loc[df_core_pangenomes.isnull().sum(axis=1)[lambda x: x > 0].index].to_string())

    # Load protein clusters
    df_protein_clusters = pd.read_csv(opts.protein_cluster_table, sep="\t", index_col=0)
    assert "average_number_of_copies_per_genome" in df_protein_clusters.columns, "--protein_cluster_table must have a column called `average_number_of_copies_per_genome` where each entry is a float ≥ 1.0"
    proteincluster_to_copy_number = df_protein_clusters["average_number_of_copies_per_genome"]
    # assert np.all(proteincluster_to_copies >= opts.minimum_number_of_copies_per_genome), "Some of the protein clusters are in less than 1 copy per genome"

    # Check overlap between inputs
    cluster_to_protein_sequence = dict()
    cluster_to_nucleotide_sequence = dict()
    proteincluster_to_genomecluster = dict()

    for id_genome_cluster, row in tqdm(df_core_pangenomes.iterrows(), "Reading panproteomes from genome clusters", total=df_core_pangenomes.shape[0], file=sys.stdout, unit=" Pangenomes"):

        # Proteins
        with get_file_object(row["proteins"], mode="read", compression="infer", safe_mode="infer", verbose=False) as f:
            for header, seq in SimpleFastaParser(f):
                id = header.split(" ")[0]
                assert id not in cluster_to_protein_sequence, f"{id} from {id_genome_cluster} is a duplicate in protein-space"
                cluster_to_protein_sequence[id] = seq
                proteincluster_to_genomecluster[id] = id_genome_cluster
        # Nucleotides
        with get_file_object(row["cds"], mode="read", compression="infer", safe_mode="infer", verbose=False) as f:
            for header, seq in SimpleFastaParser(f):
                id = header.split(" ")[0]
                assert id not in cluster_to_nucleotide_sequence, f"{id} from {id_genome_cluster} is a duplicate in nucleotide-space"
                cluster_to_nucleotide_sequence[id] = seq

    # Overlap between protein and nucleotide space
    clusters_in_protein_space = set(cluster_to_protein_sequence.keys()) 
    clusters_in_nucleotide_space = set(cluster_to_nucleotide_sequence.keys()) 
    A = len(clusters_in_protein_space - clusters_in_nucleotide_space)
    B = len(clusters_in_nucleotide_space - clusters_in_protein_space)
    assert clusters_in_protein_space == clusters_in_nucleotide_space, "Identifiers from core pangenomes table do not overlap in protein and nucleotide space.\n\nNumber of unique identifiers in protein space: {} \nNumber of unique identifiers in nucleotide space: {}".format(A, B)
    proteincluster_identifiers = clusters_in_protein_space | clusters_in_nucleotide_space
    del clusters_in_protein_space
    del clusters_in_nucleotide_space

    # Overlap between core protein clusters and copy number table
    clusters_with_average_copy_number = set(proteincluster_to_copy_number.index)
    assert proteincluster_identifiers <= clusters_with_average_copy_number, "Not all core protein clusters from --core_pangenomes_table are in --protein_clusters: {} protein clusters".format(len(proteincluster_identifiers - clusters_with_average_copy_number))
    del clusters_with_average_copy_number
    proteincluster_identifiers = list(proteincluster_identifiers)

    # Index only the core protein clusters
    proteincluster_to_copy_number = proteincluster_to_copy_number[proteincluster_identifiers]

    # Remove protein clusters that aren't within the range of accepted copy numbers
    f_passed = open(os.path.join(directories["intermediate"], "passed_qc.list"), "w")
    f_failed = open(os.path.join(directories["intermediate"], "failed_qc.list"), "w")

    for id, n in proteincluster_to_copy_number.items():
        if opts.minimum_number_of_copies_per_genome <= n <= opts.maximum_number_of_copies_per_genome:
            print(id, file=f_passed)
        else:
            print(id, file=f_failed)
            del cluster_to_protein_sequence[id]
            del cluster_to_nucleotide_sequence[id]
            del proteincluster_to_copy_number[id]
            del proteincluster_to_genomecluster[id]

    f_passed.close()
    f_failed.close()

    cluster_to_protein_sequence = pd.Series(cluster_to_protein_sequence)
    cluster_to_nucleotide_sequence = pd.Series(cluster_to_nucleotide_sequence)

    # Write sequences to use for clustering
    if opts.sequence_space == "protein":
        write_fasta(cluster_to_protein_sequence, os.path.join(directories["tmp"], "sequences.fasta"))
    if opts.sequence_space == "nucleotide":
        write_fasta(cluster_to_nucleotide_sequence, os.path.join(directories["tmp"], "sequences.fasta"))
        
     # MMSEQS2
    f_cmds = open(os.path.join(opts.output_directory, "commands.sh"), "w")

    print(format_header(" * ({}) Running MMSEQS2:".format(format_duration(t0))), file=sys.stdout)

    # Run MMSEQS2
    name = "mmseqs2__{}-space".format(opts.sequence_space)
    description = "[Program = MMSEQS2] [Sequence Space = {}]".format(opts.sequence_space)

    cmd = Command([
        os.environ["mmseqs2_wrapper.py"],
        "--fasta {}".format(os.path.join(directories["tmp"], "sequences.fasta" )),
        "--output_directory {}".format(directories["intermediate"]),
        "--algorithm {}".format(opts.algorithm),
        "--n_jobs {}".format(opts.n_jobs),
        "--minimum_identity_threshold {}".format(opts.minimum_identity_threshold),
        "--minimum_coverage_threshold {}".format(opts.minimum_coverage_threshold),
        "--mmseqs2_options='{}'" if bool(opts.mmseqs2_options) else "",
        "--cluster_prefix {}".format(opts.marker_cluster_prefix),
        "--cluster_suffix {}".format(opts.marker_cluster_suffix) if bool(opts.marker_cluster_suffix) else "",
        "--cluster_prefix_zfill {}".format(opts.marker_cluster_prefix_zfill),
        "--basename marker_clusters",
        "--identifiers {}".format(os.path.join(directories["intermediate"], "passed_qc.list")),
        "--no_sequences_and_header",
        "--representative_output_format table",

            "&&",

        "rm -rf",
        os.path.join(directories["tmp"], "sequences.fasta" ),
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

    f_cmds.close()

    # Marker gene cluster prevalence
    print(format_header(" * ({}) Calculating marker prevalence:".format(format_duration(t0))), file=sys.stdout)

    marker_to_representative = pd.read_csv(os.path.join(directories["intermediate"],"output", "representative_sequences.tsv.gz"), sep="\t", index_col=0, header=None).iloc[:,0]

    df_markers = pd.read_csv(os.path.join(directories["intermediate"],"output", "marker_clusters.tsv"), sep="\t", header=None)
    df_markers.columns = ["id_protein_cluster", "id_marker_gene_cluster"]
    df_markers.insert(0, "id_genome_cluster", df_markers["id_protein_cluster"].map(lambda x: proteincluster_to_genomecluster[x]))

    df_prevalence = get_marker_gene_cluster_prevalence(df_markers)
    df_prevalence.to_csv(os.path.join(directories["output"], "prevalence_table.tsv.gz"), sep="\t")

    maximum_prevalence = df_prevalence.max(axis=0)
    maximum_prevalence_gt1 = maximum_prevalence > 1

    marker_proteins_with_maximum_prevalence_gt1 = maximum_prevalence.index[maximum_prevalence_gt1]

    if len(marker_proteins_with_maximum_prevalence_gt1) > 1:
        warnings.warn(
            """
            When clustering at {} percent identity and {} coverage in {}-space, {}/{} of the markers have a prevalence > 1.
            You may need to increase your --minimum_identity_threshold or --minimum_coverage_threshold or switch sequence-spaces.  
            If you used VEBA for pangenomes, then the genes were clustered in protein-space and there isn't a 1-to-1 overlap between protein-space and nucleotide-space cutoffs.
            """.format(
                opts.minimum_identity_threshold,
                opts.minimum_coverage_threshold,
                opts.sequence_space,
                len(marker_proteins_with_maximum_prevalence_gt1),
                maximum_prevalence.size,
            )
        )

    with open(os.path.join(directories["intermediate"], "multiple_copy_marker_gene_clusters.list"), "w") as f:
        for id_marker_cluster in marker_proteins_with_maximum_prevalence_gt1:
            print(id_marker_cluster, file=f)

    df_prevalence = df_prevalence.drop(marker_proteins_with_maximum_prevalence_gt1, axis=1)

    df_prevalence = df_prevalence.applymap(lambda x: opts.minimum_number_of_copies_per_genome <= x <= opts.maximum_number_of_copies_per_genome).astype(int)
    marker_clusters = df_prevalence.sum(axis=0)[lambda x: x == 1].index
    markercluster_to_genomecluster = df_prevalence.loc[:,marker_clusters].idxmax(axis=0)

    # Identifier mapping for protein clusters
    print(format_header(" * ({}) Writing output tables:".format(format_duration(t0))), file=sys.stdout)

    proteincluster_to_markercluster = pd.read_csv(os.path.join(directories["intermediate"], "output", "marker_clusters.tsv"), sep="\t", index_col=0, header=None).iloc[:,0]
    df_proteinclusters = proteincluster_to_markercluster.loc[proteincluster_to_markercluster.map(lambda x: x in marker_clusters)].to_frame("id_marker_cluster")
    df_proteinclusters["id_genome_cluster"] = df_proteinclusters["id_marker_cluster"].map(lambda x: markercluster_to_genomecluster[x])
    df_proteinclusters.index.name = "id_protein_cluster"
    df_proteinclusters.to_csv(os.path.join(directories["output"], "identifier_mapping.protein_clusters.tsv"), sep="\t")

    # Identifier mapping for marker clusters
    df_markers = markercluster_to_genomecluster.to_frame("id_genome_cluster")
    df_markers["id_protein_cluster_representative"] = marker_to_representative[df_markers.index]
    df_markers.to_csv(os.path.join(directories["output"], "identifier_mapping.marker_clusters.tsv"), sep="\t")

    # Pangenome clusters
    df_pangenome_markers = markercluster_to_genomecluster.groupby(markercluster_to_genomecluster).apply(lambda x: set(x.index)).to_frame("markers")
    df_pangenome_markers.insert(0, "number_of_markers", df_pangenome_markers["markers"].map(len))
    df_pangenome_markers.index.name = "id_genome_cluster"
    df_pangenome_markers.to_csv(os.path.join(directories["output"], "pangenome_markers.tsv"), sep="\t") 

    # Marker sequences
    print(format_header(" * ({}) Writing marker gene sequences:".format(format_duration(t0))), file=sys.stdout)
    for id_genome_cluster, marker_clusters in tqdm(df_pangenome_markers["markers"].items(), total=df_pangenome_markers.shape[0], unit=" Pangenomes", file=sys.stdout):
        f_proteins = open(os.path.join(directories["marker_sequences"], f"{id_genome_cluster}.faa"), "w")
        f_cds = open(os.path.join(directories["marker_sequences"], f"{id_genome_cluster}.ffn"), "w")
        
        for id_marker_cluster in sorted(marker_clusters):
            id_protein_cluster_representative = marker_to_representative[id_marker_cluster]
            id_genome_cluster = markercluster_to_genomecluster[id_marker_cluster]
            header = "{} {} {}".format(id_marker_cluster, id_genome_cluster, id_protein_cluster_representative)

            # Proteins
            print(">{}\n{}".format(header, cluster_to_protein_sequence[id_protein_cluster_representative]), file=f_proteins)
            print(">{}\n{}".format(header, cluster_to_nucleotide_sequence[id_protein_cluster_representative]), file=f_cds)
        f_proteins.close()
        f_cds.close()

if __name__ == "__main__":
    main(sys.argv[1:])