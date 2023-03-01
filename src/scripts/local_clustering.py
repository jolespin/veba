#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, shutil, time
from multiprocessing import cpu_count
from collections import OrderedDict, defaultdict

import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser 

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.2.25"

def get_basename(x):
    _, fn = os.path.split(x)
    if fn.endswith(".gz"):
        fn = fn[:-3]
    return ".".join(fn.split(".")[:-1])

# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = {
        "edgelist_to_clusters.py",
    }

    required_executables={
        "fastANI",
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
    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -f <scaffolds.fasta> -d <metaeuk_database> -i <scaffolds_to_bins.tsv>  -o <output_directory>".format(__program__)

    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i", "--input", type=str, default="stdin",  help = "path/to/input.tsv, Format: Must include the follow columns (No header) [organism_type]<tab>[id_sample]<tab>[id_mag]<tab>[genome]<tab>[proteins] but can include additional columns to the right (e.g., [cds]<tab>[gene_models]).  Suggested input is from `compile_genomes_table.py` script. [Default: stdin]")
    parser_io.add_argument("-o","--output_directory", type=str, default="local_clustering_output", help = "path/to/project_directory [Default: local_clustering_output]")
    parser_io.add_argument("-e", "--no_singletons", action="store_true", help="Exclude singletons") #isPSLC-1_SSPC-3345__SRR178126

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
    # parser_utility.add_argument("--verbose", action='store_true')

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
        opts.n_jobs = cpu_count()
    assert opts.n_jobs >= 1, "--n_jobs must be ≥ 1 (or -1 to use all available threads)"

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

    # Load input
    t0 = time.time()
    if opts.input == "stdin":
        opts.input = sys.stdin
    df_genomes = pd.read_csv(opts.input, sep="\t", header=None)
    assert df_genomes.shape[1] >= 5, "Must include the follow columns (No header) [organism_type]<tab>[id_sample]<tab>[id_mag]<tab>[genome]<tab>[proteins] but can include additional columns to the right (e.g., [cds]<tab>[gene_models]).  Suggested input is from `compile_genomes_table.py` script.  You have provided a table with {} columns.".format(df_genomes.shape[1])
    df_genomes = df_genomes.iloc[:,:5]
    df_genomes.columns = ["organism_type", "id_sample", "id_mag", "genome", "proteins"]
    assert not np.any(df_genomes.isnull()), "Input has missing values.  Please correct this.\n{}".format(df_genomes.loc[df_genomes.isnull().sum(axis=1)[lambda x: x > 0].index].to_string())

    # Make directories
    print(format_header(" " .join(["* ({}) Creating directories:".format(format_duration(t0)), directories["intermediate"]])), file=sys.stdout)
    os.makedirs(opts.output_directory, exist_ok=True)

    for organism_type, df_1 in df_genomes.groupby("organism_type"):
        # Organism directories
        os.makedirs(os.path.join(directories["intermediate"], organism_type), exist_ok=True)

        # Sample directories
        for id_sample, df_2 in df_1.groupby("id_sample"):
            os.makedirs(os.path.join(directories["intermediate"], organism_type, id_sample), exist_ok=True)
            os.makedirs(os.path.join(directories["intermediate"], organism_type, id_sample, "clusters"), exist_ok=True)

            # Write genomes to file
            with open(os.path.join(directories["intermediate"], organism_type, id_sample, "genomes.list"), "w") as f:
                print(*sorted(df_2["genome"]), sep="\n", file=f)

            # Write genome identifiers to file
            with open(os.path.join(directories["intermediate"], organism_type, id_sample, "genome_identifiers.list"), "w") as f:
                print(*sorted(df_2["genome"].map(get_basename)), sep="\n", file=f)


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

    mag_to_numberofscaffolds = pd.Series(mag_to_numberofscaffolds)
    mag_to_numberofproteins = pd.Series(mag_to_numberofproteins)
    mag_to_sample = pd.Series(mag_to_sample)
    mag_to_organismtype = pd.Series(dict(zip(df_genomes["id_mag"], df_genomes["organism_type"])))
    scaffold_to_mag = pd.Series(scaffold_to_mag)
    scaffold_to_sample = pd.Series(scaffold_to_sample)
    protein_to_mag = pd.Series(protein_to_mag)
    protein_to_sample = pd.Series(protein_to_sample)
    protein_to_sequence = pd.Series(protein_to_sequence)

    # Commands
    f_cmds = open(os.path.join(opts.output_directory, "commands.sh"), "w")

    # FastANI
    print(format_header("* ({}) Running FastANI:".format(format_duration(t0))), file=sys.stdout)
    for fp in pv(glob.glob(os.path.join(directories["intermediate"], "*", "*", "genomes.list")), "Running FastANI"):
        fields = fp.split("/")
        organism_type = fields[-3]
        id_sample =  fields[-2]
        
        name = "fastani__{}__{}".format(organism_type, id_sample)
        description = "[Program = FastANI] [Organism_Type = {}] [Sample_ID = {}]".format(organism_type, id_sample)
        cmd = Command([
            os.environ["fastANI"],
            "-t {}".format(opts.n_jobs),
            "--rl {}".format(fp),
            "--ql {}".format(fp),
            "-o {}".format(os.path.join(os.path.split(fp)[0], "fastani_output.tsv")),
            opts.fastani_options,

                "&&",

            "cat",
            os.path.join(os.path.split(fp)[0], "fastani_output.tsv"),
            "|",
            "cut -f1-3",
            "|",
            os.environ["edgelist_to_clusters.py"],
            "--basename",
            "-t {}".format(opts.ani_threshold),
            "--no_singletons" if bool(opts.no_singletons) else "",
            "--cluster_prefix {}__{}{}".format(id_sample, organism_type[0].upper(), opts.genome_cluster_prefix),
            "--cluster_suffix {}".format(opts.genome_cluster_suffix) if bool(opts.genome_cluster_suffix) else "",
            "--cluster_prefix_zfill {}".format(opts.genome_cluster_prefix_zfill),
            "-o {}".format(os.path.join(os.path.split(fp)[0], "genome_clusters.tsv")),
            "--identifiers {}".format(os.path.join(directories["intermediate"], organism_type, id_sample, "genome_identifiers.list")),

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

    # MMSEQS2
    print(format_header(" * ({}) Running MMSEQS2:".format(format_duration(t0))), file=sys.stdout)
    mag_to_genomecluster = dict()
    protein_to_proteincluster = dict()
    for fp in pv(glob.glob(os.path.join(directories["intermediate"], "*", "*", "genome_clusters.tsv")), "Running MMSEQS2"):
        fields = fp.split("/")
        organism_type = fields[-3]
        id_sample =  fields[-2]

        mag_to_genomecluster_within_sample = pd.read_csv(fp, sep="\t", index_col=0, header=None).iloc[:,0]
        mag_to_genomecluster.update(mag_to_genomecluster_within_sample)

        for id_genomecluster, data in mag_to_genomecluster_within_sample.groupby(mag_to_genomecluster_within_sample):
            genomecluster_directory = os.path.join(directories["intermediate"], organism_type, id_sample, "clusters", id_genomecluster)
            os.makedirs(genomecluster_directory, exist_ok=True)
            # if not os.path.exists(os.path.join(genomecluster_directory, "protein_clusters.tsv")):

            # Get MAGs and proteins
            mags = data.index
            proteins = protein_to_mag[protein_to_mag.map(lambda x: x in mags)].index 
            with open(os.path.join(genomecluster_directory, "genomes.list" ), "w") as f:
                print(*sorted(mags), sep="\n", file=f)
            with open(os.path.join(genomecluster_directory, "protein_identifiers.list"), "w") as f:
                print(*proteins, sep="\n", file=f)
            write_fasta(protein_to_sequence[proteins], os.path.join(genomecluster_directory, "proteins.faa" ))

            # Run MMSEQS2
            name = "mmseqs2__{}__{}".format(organism_type, id_genomecluster)
            description = "[Program = MMSEQS2] [Organism_Type = {}] [Sample_ID = {}] [Genome_Cluster = {}]".format(organism_type, id_sample, id_genomecluster)
            cmd = Command([
                os.environ["mmseqs"],
                "easy-cluster",
                os.path.join(genomecluster_directory, "proteins.faa" ),
                os.path.join(genomecluster_directory, "mmseqs2"),
                directories["tmp"],
                "--threads {}".format(opts.n_jobs),
                "--min-seq-id {}".format(opts.minimum_identity_threshold/100),
                "-c {}".format(opts.minimum_coverage_threshold),
                "--cov-mode 1",
                "--dbtype 1",
                opts.mmseqs2_options,

                    "&&",

                os.environ["edgelist_to_clusters.py"],
                "-i {}".format(os.path.join(genomecluster_directory, "mmseqs2_cluster.tsv")),
                "--no_singletons" if bool(opts.no_singletons) else "",
                "--cluster_prefix {}_{}".format(id_genomecluster, opts.protein_cluster_prefix),
                "--cluster_suffix {}".format(opts.protein_cluster_suffix) if bool(opts.protein_cluster_suffix) else "",
                "--cluster_prefix_zfill {}".format(opts.protein_cluster_prefix_zfill),
                "-o {}".format(os.path.join(genomecluster_directory, "protein_clusters.tsv")),
                "--identifiers {}".format(os.path.join(genomecluster_directory, "protein_identifiers.list")),

                    "&&",

                "rm -rf",
                os.path.join(genomecluster_directory, "mmseqs2_all_seqs.fasta"),
                os.path.join(genomecluster_directory, "proteins.faa"),
                os.path.join(directories["tmp"], "*"),

                    "&&",

                "gzip",
                os.path.join(genomecluster_directory, "mmseqs2_rep_seq.fasta"),

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
            protein_to_proteinluster_within_sample = pd.read_csv(os.path.join(genomecluster_directory, "protein_clusters.tsv"), sep="\t", index_col=0, header=None).iloc[:,0]
            protein_to_proteincluster.update(protein_to_proteinluster_within_sample.to_dict())
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
    for organism_type, df_1 in df_mags.groupby("organism_type"):
        for id_sample, df_2 in df_1.groupby("sample_of_origin"):
            feature_to_cluster = df_2["id_genome_cluster"].dropna()
            fcr_data[(organism_type, id_sample)]["genomic_fcr"] = 1 - (feature_to_cluster.nunique()/feature_to_cluster.size)

    for organism_type, df_1 in df_proteins.groupby("organism_type"):
        for id_sample, df_2 in df_1.groupby("sample_of_origin"):
            feature_to_cluster = df_2["id_protein_cluster"].dropna()
            fcr_data[(organism_type, id_sample)]["functional_fcr"] = 1 - (feature_to_cluster.nunique()/feature_to_cluster.size)
    df_fcr = pd.DataFrame(fcr_data).T.sort_index()
    df_fcr.index.names = ["organism_type", "id_sample"]
    
    # Writing output files
    print(format_header(" * ({}) Writing Output Tables:".format(format_duration(t0))), file=sys.stdout)
    df_mags.to_csv(os.path.join(directories["output"], "identifier_mapping.genomes.tsv"), sep="\t")
    df_scaffolds.to_csv(os.path.join(directories["output"], "identifier_mapping.scaffolds.tsv"), sep="\t")
    df_proteins.to_csv(os.path.join(directories["output"], "identifier_mapping.proteins.tsv"), sep="\t")
    df_genomeclusters.to_csv(os.path.join(directories["output"], "genome_clusters.tsv"), sep="\t")
    df_proteinclusters.to_csv(os.path.join(directories["output"], "protein_clusters.tsv"), sep="\t")
    df_fcr.to_csv(os.path.join(directories["output"], "feature_compression_ratios.tsv"), sep="\t")
    df_mags["id_genome_cluster"].to_frame().dropna(how="any", axis=0).to_csv(os.path.join(directories["output"], "mags_to_slcs.tsv"), sep="\t", header=None)
    df_scaffolds["id_genome_cluster"].to_frame().dropna(how="any", axis=0).to_csv(os.path.join(directories["output"], "scaffolds_to_slcs.tsv"), sep="\t", header=None)
    df_proteins["id_protein_cluster"].to_frame().dropna(how="any", axis=0).to_csv(os.path.join(directories["output"], "proteins_to_orthogroups.tsv"), sep="\t", header=None) # Change labels?
    print(*map(lambda fp: " * {}".format(fp), glob.glob(os.path.join(directories["output"],"*.tsv"))), sep="\n", file=sys.stdout )


if __name__ == "__main__":
    main(sys.argv[1:])
