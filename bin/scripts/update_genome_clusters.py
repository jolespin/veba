#!/usr/bin/env python
import sys
import os
import argparse
import tempfile
from collections import defaultdict
from tqdm import tqdm
import pandas as pd
from pyexeggutor import (
    open_file_reader,
    open_file_writer,
    read_json,
    write_json,
    build_logger,
    reset_logger,
    format_duration,
    format_header,
    format_bytes,
    get_timestamp,
    get_directory_tree,
    get_directory_size,
    get_filepath_basename,
    get_md5hash_from_directory,
    add_executables_to_environment,
    RunShellCommand,
)
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.11.21"
        
# Run skani dist non-viral
def run_skani_dist_nonviral_reference_based(logger, log_directory, opts, **arguments):
    
    with (
        # Temporary files
        tempfile.NamedTemporaryFile(mode="w") as f_ql,
        tempfile.NamedTemporaryFile(mode="w") as f_rl,
        # File readers
        open_file_reader(opts.query_genomes) as f_query,
        open_file_reader(opts.reference_genomes_with_clusters) as f_reference,
        ):
        
        # Write query genome list
        logger.info(f"[NamedTemporaryFile] Writing query genome list to temporary file: {f_ql.name}")
        for line in f_query:
            line = line.strip()
            if not line.startswith("#"):
                organism_type, id_genome_query, path_genome = line.split("\t")
                if organism_type == arguments["organism_type"]:
                    print(path_genome, file=f_ql)
        f_ql.flush()
        
        # Write reference genome list
        logger.info(f"[NamedTemporaryFile] Writing reference genome list to temporary file: {f_rl.name}")
        for line in f_reference:
            line = line.strip()
            if not line.startswith("#"):
                organism_type, id_genome_reference, id_genome_cluster, path_genome = line.split("\t")
                if organism_type == arguments["organism_type"]:
                    print(path_genome, file=f_rl)
        f_rl.flush()
        
        # Command
        cmd = RunShellCommand(
            command=[
                os.environ["skani"],
                "dist",
                "--ql",
                f_ql.name,
                "--rl",
                f_rl.name,
                "-o",
                os.path.join(arguments["working_directory"], f"skani-dist_reference.tsv"),
                "--ci" if not opts.skani_no_confidence_interval else "",
                "--min-af",
                opts.skani_minimum_af,
                "-s",
                opts.skani_target_ani,
                "-c",
                opts.skani_nonviral_compression_factor,
                "-m",
                opts.skani_nonviral_marker_kmer_compression_factor,
                f"--{opts.skani_nonviral_preset}" if opts.skani_nonviral_preset else "",
                opts.skani_nonviral_options,
            ],
            name="skani-dist_nonviral_reference_based",
            validate_input_filepaths=[
                f_ql.name,
                f_rl.name,
            ],
            validate_output_filepaths=[
                os.path.join(arguments["working_directory"],  f"skani-dist_reference.tsv"),
            ],
        )
        
        # Run
        logger.info(f"[{cmd.name}] running command: {cmd.command}")
        cmd.run()
        logger.info(f"[{cmd.name}] duration: {cmd.duration_}")
        logger.info(f"[{cmd.name}] peak memory: {format_bytes(cmd.peak_memory_)}")

        # Dump
        logger.info(f"[{cmd.name}] dumping stdout, stderr, and return code: {log_directory}")
        cmd.dump(log_directory)
        
        # Validate
        logger.info(f"[{cmd.name}] checking return code status: {cmd.returncode_}")
        cmd.check_status()
        return cmd
    
# Run skani dist viral
def run_skani_dist_viral_reference_based(logger, log_directory, opts, **arguments):
    with (
        tempfile.TemporaryFile() as f_ql,
        tempfile.TemporaryFile() as f_rl,
        open_file_reader(opts.query_genomes) as f_query,
        open_file_reader(opts.reference_genomes_with_clusters) as f_reference,
        ):
        
        # Write query genome list
        logger.info(f"[TemporaryFile] Witing query genome list to temporary file: {f_ql.name}")
        for line in f_query:
            line = line.strip()
            if not line.startswith("#"):
                organism_type, id_genome_query, path_genome = line.split("\t")
                print(path_genome, file=f_ql)
        f_ql.flush()
        
        # Write reference genome list
        logger.info(f"[NamedTemporaryFile] Writing reference genome list to temporary file: {f_rl.name}")
        for line in f_reference:
            line = line.strip()
            if not line.startswith("#"):
                organism_type, id_genome_reference, id_genome_cluster, path_genome = line.split("\t")
                if organism_type == arguments["organism_type"]:
                    print(path_genome, file=f_rl)
        f_rl.flush()
        
        # Command
        cmd = RunShellCommand(
            command=[
                os.environ["skani"],
                "dist",
                "--ql",
                f_ql.name,
                "--rl",
                f_rl.name,
                "-o",
                os.path.join(arguments["working_directory"], f"skani-dist_reference.tsv"),
                "--ci" if not opts.skani_no_confidence_interval else "",
                "--sparse",
                "--min-af",
                opts.skani_minimum_af,
                "-s",
                opts.skani_target_ani,
                "-c",
                opts.skani_viral_compression_factor,
                "-m",
                opts.skani_viral_marker_kmer_compression_factor,
                f"--{opts.skani_viral_preset}" if opts.skani_viral_preset else "",
                opts.skani_viral_options,
            ],
            name="skani-dist_viral_reference_based",
            validate_input_filepaths=[
                f_ql.name,
                f_rl.name,
            ],
            validate_output_filepaths=[
                os.path.join(arguments["working_directory"], f"skani-dist_reference.tsv"),
            ],
        )
        
        # Run
        logger.info(f"[{cmd.name}] running command: {cmd.command}")
        cmd.run()
        logger.info(f"[{cmd.name}] duration: {cmd.duration_}")
        logger.info(f"[{cmd.name}] peak memory: {format_bytes(cmd.peak_memory_)}")

        # Dump
        logger.info(f"[{cmd.name}] dumping stdout, stderr, and return code: {log_directory}")
        cmd.dump(log_directory)
        
        # Validate
        logger.info(f"[{cmd.name}] checking return code status: {cmd.returncode_}")
        cmd.check_status()
        return cmd
    
# Run skani triangle non-viral
def run_skani_triangle_nonviral_query_based(logger, log_directory, opts, **arguments):
    
    with (
        # Temporary files
        tempfile.NamedTemporaryFile(mode="w") as f_ql,
        # File readers
        open_file_reader(opts.query_genomes) as f_query,
        ):
        
        # Write query genome list
        logger.info(f"[NamedTemporaryFile] Writing query genome list to temporary file: {f_ql.name}")
        for line in f_query:
            line = line.strip()
            if not line.startswith("#"):
                organism_type, id_genome_query, path_genome = line.split("\t")
                conditions = [
                    organism_type == arguments["organism_type"],
                    id_genome_query in arguments["query_genomes"],
                ]
                if all(conditions):
                    print(path_genome, file=f_ql)
        f_ql.flush()
        
        # Command
        cmd = RunShellCommand(
            command=[
                os.environ["skani"],
                "triangle",
                "-l",
                f_ql.name,
                "-o",
                os.path.join(arguments["working_directory"], "skani-triangle_query.tsv"),
                "--ci" if not opts.skani_no_confidence_interval else "",
                "--sparse",
                "--min-af",
                opts.skani_minimum_af,
                "-s",
                opts.skani_target_ani,
                "-c",
                opts.skani_nonviral_compression_factor,
                "-m",
                opts.skani_nonviral_marker_kmer_compression_factor,
                f"--{opts.skani_nonviral_preset}" if opts.skani_nonviral_preset else "",
                opts.skani_nonviral_options,
            ],
            name="skani-triangle_nonviral_query_based",
            validate_input_filepaths=[
                f_ql.name,
            ],
            validate_output_filepaths=[
                os.path.join(arguments["working_directory"], "skani-triangle_query.tsv"),
            ],
        )
        
        # Run
        logger.info(f"[{cmd.name}] running command: {cmd.command}")
        cmd.run()
        logger.info(f"[{cmd.name}] duration: {cmd.duration_}")
        logger.info(f"[{cmd.name}] peak memory: {format_bytes(cmd.peak_memory_)}")

        # Dump
        logger.info(f"[{cmd.name}] dumping stdout, stderr, and return code: {log_directory}")
        cmd.dump(log_directory)
        
        # Validate
        logger.info(f"[{cmd.name}] checking return code status: {cmd.returncode_}")
        cmd.check_status()
        return cmd
    
# Run skani triangle viral
def run_skani_triangle_viral_query_based(logger, log_directory, opts, **arguments):
    with (
        # Temporary files
        tempfile.NamedTemporaryFile(mode="w") as f_ql,
        # File readers
        open_file_reader(opts.query_genomes) as f_query,
        ):
        
        # Write query genome list
        logger.info(f"[NamedTemporaryFile] Writing query genome list to temporary file: {f_ql.name}")
        for line in f_query:
            line = line.strip()
            if not line.startswith("#"):
                organism_type, id_genome_query, path_genome = line.split("\t")
                conditions = [
                    organism_type == arguments["organism_type"],
                    id_genome_query in arguments["query_genomes"],
                ]
                if all(conditions):
                    print(path_genome, file=f_ql)
        f_ql.flush()
        
        # Command
        cmd = RunShellCommand(
            command=[
                os.environ["skani"],
                "triangle",
                "-l",
                f_ql.name,
                "-o",
                os.path.join(arguments["working_directory"], "skani-triangle_query.tsv"),
                "--ci" if not opts.skani_no_confidence_interval else "",
                "--sparse",
                "--min-af",
                opts.skani_minimum_af,
                "-s",
                opts.skani_target_ani,
                "-c",
                opts.skani_viral_compression_factor,
                "-m",
                opts.skani_viral_marker_kmer_compression_factor,
                f"--{opts.skani_viral_preset}" if opts.skani_viral_preset else "",
                opts.skani_viral_options,
            ],
            name="skani-triangle_viral_query_based",
            validate_input_filepaths=[
                f_ql.name,
            ],
            validate_output_filepaths=[
                os.path.join(arguments["working_directory"], "skani-triangle_query.tsv"),
            ],
        )
        
        # Run
        logger.info(f"[{cmd.name}] running command: {cmd.command}")
        cmd.run()
        logger.info(f"[{cmd.name}] duration: {cmd.duration_}")
        logger.info(f"[{cmd.name}] peak memory: {format_bytes(cmd.peak_memory_)}")

        # Dump
        logger.info(f"[{cmd.name}] dumping stdout, stderr, and return code: {log_directory}")
        cmd.dump(log_directory)
        
        # Validate
        logger.info(f"[{cmd.name}] checking return code status: {cmd.returncode_}")
        cmd.check_status()
        return cmd
    
# Run skani triangle viral
def run_edgelist_to_clusters(logger, log_directory, opts, **arguments):

    with (
        # Temporary files
        tempfile.NamedTemporaryFile(mode="w") as f_genomes,
        ):
        
        # Write query genome list
        logger.info(f"[NamedTemporaryFile] Writing query genome list to temporary file: {f_genomes.name}")
        for line in arguments["query_genomes"]:
            print(line, file=f_genomes)
        f_genomes.flush()
    
        # Command
        cmd = RunShellCommand(
            command=[
                "cat",
                os.path.join(arguments["working_directory"], "skani-triangle_query.tsv"),
                "|",
                "cut",
                "-f",
                "1-5",
                "|",
                "tail",
                "-n",
                "+2",
                "|",
                os.environ["edgelist_to_clusters.py"],
                "--basename",
                "-t",
                opts.ani_threshold,
                "-a",
                opts.af_threshold,
                "--af_mode",
                opts.af_mode,
                "--no_singletons" if bool(opts.no_singletons) else "",
                "--cluster_prefix",
                "{}{}".format(arguments["organism_type"][0].upper(), opts.genome_cluster_prefix),
                "--cluster_suffix {}".format(opts.genome_cluster_suffix) if bool(opts.genome_cluster_suffix) else "",
                "--cluster_prefix_zfill",
                opts.genome_cluster_prefix_zfill,
                "-o",
                os.path.join(arguments["working_directory"], "skani-triangle_query.clusters.tsv"),
                "--identifiers",
                f_genomes.name,
                "--export_graph",
                os.path.join(arguments["working_directory"], "networkx_graph.pkl"),
                "--export_dict",
                os.path.join(arguments["working_directory"], "dict.pkl"),
                "--export_representatives",
                os.path.join(arguments["working_directory"], "representatives.tsv"),
                "--cluster_label_mode",
                opts.cluster_label_mode,
            ],
            name="skani-triangle_viral_query_based",
            validate_input_filepaths=[
                f_genomes.name,
                os.path.join(arguments["working_directory"], "skani-triangle_query.tsv"),
            ],
            validate_output_filepaths=[
                os.path.join(arguments["working_directory"], "skani-triangle_query.clusters.tsv"),
                os.path.join(arguments["working_directory"], "networkx_graph.pkl"),
                os.path.join(arguments["working_directory"], "dict.pkl"),
                os.path.join(arguments["working_directory"], "representatives.tsv"),
            ],
        )
        
        # Run
        logger.info(f"[{cmd.name}] running command: {cmd.command}")
        cmd.run()
        logger.info(f"[{cmd.name}] duration: {cmd.duration_}")
        logger.info(f"[{cmd.name}] peak memory: {format_bytes(cmd.peak_memory_)}")

        # Dump
        logger.info(f"[{cmd.name}] dumping stdout, stderr, and return code: {log_directory}")
        cmd.dump(log_directory)
        
        # Validate
        logger.info(f"[{cmd.name}] checking return code status: {cmd.returncode_}")
        cmd.check_status()
        return cmd
    


def main(args=None):
    # Options
    # =======
    # Path info
    python_executable = sys.executable
    bin_directory = "/".join(python_executable.split("/")[:-1])
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, sys.version.split(" ")[0], python_executable, script_filename)
    usage = f"{__program__} --fasta path/to/cds.fasta --feature_mapping path/to/features.tsv --genomes path/to/genomes.tsv  --index_directory path/to/leviathan_index/"
    epilog = "Copyright 2024 New Atlantis Labs (jolespin@newatlantis.io)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-q","--query_genomes", type=str,   help = "path/to/genomes.tsv [organism_type, id_genome, path/to/genome] (No header)")
    parser_io.add_argument("-r","--reference_genomes_with_clusters", type=str,   help = "path/to/genomes.tsv [organism_type, id_genome, id_cluster, path/to/genome] (No header)")
    parser_io.add_argument("-o","--output_directory", type=str, default="clustering_output", help = "path/to/output_directory/ [Default: clustering_output]")
    parser_io.add_argument("-e", "--no_singletons", action="store_true", help="Exclude singletons")
    parser_io.add_argument("--no_references", action="store_true", help="Exclude reference genomes from cluster output")

    # # Utilities
    # parser_utility = parser.add_argument_group('Utility arguments')
    # parser_utility.add_argument("-p","--n_jobs", type=int, default=1,  help = "Number of threads to use.  Use -1 for all available. [Default: 1]")

    # ANI
    parser_genome_clustering = parser.add_argument_group('Genome clustering arguments')
    # parser_genome_clustering.add_argument("-G", "--genome_clustering_algorithm", type=str,  choices={"skani"}, default="skani", help="Program to use for ANI calculations.  `skani` is faster and more memory efficient. For v1.0.0 - v1.3.x behavior, use `fastani`. [Default: skani]")
    parser_genome_clustering.add_argument("-A", "--ani_threshold", type=float, default=95.0, help="Species-level cluster (SLC) ANI threshold (Range (0.0, 100.0]) [Default: 95.0]")
    parser_genome_clustering.add_argument("-F", "--af_threshold", type=float, default=50.0, help="Species-level cluster (SLC) alignment fraction threshold. Only available if `skani` is used as --genome_clustering_algorithm.  Recommended settings: 15 (relaxed), 30 (OceanDNA), and 50 (GTDB-Tk). The higher the value the more singleton clusters will be returned.  (Range (0.0, 100.0]) [Default: 50.0]")
    parser_genome_clustering.add_argument("--af_mode", type=str, default="relaxed",  choices={"relaxed", "strict"}, help = "Minimum alignment fraction mode with either `relaxed = max([AF_ref, AF_query]) > minimum_af` or `strict = (AF_ref > minimum_af) & (AF_query > minimum_af)` [Default: relaxed]") 
    parser_genome_clustering.add_argument("--genome_cluster_prefix", type=str, default="SLC-", help="Cluster prefix [Default: 'SLC-")
    parser_genome_clustering.add_argument("--genome_cluster_suffix", type=str, default="", help="Cluster suffix [Default: '")
    parser_genome_clustering.add_argument("--genome_cluster_prefix_zfill", type=int, default=0, help="Cluster prefix zfill. Use 7 to match identifiers from OrthoFinder.  Use 0 to add no zfill. [Default: 0]") #7
    parser_genome_clustering.add_argument("--cluster_label_mode", type=str, default="md5", choices={"numeric", "random", "pseudo-random", "md5", "nodes"}, help="Cluster label. [Default: 'md5']")

    # skani
    parser_skani = parser.add_argument_group('Skani dist/triangle arguments')
    # parser_skani.add_argument("--skani_executable", type=str, help="skani executable [Default: $PATH]")

    parser_skani.add_argument("--skani_target_ani",  type=float, default=80, help="skani | If you set --skani_target_ani to --ani_threshold, you may screen out genomes ANI ≥ --ani_threshold [Default: 80]")
    parser_skani.add_argument("--skani_minimum_af",  type=float, default=15, help="skani | Minimum aligned fraction greater than this value. If you set --skani_minimum_af to --af_threshold, you may screen out genomes AF ≥ --af_threshold [Default: 15]")
    parser_skani.add_argument("--skani_no_confidence_interval",  action="store_true", help="skani | Output [5,95] ANI confidence intervals using percentile bootstrap on the putative ANI distribution")

    parser_skani = parser.add_argument_group('[Prokaryotic & Eukaryotic] Skani dist/triangle arguments')
    parser_skani.add_argument("--skani_nonviral_preset", type=str, default="medium", choices={"fast", "medium", "slow", "none"}, help="skani [Prokaryotic & Eukaryotic]| Use `none` if you are setting skani -c (compression factor) {fast, medium, slow, none} [Default: medium]")
    parser_skani.add_argument("--skani_nonviral_compression_factor", type=int, default=125,  help="skani [Prokaryotic & Eukaryotic]|  Compression factor (k-mer subsampling rate).	[Default: 125]")
    parser_skani.add_argument("--skani_nonviral_marker_kmer_compression_factor", type=int, default=1000,  help="skani [Prokaryotic & Eukaryotic] | Marker k-mer compression factor. Markers are used for filtering. [Default: 1000]")
    parser_skani.add_argument("--skani_nonviral_options", type=str, default="", help="skani [Prokaryotic & Eukaryotic] | More options for `skani triangle` (e.g. --arg 1 ) [Default: '']")

    parser_skani = parser.add_argument_group('[Viral] Skani dist/triangle arguments')
    parser_skani.add_argument("--skani_viral_preset", type=str, default="slow", choices={"fast", "medium", "slow", "none"}, help="skani | Use `none` if you are setting skani -c (compression factor) {fast, medium, slow, none} [Default: slow]")
    parser_skani.add_argument("--skani_viral_compression_factor", type=int, default=30,  help="skani [Viral] | Compression factor (k-mer subsampling rate).	[Default: 30]")
    parser_skani.add_argument("--skani_viral_marker_kmer_compression_factor", type=int, default=200,  help="skani [Viral] | Marker k-mer compression factor. Markers are used for filtering. Consider decreasing to ~200-300 if working with small genomes (e.g. plasmids or viruses). [Default: 200]")
    parser_skani.add_argument("--skani_viral_options", type=str, default="", help="skani [Viral] | More options for `skani triangle` (e.g. --arg 1 ) [Default: '']")
  
    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # logger
    logger = build_logger("Genome Cluster Updater")

    # Commands
    logger.info(f"Command: {sys.argv}")
     
    # # Threads
    # if opts.n_jobs == -1:
    #     from multiprocessing import cpu_count 
    #     opts.n_jobs = cpu_count()
    #     logger.info(f"Setting --n_jobs to maximum threads {opts.n_jobs}")

    # assert opts.n_jobs >= 1, "--n_jobs must be ≥ 1.  To select all available threads, use -1."
    
    # Executables
    add_executables_to_environment([
        "skani", 
        "edgelist_to_clusters.py",
        ])
    
    # Output
    os.makedirs(opts.output_directory, exist_ok=True)
    os.makedirs(os.path.join(opts.output_directory, "output"), exist_ok=True)
    os.makedirs(os.path.join(opts.output_directory, "intermediate"), exist_ok=True)
    os.makedirs(os.path.join(opts.output_directory, "logs"), exist_ok=True)
    # os.makedirs(os.path.join(opts.output_directory, "tmp"), exist_ok=True)
    
    # Checks
    if opts.cluster_label_mode ==  "numeric":
        raise ValueError("--cluster_label_mode numeric is not currently support because if there any genomes that need new clusters, there will be duplicate labels with current implementation.")
    
    # Query organism_types
    logger.info(f"Reading query genomes: {opts.query_genomes}")
    query_organism_types = set()
    query_organism_type_to_genomes = defaultdict(set)
    with open_file_reader(opts.query_genomes) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                organism_type, id_genome_query, path_genome = line.split("\t")
                query_organism_types.add(organism_type)
                query_organism_type_to_genomes[organism_type].add(id_genome_query)
                
    # Reference clusters
    logger.info(f"Reading reference genome clusters: {opts.reference_genomes_with_clusters}")
    genome_to_cluster = dict()
    reference_organism_types = set()
    with open_file_reader(opts.reference_genomes_with_clusters) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                organism_type, id_genome_reference, id_genome_cluster, path_genome = line.split("\t")
                genome_to_cluster[id_genome_reference] = id_genome_cluster
                reference_organism_types.add(organism_type)

    # Phase I          
    logger.info(f"Phase I: Identifying representatives in existing clusters")

    # Identifying representatives in existing clusters
    for organism_type in query_organism_types:
        if organism_type not in {"prokaryotic", "eukaryotic", "viral"}:
            msg = f"Invalid organism_type: {organism_type}"
            logger.critical(msg)
            raise ValueError(msg)

        arguments = {
            "organism_type": organism_type,
            "working_directory":os.path.join(opts.output_directory, "intermediate", organism_type),
        }
        os.makedirs(arguments["working_directory"], exist_ok=True)
        
        if organism_type in {"prokaryotic", "eukaryotic"}:
            run_skani_dist_nonviral_reference_based(
                logger=logger, 
                log_directory=os.path.join(opts.output_directory, "logs"), 
                opts=opts, 
                **arguments,
            )
        if organism_type in {"viral"}:
            run_skani_dist_viral_reference_based(
                logger=logger, 
                log_directory=os.path.join(opts.output_directory, "logs"), 
                opts=opts, 
                **arguments,
            )
        
    # Update genome clusters
    organism_type_to_genome_to_cluster = defaultdict(dict)
    for organism_type in query_organism_types:
        with open_file_reader(os.path.join(opts.output_directory, "intermediate", organism_type, f"skani-dist_reference.tsv")) as f:
            logger.info(f"Reading skani output: {os.path.join(opts.output_directory, 'intermediate', f'skani-dist_reference.tsv')}")
            next(f)
            logger.info(f"Updating genome clusters using mode: {opts.af_mode}")
            if opts.af_mode == "relaxed":
                for line in tqdm(f, desc=f"Updating genome clusters using mode: {opts.af_mode}"):
                    line = line.strip()
                    if line:
                        id_reference_genome, id_query_genome, ani, af_ref, af_query, *_ = line.split("\t")
                        ani, af_ref, af_query = map(float, [ani, af_ref, af_query])
                        conditions = [
                            max(af_ref, af_query) >= opts.af_threshold, 
                            ani >= opts.ani_threshold,
                        ]
                        if all(conditions):
                            id_reference_genome, id_query_genome = map(get_filepath_basename, [id_reference_genome,id_query_genome])
                            id_genome_cluster = genome_to_cluster[id_reference_genome]
                            logger.info(f"Adding {id_query_genome}[AF: {af_query}] to genome cluster {id_genome_cluster} matching {id_reference_genome}[AF: {af_ref}] with ANI={ani}")
                            if id_query_genome in organism_type_to_genome_to_cluster[organism_type]:
                                current_cluster_label = organism_type_to_genome_to_cluster[organism_type][id_query_genome]
                                if current_cluster_label != id_genome_cluster:
                                    raise ValueError(f"Genome {id_query_genome} is in multiple clusters: {current_cluster_label} and {id_genome_cluster} which is likely a result of using different thresholds or modes during clustering")
                            organism_type_to_genome_to_cluster[organism_type][id_query_genome] = id_genome_cluster
            elif opts.af_mode == "strict":
                for line in tqdm(f, desc=f"Updating genome clusters using mode: {opts.af_mode}"):
                    line = line.strip()
                    if line:
                        id_reference_genome, id_query_genome, ani, af_ref, af_query, *_ = line.split("\t")
                        conditions = [
                            af_ref >=opts.af_threshold, 
                            af_query >= opts.af_threshold,
                            ani >= opts.ani_threshold,
                        ]
                        if all(conditions):
                            id_genome_cluster = genome_to_cluster[id_reference_genome]
                            logger.info(f"Adding {id_query_genome}[AF: {af_query}] to genome cluster {id_genome_cluster} matching {id_reference_genome}[AF: {af_ref}] with ANI={ani}")
                            id_genome_cluster = genome_to_cluster[id_reference_genome]
                            if id_query_genome in organism_type_to_genome_to_cluster[organism_type]:
                                current_cluster_label = organism_type_to_genome_to_cluster[organism_type][id_query_genome]
                                if current_cluster_label != id_genome_cluster:
                                    raise ValueError(f"Genome {id_query_genome} is in multiple clusters: {current_cluster_label} and {id_genome_cluster} which is likely a result of using different thresholds or modes during clustering")
                            organism_type_to_genome_to_cluster[organism_type][id_query_genome] = id_genome_cluster
                                
    # Write genome clusters
    output = list()
    for organism_type, genome_to_cluster in organism_type_to_genome_to_cluster.items():
        for id_genome, id_genome_cluster in genome_to_cluster.items():
            output.append([id_genome, id_genome_cluster])
    logger.info(f"Identified N={len(output)} genomes with cluster representatives")


    # Phase II   
    logger.info(f"Phase II: Assigning remaining genomes to new clusters")
    
    C=0
    for organism_type, query_genomes in query_organism_type_to_genomes.items():
        query_genomes_without_clusters = query_genomes - set(genome_to_cluster.keys())
        if query_genomes_without_clusters:
            logger.info(f"Assigning N={len(query_genomes_without_clusters)} {organism_type} genomes to new clusters")
            arguments = {
            "organism_type": organism_type,
            "working_directory":os.path.join(opts.output_directory, "intermediate", organism_type),
            "query_genomes": query_genomes_without_clusters,
            }
        
            if organism_type in {"prokaryotic", "eukaryotic"}:
                run_skani_triangle_nonviral_query_based(
                    logger=logger, 
                    log_directory=os.path.join(opts.output_directory, "logs"), 
                    opts=opts, 
                    **arguments,
                )
            if organism_type in {"viral"}:
                run_skani_triangle_viral_query_based(
                    logger=logger, 
                    log_directory=os.path.join(opts.output_directory, "logs"), 
                    opts=opts, 
                    **arguments,
                )
            # Compile clusters
            run_edgelist_to_clusters(
                logger=logger, 
                log_directory=os.path.join(opts.output_directory, "logs"), 
                opts=opts, 
                **arguments,
                )
            with open_file_reader(os.path.join(opts.output_directory, "intermediate", organism_type, "skani-triangle_query.clusters.tsv")) as f:
                logger.info(f"Reading query clustering output: {os.path.join(opts.output_directory, 'intermediate',organism_type, "skani-triangle_query.clusters.tsv")}")
                for line in tqdm(f, desc=f"Updating genome clusters using mode: {opts.af_mode}"):
                    line = line.strip()
                    if line:
                        id_query_genome, id_genome_cluster = line.split("\t")                                           
                        organism_type_to_genome_to_cluster[organism_type][id_query_genome] = id_genome_cluster
                        output.append([id_query_genome, id_genome_cluster])
                        C += 1
    logger.info(f"Assigned N={C} genomes to new clusters")
    if not opts.no_references:
        C = 0
        logger.info("Appended reference genome clusters to output")
        with open_file_reader(opts.reference_genomes_with_clusters) as f:
            for line in f:
                line = line.strip()
                if not line.startswith("#"):
                    organism_type, id_genome_reference, id_genome_cluster, path_genome = line.split("\t")
                    organism_type_to_genome_to_cluster[organism_type][id_genome_reference] = id_genome_cluster
                    output.append([id_genome_reference, id_genome_cluster])
                    C += 1
        logger.info(f"Appended N={C} reference genomes and clusters")

    
    df_output = pd.DataFrame(output, columns=["id_genome", "id_genome_cluster"])
    df_output = df_output.sort_values(ascending=[True, True], by=[ "id_genome_cluster", "id_genome"])
    df_output.to_csv(os.path.join(opts.output_directory, "output", "genome_to_cluster.tsv"), sep="\t", header=None, index=None)
            
    # ========
    # Hash
    # ========   
    md5hash_filepath = os.path.join(opts.output_directory, "output", "md5hashes.json")
    logger.info(f"Calculating md5 hashes: {md5hash_filepath}")
    file_to_hash = get_md5hash_from_directory(os.path.join(opts.output_directory, "intermediate"))
    write_json(file_to_hash, md5hash_filepath)

    # ========
    # Complete
    # ========    
    logger.info(f"Completed building updating genome clusters: {opts.output_directory}")
    logger.info(f"Directory size of output directory: {format_bytes(get_directory_size(opts.output_directory))}")
    logger.info("Updated genome clusters: {}".format(os.path.join(opts.output_directory, "output", "genome_to_cluster.tsv")))

if __name__ == "__main__":
    main()
    
    

    
