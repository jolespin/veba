#!/usr/bin/env python
import sys, os, argparse, gzip, pickle, hashlib
from collections import defaultdict
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
from Bio.SeqIO.FastaIO import SimpleFastaParser
from pyexeggutor import open_file_writer


__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2025.4.8"

def generate_unique_cluster_base_label(nodes:set, mode:str, index:int, cluster_prefix_zfill:int):
    ALPHANUMERIC=list("1234567890abcdefghijklmnopqrstuvwxyz")
    # {"numeric", "md5", "pseudo-random", "random"}
    if mode == "numeric":
        return str(index).zfill(cluster_prefix_zfill)
    elif mode == "nodes":
        return "|".join(sorted(nodes))
    elif mode == "md5":
        return hashlib.md5(str(sorted(nodes)).encode()).hexdigest()
    elif mode == "pseudo-random":
        return "".join(np.random.RandomState(seed=index).choice(ALPHANUMERIC, 32))
    elif mode == "random":
        return "".join(np.random.choice(ALPHANUMERIC, 32))

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <edgelist.tsv [No header]> -o <output.tsv>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline

    parser.add_argument("-i","--input", type=str, default="stdin", help = \
                        "path/to/edgelist.tsv, No header. Accepted formats:\n"
                        "[id_1]<tab>[id_2]\n"
                        "[id_1]<tab>[id_2]<tab>[weight]\n"
                        "[id_1]<tab>[id_2]<tab>[weight_1]<tab>[weight_2]\n"
                        "[id_1]<tab>[id_2]<tab>[weight]<tab>[alignment_fraction_reference]<tab>[alignment_fraction_query]\n"
                        "[Default: stdin]") #  
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/clusters.tsv [Default: stdout]")
    parser.add_argument("-t","--threshold", type=float, default=0.0,  help = "Minimum weight threshold(1). [Default: 0.0]") 
    parser.add_argument("--threshold2", type=float, default=0.0,  help = "Minimum weight threshold2. [Default: 0.0]") 
    parser.add_argument("-a","--minimum_af", type=float, default=0.0,  help = "Minimum alignment fraction. [Default: 0.0]") 
    parser.add_argument("-m","--af_mode", type=str, default="relaxed",  choices={"relaxed", "strict"}, help = "Minimum alignment fraction mode with either `relaxed = max([AF_ref, AF_query]) > minimum_af` or `strict = (AF_ref > minimum_af) & (AF_query > minimum_af)`. `strict` will be biased against fragmented or partial genomes likely derived from metagenomes [Default: relaxed]") 
    parser.add_argument("-n", "--no_singletons", action="store_true", help = "Don't include self-interactions. Self-interactions will ensure unclustered genomes make it into the output")
    parser.add_argument("-b", "--basename", action="store_true", help = "Removes filepath prefix and extension.  Support for gzipped filepaths.")
    parser.add_argument("--identifiers", type=str, help = "Identifiers to include.  If missing identifiers and singletons are allowed, then they will be included as singleton clusters with weight of np.nan")

    parser_labels = parser.add_argument_group('Label arguments')
    parser_labels.add_argument("-p", "--cluster_prefix", type=str, default="c-", help="Cluster prefix [Default: 'c-']")
    parser_labels.add_argument("-z", "--cluster_prefix_zfill", type=int, default=0, help="Cluster prefix zfill. Use 7 to match identifiers from OrthoFinder.  Use 0 to add no zfill. Only applicable when --cluster_label_mode numeric. [Default: 0]") #7
    parser_labels.add_argument("-s", "--cluster_suffix", type=str, default="", help="Cluster suffix [Default: '']")
    parser_labels.add_argument("-c", "--cluster_label_mode", type=str, default="md5", choices={"numeric", "random", "pseudo-random", "md5", "nodes"}, help="Cluster label. [Default: 'md5']")

    parser_fasta = parser.add_argument_group('Fasta arguments')
    parser_fasta.add_argument("-f","--fasta", type=str,  help = "path/to/sequences.fasta")
    parser_fasta.add_argument("--output_fasta_directory", type=str, default="clusters", help = "path/to/clusters/cluster_x.fasta [Default: clusters]")
    parser_fasta.add_argument("-x", "--output_fasta_extension", type=str, default="fasta", help = "path/to/clusters/cluster_x.[extension] [Default: fasta]")

    parser_export = parser.add_argument_group('Export arguments')
    parser_export.add_argument("-g", "--export_graph", type=str,   help = "prefix/to/graph pickled output files: nx.Graph suggested prefix is .graph.pkl")
    parser_export.add_argument("-d", "--export_dict", type=str,   help = "prefix/to/dict pickled output file: {id_cluster: {id_a, id_b, ...}} suggested prefix is .dict.pkl")
    parser_export.add_argument("-r", "--export_representatives", type=str,   help = "prefix/to/representatives.tsv table")

    # parser.add_argument("--export_edgelist", type=str,   help = "prefix/to/edgelist output files")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    
    # Input
    if opts.input == "stdin":
        opts.input = sys.stdin 

    # Output
    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Edge list
    try:
        df_edgelist = pd.read_csv(opts.input, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        df_edgelist = pd.DataFrame(columns=["query", "reference"])

    assert df_edgelist.shape[1] in  {2,3,4, 5}, "Must have 2, 3, 4, or 5 columns.  {} provided.".format(df_edgelist.shape[1])
    if opts.basename:
        def get_basename(x):
            _, fn = os.path.split(x)
            if fn.endswith(".gz"):
                fn = fn[:-3]
            return ".".join(fn.split(".")[:-1])
        
        # Support for .map superceding .applymap
        if hasattr(df_edgelist, "map"):
            df_edgelist.iloc[:,:2] = df_edgelist.iloc[:,:2].map(get_basename)
        else:
            df_edgelist.iloc[:,:2] = df_edgelist.iloc[:,:2].applymap(get_basename)

    # Identifiers from edgelist
    if not df_edgelist.empty:
        edgelist = df_edgelist.iloc[:,:2].values.tolist()
        identifiers = set.union(*map(set, edgelist))
    else:
        edgelist = list()
        identifiers = set()

    all_identifiers = identifiers
    if opts.identifiers:
        all_identifiers = set()
        with open(opts.identifiers, "r") as f:
            for line in f.readlines():
                id = line.strip()
                all_identifiers.add(id)

    # Read in fasta
    if opts.fasta:
        id_to_sequence = dict()
        if opts.fasta.endswith(".gz"):
            f = gzip.open(opts.fasta, "rt")
        else:
            f = open(opts.fasta, "r")
        for header, seq in tqdm(SimpleFastaParser(f), "Reading fasta file: {}".format(opts.fasta)):
            id = header.split(" ")[0]
            id_to_sequence[id] = seq 

        f.close()
        assert set(id_to_sequence.keys()) >= identifiers, "Not all of the sequences in --input are available in --fasta.  Either add the sequences to --fasta file or remove --fasta argument."
        os.makedirs(opts.output_fasta_directory, exist_ok=True)

    # Construct graph
    graph = nx.Graph()

    # Unweighted
    if df_edgelist.shape[1] == 2:
        for i, (id_query, id_target) in tqdm(df_edgelist.iterrows(), "Reading edgelist: {}".format(opts.input), total=df_edgelist.shape[0]):
            if {id_query, id_target}.issubset(all_identifiers):
                graph.add_edge(id_query, id_target, weight=1)
        if not opts.no_singletons:
            for id in all_identifiers:
                if id not in graph.nodes():
                    graph.add_edge(id, id, weight=np.nan)

    # Weighted
    if df_edgelist.shape[1] == 3:
        for i, (id_query, id_target, w) in tqdm(df_edgelist.iterrows(), "Reading edgelist with weights (≥ {}): {}".format(opts.threshold, opts.input), total=df_edgelist.shape[0]):
            if w >= opts.threshold:
                if {id_query, id_target}.issubset(all_identifiers):
                    graph.add_edge(id_query, id_target, weight=w)
        if not opts.no_singletons:
            for id in all_identifiers:
                if id not in graph.nodes():
                    graph.add_edge(id, id, weight=np.nan)

    # Weighted with alignment fraction
    if df_edgelist.shape[1] == 4:
        for i, (id_query, id_target, w1, w2) in tqdm(df_edgelist.iterrows(), "Reading edgelist with weights (≥ {}) and weights2 (≥ {}): {}".format(opts.threshold, opts.threshold2, opts.input), total=df_edgelist.shape[0]):
            if all([
                w1 >= opts.threshold,
                w2 >= opts.threshold2,
                ]):
                if {id_query, id_target}.issubset(all_identifiers):
                    graph.add_edge(id_query, id_target, weight=w1, weight2=w2)
        if not opts.no_singletons:
            for id in all_identifiers:
                if id not in graph.nodes():
                    graph.add_edge(id, id, weight=np.nan, weight2=np.nan)
                    
    # Weighted with alignment fraction
    if df_edgelist.shape[1] == 5:
        if opts.af_mode == "relaxed":
            for i, (id_query, id_target, w, af_ref, af_query) in tqdm(df_edgelist.iterrows(), "Reading edgelist with weights (≥ {}) and alignment fractions (≥ {}) using {} mode: {}".format(opts.threshold, opts.minimum_af, opts.af_mode, opts.input), total=df_edgelist.shape[0]):
                if all([
                    w >= opts.threshold,
                    max([af_ref, af_query]) >= opts.minimum_af,
                    ]):
                    if {id_query, id_target}.issubset(all_identifiers):
                        graph.add_edge(id_query, id_target, weight=w, alignment_fraction_reference=af_ref, alignment_fraction_query=af_query)
                        
        if opts.af_mode == "strict":
            for i, (id_query, id_target, w, af_ref, af_query) in tqdm(df_edgelist.iterrows(), "Reading edgelist with weights (≥ {}) and alignment fractions (≥ {}) using {} mode: {}".format(opts.threshold, opts.minimum_af, opts.af_mode, opts.input), total=df_edgelist.shape[0]):
                if all([
                    w >= opts.threshold,
                    af_ref >= opts.minimum_af,
                    af_query >= opts.minimum_af,
                    ]):
                    if {id_query, id_target}.issubset(all_identifiers):
                        graph.add_edge(id_query, id_target, weight=w, alignment_fraction_reference=af_ref, alignment_fraction_query=af_query)
                        
        if not opts.no_singletons:
            for id in all_identifiers:
                if id not in graph.nodes():
                    graph.add_edge(id, id, weight=np.nan, alignment_fraction=100.0)
            
    # Get connected components
    node_to_cluster = dict()
    cluster_to_nodes = dict()

    for i, nodes in tqdm(enumerate(sorted(nx.connected_components(graph), key=len, reverse=True), start=1), "Organizing clusters", unit=" clusters"):
        id_cluster = generate_unique_cluster_base_label(nodes=nodes, mode=opts.cluster_label_mode, index=i, cluster_prefix_zfill=opts.cluster_prefix_zfill)

        # Add cluster prefix and suffix
        if bool(opts.cluster_prefix):
            id_cluster =  "{}{}".format(opts.cluster_prefix, id_cluster)
        if bool(opts.cluster_suffix):
            id_cluster =  "{}{}".format(id_cluster, opts.cluster_suffix)
                    
        # Get subgraph
        if len(nodes) > 1:
            subgraph = graph.subgraph(nodes)
            node_to_degree = dict(nx.degree(subgraph, weight="weight"))
            representative = max(node_to_degree, key=node_to_degree.get)
        else:
            node_to_degree = dict(zip(nodes, [np.nan]))
            representative = list(nodes)[0]

        for id_node in nodes:
            node_to_cluster[id_node] = id_cluster
            graph.nodes[id_node]["id_cluster"] = id_cluster
            graph.nodes[id_node]["intra-cluster_connectivity"] = node_to_degree[id_node]
            graph.nodes[id_node]["representative"] = id_node == representative

        cluster_to_nodes[id_cluster] = set(nodes)
    node_to_cluster = pd.Series(node_to_cluster, name="Clusters")
    node_to_cluster.to_frame().to_csv(opts.output, sep="\t", header=None)

    # Export pickle
    if opts.export_graph is not None:
        with open_file_writer("{}".format(opts.export_graph)) as f:
            pickle.dump(graph, f)

    if opts.export_dict is not None:
        with open_file_writer("{}".format(opts.export_dict)) as f:
            pickle.dump(cluster_to_nodes, f)

    # Export representatives
    if opts.export_representatives is not None:
        with open_file_writer("{}".format(opts.export_representatives)) as f:
            print("id_node", "id_cluster", "intra-cluster_connectivity", "representative", sep="\t", file=f)
            for id_node, node_metadata in graph.nodes(data=True):
                k = node_metadata["intra-cluster_connectivity"]
                print(
                    id_node, 
                    node_metadata["id_cluster"], 
                    k if pd.notnull(k) else "", 
                    node_metadata["representative"], 
                    sep="\t", 
                    file=f,
                )


    # # Export edgelist
    # if opts.export_edgelist is not None:
    #     nx.write_weighted_edgelist(graph, "{}.edgelist.tsv".format(opts.export_pickle), delimiter="\t")

    if opts.fasta:
        for id_cluster, nodes in tqdm(cluster_to_nodes.items(), "Writing fasta files for each cluster: {}".format(opts.output_fasta_directory)):
            with open(os.path.join(opts.output_fasta_directory, "{}.{}".format(id_cluster, opts.output_fasta_extension)), "w") as f:
                for id_node in sorted(nodes):
                    print(">{} {}\n{}".format(id_node, id_cluster, id_to_sequence[id_node]), file=f)

if __name__ == "__main__":
    main()
    
                