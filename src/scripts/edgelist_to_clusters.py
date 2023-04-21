#!/usr/bin/env python
import sys, os, argparse, gzip, pickle
from collections import defaultdict
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
from Bio.SeqIO.FastaIO import SimpleFastaParser

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.4.17"

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

    parser.add_argument("-i","--input", type=str, default="stdin", help = "path/to/edgelist.tsv, No header. [id_1]<tab>[id_2] or [id_1]<tab>[id_2]<tab>[weight] [Default: stdin]") #  
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/clusters.tsv [Default: stdout]")
    parser.add_argument("-t","--threshold", type=float, default=0.5,  help = "Minimum weight threshold. [Default: 0.5]") 
    parser.add_argument("-n", "--no_singletons", action="store_true", help = "Don't include self-interactions. Self-interactions will ensure unclustered genomes make it into the output")
    parser.add_argument("-b", "--basename", action="store_true", help = "Removes filepath prefix and extension.  Support for gzipped filepaths.")
    parser.add_argument("--identifiers", type=str, help = "Identifiers to include.  If missing identifiers and singletons are allowed, then they will be included as singleton clusters with weight of np.inf")

    parser_labels = parser.add_argument_group('Label arguments')
    parser_labels.add_argument("-p", "--cluster_prefix", type=str, default="c", help="Cluster prefix [Default: 'c']")
    parser_labels.add_argument("-z", "--cluster_prefix_zfill", type=int, default=0, help="Cluster prefix zfill. Use 7 to match identifiers from OrthoFinder.  Use 0 to add no zfill. [Default: 0]") #7
    parser_labels.add_argument("-s", "--cluster_suffix", type=str, default="", help="Cluster suffix [Default: '']")

    parser_fasta = parser.add_argument_group('Fasta arguments')
    parser_fasta.add_argument("-f","--fasta", type=str,  help = "path/to/sequences.fasta")
    parser_fasta.add_argument("--output_fasta_directory", type=str, default="clusters", help = "path/to/clusters/cluster_x.fasta [Default: clusters]")
    parser_fasta.add_argument("-x", "--output_fasta_extension", type=str, default="fasta", help = "path/to/clusters/cluster_x.[extension] [Default: fasta]")


    parser_export = parser.add_argument_group('Export arguments')
    parser_export.add_argument("-g", "--export_graph", type=str,   help = "prefix/to/graph pickled output files: nx.Graph suggested prefix is .graph.pkl")
    parser_export.add_argument("-d", "--export_dict", type=str,   help = "prefix/to/dict pickled output file: {id_cluster: {id_a, id_b, ...}} suggested prefix is .dict.pkl")

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
    df_edgelist = pd.read_csv(opts.input, sep="\t", header=None)
    assert df_edgelist.shape[1] in  {2,3}, "Must have 2 or 3 columns.  {} provided.".format(df_edgelist.shape[1])
    if opts.basename:
        def get_basename(x):
            _, fn = os.path.split(x)
            if fn.endswith(".gz"):
                fn = fn[:-3]
            return ".".join(fn.split(".")[:-1])
        df_edgelist.iloc[:,:2] = df_edgelist.iloc[:,:2].applymap(get_basename)

    edgelist = df_edgelist.iloc[:,:2].values.tolist()

    identifiers = set.union(*map(set, edgelist))

    all_identifiers = identifiers
    if opts.identifiers:
        all_identifiers = set()
        with open(opts.identifiers, "r") as f:
            for line in f.readlines():
                id = line.strip()
                all_identifiers.add(id)

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
            graph.add_edge(id_query, id_target)
        if not opts.no_singletons:
            for id in all_identifiers:
                if id not in graph.nodes():
                    graph.add_edge(id, id)

    # Weighted
    if df_edgelist.shape[1] == 3:
        for i, (id_query, id_target, w) in tqdm(df_edgelist.iterrows(), "Reading edgelist: {}".format(opts.input), total=df_edgelist.shape[0]):
            if w >= opts.threshold:
                graph.add_edge(id_query, id_target, weight=w)
        if not opts.no_singletons:
            for id in all_identifiers:
                if id not in graph.nodes():
                    graph.add_edge(id, id, weight=np.inf)
            
    # Get connected components
    node_to_cluster = dict()
    cluster_to_nodes = dict()

    for id_cluster, nodes in tqdm(enumerate(sorted(nx.connected_components(graph), key=len, reverse=True), start=1), "Organizing clusters"):
        # Add cluster prefix and suffix
        if bool(opts.cluster_prefix):
            id_cluster =  "{}{}".format(opts.cluster_prefix, str(id_cluster).zfill(opts.cluster_prefix_zfill))
        if bool(opts.cluster_suffix):
            id_cluster =  "{}{}".format(id_cluster, opts.cluster_suffix)
        for id_node in nodes:
            node_to_cluster[id_node] = id_cluster
            graph.nodes[id_node]["id_cluster"] = id_cluster
        cluster_to_nodes[id_cluster] = set(nodes)
    node_to_cluster = pd.Series(node_to_cluster, name="Clusters")
    node_to_cluster.to_frame().to_csv(opts.output, sep="\t", header=None)

    # Export pickle
    if opts.export_graph is not None:
        with open("{}".format(opts.export_graph), "wb") as f:
            pickle.dump(graph, f)

    if opts.export_dict is not None:
        with open("{}".format(opts.export_dict), "wb") as f:
            pickle.dump(cluster_to_nodes, f)

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
    
                