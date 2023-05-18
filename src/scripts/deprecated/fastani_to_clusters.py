#!/usr/bin/env python
import sys, os, argparse, operator
from collections import defaultdict
import pandas as pd
import networkx as nx

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.1.20"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <input.gff> -a gene_id -o <output.gff>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", type=str, default="stdin", help = "path/to/fastani.tsv [Default: stdin]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/clusters.tsv [Default: stdout]")
    parser.add_argument("-t","--threshold", type=str, default="95",  help = "FastANI threshold.  Can be either scalar or comma-separated list (e.g., 95,96.5) [Default: 95]")
    parser.add_argument("-g","--genomes", type=str,  help = "path/to/genomes.list. If filepath, then the basename is taken")

    parser.add_argument("-p", "--cluster_prefix", type=str, default="SLC", help="Cluster prefix [Default: 'SLC")
    parser.add_argument("-s", "--cluster_suffix", type=str, default="", help="Cluster suffix [Default: '")

    parser.add_argument("-x","--fasta_extension", type=str, default="fa",  help = "Fasta extension [Default: 'fa']")
    parser.add_argument("--export_pickle", type=str,   help = "prefix/to/pickle output files")
    parser.add_argument("--export_edgelist", type=str,   help = "prefix/to/edgelist output files")
    # parser.add_argument("-m", "--metadata", type=str, help="Metadata to append to table")
    parser.add_argument("--include", type=str,  help = "path/to/list of identifiers without header (Note, this happens after name trimming)")
    parser.add_argument("--interval_type", type=str, default="open",  help = "Interval type either 'open' (>) or 'closed' (>=)  [Default: open]")
    parser.add_argument("--no_self", action="store_true", help = "Don't include self-interactions. Self-interactions will ensure unclustered genomes make it into the output")
    parser.add_argument("--no_trim", action="store_true", help = "Don't trim the name.  Trimming includes removing prefix paths and file extension")
    parser.add_argument("--no_header", action="store_true", help = "Don't add header to output.  Not recommended when running multiple ANIs")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Input
    if opts.input == "stdin":
        f_in = sys.stdin 
    else:
        f_in = open(opts.input, "r")

    # Output
    if opts.output == "stdout":
        f_out = sys.stdout 
    else:
        f_out = open(opts.output, "w")
    
    # Pickles        

    assert opts.interval_type in {"closed", "open"}, "--interval_type must be either 'closed' or 'open'"
    if opts.interval_type == "open":
        f_interval = operator.gt
    else:
        f_interval = operator.gte

    opts.threshold = list(map(lambda x: float(x.strip()), opts.threshold.split(",")))
    df_fastani = pd.read_csv(opts.input, sep="\t", header=None)

    if not opts.no_trim:
        f_trim = lambda x:x.split("/")[-1][:-1*(len(opts.fasta_extension) + 1)]
        for j in [0,1]:
            df_fastani.iloc[:,j] = df_fastani.iloc[:,j].map(f_trim)

    ani_values = df_fastani.set_index([0,1]).iloc[:,0]
    identifiers = set.union(*ani_values.index.map(set))

    if opts.genomes:
        genome_set = set()
        with open(opts.genomes, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    if "/" in line:
                        id_genome = line.split("/")[-1]
                        if id_genome.endswith(".{}".format(opts.fasta_extension)):
                            id_genome = id_genome[:-1*(len(opts.fasta_extension) + 1)]
                            genome_set.add(id_genome)
        identifiers |= genome_set

    inclusion_set = identifiers
    if opts.include:
        inclusion_set = set()
        with open(opts.include, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    inclusion_set.add(line)
        assert inclusion_set <= identifiers, "There are {} identifiers that are not in --input: {}".format(len(inclusion_set - identifiers))
    
    connected_components = defaultdict(dict)
    for tol in opts.threshold:
        tol_label = "Cluster({})".format(tol)
        # Construct graph
        graph = nx.Graph()
        for (id_query, id_target), ani in ani_values.items():
            if all([
                id_query in inclusion_set,
                id_target in inclusion_set,
            ]):
                if f_interval(ani, tol):
                    graph.add_edge(id_query, id_target, weight=ani)
        if not opts.no_self:
            for id_genome in inclusion_set:
                graph.add_edge(id_genome, id_genome, weight=100.0)
                
        # Get connected components
        for id_cluster, nodes in enumerate(sorted(nx.connected_components(graph), key=len, reverse=True), start=1):
            # Add cluster prefix and suffix
            if bool(opts.cluster_prefix):
                id_cluster =  "{}{}".format(opts.cluster_prefix, id_cluster)
            if bool(opts.cluster_suffix):
                id_cluster =  "{}{}".format( id_cluster, opts.cluster_suffix)
            for node in nodes:
                connected_components[tol_label][node] = id_cluster

        # Export pickle
        if opts.export_pickle is not None:
            nx.write_gpickle(graph, "{}-ani_{}.graph.pkl".format(opts.export_pickle, tol))
        # Export edgelist
        if opts.export_edgelist is not None:
            nx.write_weighted_edgelist(graph, "{}-ani_{}.edgelist.tsv".format(opts.export_pickle, tol), delimiter="\t")


        
    connected_components = pd.DataFrame(connected_components)
    connected_components = connected_components.sort_values(list(connected_components.columns))
    if opts.no_header:
        connected_components.to_csv(f_out, sep="\t", header=None)
    else:
        connected_components.to_csv(f_out, sep="\t")

    # Close files
    if f_in is not sys.stdin:
        f_in.close()
    if f_out is not sys.stdout: 
        f_out.close()


if __name__ == "__main__":
    main()
    
                
#  dict(zip(faa.index.map(lambda x: x.split("ID=")[1].split(";")[0]), faa.index.map(lambda x: x.split(" ")[0])))
