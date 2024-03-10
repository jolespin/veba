#!/usr/bin/env python
import sys, os, argparse, gzip 
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.11.10"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__

    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <input.fasta> -o <output.fasta>)".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-c","--genome_to_cluster", type=str, help = "Input fasta file")
    parser.add_argument("-g","--graph", type=str, help = "Input fasta file")
    parser.add_argument("-o","--output", default="stdout", type=str, help = "Output fasta file")
    parser.add_argument("-m","--maximum_weight", default=100, type=str, help = "Output fasta file")
    parser.add_argument("--genome_statistics",  type=str, help = "Output fasta file")
    parser.add_argument("--sort_by",  type=str, help = "Output fasta file")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename




G = nx.path_graph(4)  # or DiGraph, MultiGraph, MultiDiGraph, etc
H = G.subgraph([0, 1, 2])
list(H.edges)
[(0, 1), (1, 2)]





if __name__ == "__main__":
    main()
    
                

