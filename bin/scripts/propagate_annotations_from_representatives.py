#!/usr/bin/env python
import sys, os, glob, argparse, gzip
from collections import OrderedDict
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.5.14"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <input> -c <proteins_to_orthogroups.tsv> -o <output.tsv>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", type=str, default="stdin", help = "path/to/annotations.tsv.gz [Default: stdin]")
    parser.add_argument("-c","--protein_clusters", type=str, required=True, help = "Tab-seperated value table of [id_protein]<tab>[id_protein_cluster].  Use this if the --proteins are representative sequences")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/annotations.proteins.tsv.gz [Default: stdout]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    if opts.input == "stdin":
        opts.input = sys.stdin 

    if opts.output == "stdout":
        opts.output = sys.stdout

    # Input
    df_annotations_clusters = pd.read_csv(opts.input, sep="\t", header=[0,1], index_col=0)

    protein_to_cluster = pd.read_csv(opts.protein_clusters, sep="\t", index_col=0, header=None).iloc[:,0]

    # Output
    df_annotations_proteins = df_annotations_clusters.reindex(protein_to_cluster.values)
    df_annotations_proteins.index = protein_to_cluster.index
    df_annotations_proteins.insert(loc=0, column=("Identifiers", "id_protein_cluster"), value=protein_to_cluster)
    df_annotations_proteins.index.name = "id_protein"
    df_annotations_proteins.to_csv(opts.output, sep="\t")

if __name__ == "__main__":
    main()
    
                

