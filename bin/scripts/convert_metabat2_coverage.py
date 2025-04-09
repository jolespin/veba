#!/usr/bin/env python
from __future__ import print_function, division
import sys
import os
import argparse
import warnings
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2025.4.8"
        
def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "cat <input.fasta> | {} --prefix <prefix> --suffix <suffix> --minimum_length <int> > <output.fasta>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", type=str, default="stdin", help = "path/to/coverage_metabat2.tsv")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/coverage.tsv")
    parser.add_argument("-n","--index_name", type=str, default="contigname", help = "Name of first column [Default: contigname]")
    parser.add_argument("--no_header", action="store_true", help = "Do not include header")
    parser.add_argument("--identifiers", type=str,  help = "Contig identifier order and subset. Fill missing with 0.")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Input
    if opts.input == "stdin":
        opts.input = sys.stdin 
    df_input = pd.read_csv(opts.input, sep="\t", index_col=0)
    df_input = df_input.iloc[:,2:]
    mask = df_input.columns.map(lambda x: x.endswith((".bam",".sam")))
    df_input = df_input.loc[:,mask]
    df_input.index.name = opts.index_name

    # Identifiers
    if opts.identifiers:
        identifiers = list()
        with open(opts.identifiers) as f:
            for line in f:
                line = line.strip()
                if line:
                    identifiers.append(line)
        A = set(df_input.index)
        B = set(identifiers)
        if not A >= B:
            warnings.warn("Identifiers not found in input. Filling with 0. Missing: {}".format(A-B))
        df_input = df_input.reindex(identifiers)
        df_input = df_input.fillna(0)
        
    # Output
    if opts.output == "stdout":
        opts.output = sys.stdout 
        
    df_input.to_csv(opts.output, sep="\t", header=not bool(opts.no_header))

if __name__ == "__main__":
    main()
