#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse
import pandas as pd

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.04.20"

#
def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <index.list> -t <table.tsv> -o <output.tsv>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--index", default="stdin", type=str, help = "path/to/index.list [Default: stdin]")
    parser.add_argument("-t","--table", required=True, type=str, help = "path/to/table.tsv")
    parser.add_argument("-o","--output_table", default="stdout", type=str, help = "path/to/output_table.tsv [Default: stdout]")
    parser.add_argument("-a","--axis", type=int, default=0, help = "index:axis=0, columns:axis=1")
    parser.add_argument("--sep", type=str, default="\t", help = "Separator [Default: <tab>]")
    parser.add_argument("--skiprows", type=int, help = "Skiprows")
    parser.add_argument("-v", "--inverse", action="store_true", help = "Inverse")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert opts.axis in {0,1}, "--axis must be in {0,1}"

    # I/O
    if opts.index == "stdin":
        f_index = sys.stdin 
    else:
        f_index = open(opts.index, "r")
    index = list() 
    for line in f_index:
        line = line.strip()
        if line:
            index.append(line)
    if not opts.index == "stdin":
        f_index.close() 

    if opts.output_table == "stdout":
        opts.output_table = sys.stdout 

    # Read Table
    df = pd.read_csv(opts.table, sep=opts.sep, index_col=0, skiprows=opts.skiprows)
    if not opts.inverse:
        if opts.axis == 0:
            assert set(index) <= set(df.index), "--index isn't a subset of --table index"
            df = df.loc[index]
        if opts.axis == 1:
            assert set(index) <= set(df.columns), "--index isn't a subset of --table columns"
            df = df.loc[:,index]
    else:
        if opts.axis == 0:
            assert set(index) <= set(df.index), "--index isn't a subset of --table index"
            df = df.drop(index, axis=0)
        if opts.axis == 1:
            assert set(index) <= set(df.columns), "--index isn't a subset of --table columns"
            df = df.drop(index, axis=1)

    # Write table
    df.to_csv(opts.output_table, sep=opts.sep)

if __name__ == "__main__":
    main()
