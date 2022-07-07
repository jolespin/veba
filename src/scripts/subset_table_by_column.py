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
    usage = "{} -s <set.list> -t <table.tsv> -o <output.tsv> -c <column_name>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-s","--set", default="stdin", required=True, type=str, help = "path/to/index.list [Default: stdin]")
    parser.add_argument("-t","--table", required=True, type=str, help = "path/to/table.tsv")
    parser.add_argument("-o","--output_table", default="stdout", type=str, help = "path/to/output_table.tsv [Default: stdout]")
    parser.add_argument("-c", "--column", type=str, help = "Column to look for index")
    parser.add_argument("-i", "--icolumn", type=str, help = "Column position to look for index")
    parser.add_argument("--sep", type=str, default="\t", help = "Separator [Default: <tab>]")
    parser.add_argument("--skiprows", type=int, help = "Skiprows")
    parser.add_argument("-v", "--inverse", action="store_true", help = "Inverse")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert any([opts.column, opts.icolumn]), "Must provide either --column or --icolumn"
    assert bool(opts.column) != bool(opts.icolumn), "Cannot provide both --column and --icolumn"

    # I/O
    if opts.index == "stdin":
        f_set = sys.stdin 
    else:
        f_set = open(opts.set, "r")

    query_set = set() 
    for line in f_index:
        line = line.strip()
        if line:
            query_set.add(line)
    if not opts.index == "stdin":
        f_set.close() 
    
    if opts.output_table == "stdout":
        opts.output_table = sys.stdout 

    # Read Table
    df_input = pd.read_csv(opts.table, sep=opts.sep, index_col=0, skiprows=opts.skiprows)

    # Mask function
    func_mask = {True:lambda x: x not in query_set, False:lambda x: x in query_set}[bool(opts.inverse)]
    
    if opts.icolumn:
        mask = df_input.iloc[:,opts.icolumn].map(func_mask)
    else:
        mask = df_input.loc[:,opts.column].map(func_mask)

    df_output = df_input.loc[mask]


    # Write table
    df_output.to_csv(opts.output_table, sep=opts.sep)

if __name__ == "__main__":
    main()
