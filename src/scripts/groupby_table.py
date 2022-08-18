#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse
import pandas as pd

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.08.17"

#
def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -m <mapping.tsvt> -t <table.tsv> -o <output.tsv>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-t","--table", required="stdin", type=str, help = "path/to/table.tsv, No header. [Default: stdin]")
    parser.add_argument("-m","--mapping", required = True, type=str, help = "path/to/mapping.tsv [id_key]<tab>[id_group], No header")
    parser.add_argument("-o","--output_table", default="stdout", type=str, help = "path/to/output_table.tsv, [Default: stdout]")
    parser.add_argument("-a","--axis", type=int, default=0, help = "index:axis=0, columns:axis=1")
    parser.add_argument("--sep", type=str, default="\t", help = "Separator [Default: <tab>]")
    parser.add_argument("--skiprows", type=int, help = "Skiprows for --table")
    # parser.add_argument("--table_header", action="store_true",  help = "--table header")
    # parser.add_argument("--mapping_header", action="store_true",  help = "--mapping header")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert opts.axis in {0,1}, "--axis must be in {0,1}"

    # I/O
    if opts.table == "stdin":
        opts.table = sys.stdin 
   
    if opts.output_table == "stdout":
        opts.output_table = sys.stdout 

    # Read Table
    df_table = pd.read_csv(opts.table, sep=opts.sep, index_col=0, skiprows=opts.skiprows)#, header=bool(opts.table_header))
    key_to_group = pd.read_csv(opts.mapping, sep=opts.sep, index_col=0).iloc[:,0]#, header=bool(opts.mapping_header)).iloc[:,0]

    df_output = df_table.groupby(key_to_group, axis=opts.axis).sum()

    # Write table
    df_output.to_csv(opts.output_table, sep=opts.sep)#, header=bool(opts.table_header))

if __name__ == "__main__":
    main()
