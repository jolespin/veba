#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse
import pandas as pd

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.2.15"

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
    parser.add_argument('input', type=str, nargs='?', default="stdin", help='input table [Default: stdin]')
    parser.add_argument("-o","--output", default="stdout", type=str, help = "path/to/output_table.tsv [Default: stdout]")
    parser.add_argument("-f", "--columns", type=str, help = "Column positions separated by commas")
    parser.add_argument("-d", "--sep", type=str, default="\t", help = "Separator [Default: <tab>]")
    parser.add_argument("-i", "--index_column", type=int, default=0, help = "Index column [Default: 0]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # I/O
    if opts.input == "stdin":
        opts.input = sys.stdin 

    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Read Table
    df = pd.read_csv(opts.input, sep=opts.sep, index_col=opts.index_column)

    # Subset
    if opts.columns:
        columns = list()
        for field in opts.columns.strip().split(","):
            columns.append(field)
        assert set(columns) <= set(df.columns), "The following fields are not columns in the table:{}".format("\n".join(list(set(columns) - set(df.columns))))
        df = df.loc[:,columns]

    # Write table
    df.to_csv(opts.output, sep=opts.sep)

if __name__ == "__main__":
    main()
