#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse
import pandas as pd

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.2.9"

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
    parser.add_argument("-f", "--columns", type=str, help = "Column positions separated by commas [Integers]")
    parser.add_argument("-d", "--sep", type=str, default="\t", help = "Separator [Default: <tab>]")
    parser.add_argument("-0", "--python_indexing", action="store_true", help = "Python indexing which starts at 0.  Negative values can be used here which take from the end.")

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
    df = pd.read_csv(opts.input, sep=opts.sep, index_col=None, header=None)

    # Subset
    if opts.columns:
        columns = list()
        if opts.python_indexing:
            for i in map(int, opts.columns.strip().split(",")):
                columns.append(i)
        else:
            for i in map(int, opts.columns.strip().split(",")):
                assert i > 0, "All columns must be ≥ 1 if --start_at_1 is selected"
                columns.append(i - 1)
        df = df.iloc[:,columns]

    # Write table
    df.to_csv(opts.output, sep=opts.sep, index=None, header=None)

if __name__ == "__main__":
    main()
