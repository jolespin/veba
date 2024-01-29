#!/usr/bin/env python
import sys, os, glob, argparse, warnings
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.03.24"

def main(argv=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -t <table> -c <column> -d <delimiter> -o <output.tsv> --header_table 0 --header_column 0".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument('-t', '--table', type=str, default="stdin",   help='path/to/table [Default: stdin]')
    parser.add_argument('-c', '--column_values', type=str,  required=True, help='path/to/column')
    parser.add_argument('-n', '--column_name', type=str,  required=True, help='Insert column name')
    parser.add_argument('-i', '--column_index', type=int,  required=True, help='Insert column position')
    parser.add_argument("--header_table", type=int,  help = "Header. Choose [Default: None]")
    parser.add_argument("--header_column", type=int,  help = "Header. Choose [Default: None]")

    parser.add_argument("-d","--delimiter", type=str,  default="\t", help = "delimter [Default: <tab>]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")

    # Options
    opts = parser.parse_args(argv)
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Input
    if opts.table == "stdin":
        opts.table = sys.stdin

    # Open output file
    if opts.output == "stdout":
        opts.output = sys.stdout 

    df = pd.read_csv(opts.table, sep=opts.delimiter, index_col=0, header=opts.header_table )

    new_column = pd.read_csv(opts.column_values, sep=opts.delimiter, index_col=0, header=opts.header_column ).iloc[:,0]


    # Output
    df.insert(loc=opts.column_index, column=opts.column_name, value=new_column)
    df.to_csv(opts.output, sep=opts.delimiter, header=opts.header_table)


if __name__ == "__main__":
    main()
    
                

