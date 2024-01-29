#!/usr/bin/env python
import sys, os, glob, argparse, warnings
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.1.31"

def main(argv=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <intput.tsv> -a <axis> -d <delimiter> -o <output.tsv> --how any".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", type=str, default="stdin", help = "path/to/input.tsv [Default: stdin]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")
    parser.add_argument("-a","--axis", type=int,  default=0, help = "axis to remove rows(0) or columns(1) [Default: 0]")
    parser.add_argument("-m", "--how", default="any", type=str, help = "{any,all [Default: any]")
    parser.add_argument("-n", "--no_header", action="store_true", help = "No header")
    parser.add_argument("-d","--delimiter", type=str,  default="\t", help = "delimter [Default: <tab>]")
    parser.add_argument("-I", "--index_column", type=str,  default=0, help = "Index column.  Use -1 for None [Default: 0]")

    # Options
    opts = parser.parse_args(argv)
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert opts.axis in {0,1}
    assert opts.how in {"any","all"}
    assert opts.index_column >= -1

    # Open input file
    if opts.input == "stdin":
        opts.input = sys.stdin 

    # Open output file
    if opts.output == "stdout":
        opts.output = sys.stdout    

    if opts.index_column == -1:
        opts.index_column = None

    # Format header
    df = pd.read_csv(opts.input, sep=opts.delimiter, index_col=opts.index_column, header=None if opts.no_header else 0)

    df = df.dropna(how=opts.how, axis=opts.axis)

    df.to_csv(opts.output, sep=opts.delimiter)


if __name__ == "__main__":
    main()
    
                

