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
    usage = "{} -a <axis> -d <delimiter> -o <output.tsv> [table_1] [table_2]".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument('dataframes', type=str, nargs='+', help='dataframes')
    parser.add_argument("-a","--axis", type=int,  default=0, help = "axis to merge dataframes")
    parser.add_argument("--header_index", type=str,  default="0", help = "Comma-separeted list for header (e.g., 0,1) [Default: 0]")
    parser.add_argument("-d","--delimiter", type=str,  default="\t", help = "delimter [Default: <tab>]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")
    parser.add_argument("-e", "--allow_empty_or_missing_files", action="store_true", help = "Allow empty or missing files")

    # Options
    opts = parser.parse_args(argv)
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert opts.axis in {0,1}

    # Open output file
    if opts.output == "stdout":
        opts.output = sys.stdout 

    opts.header_index = list(map(int, opts.header_index.split(",")))

    dataframes = list() 
    for fp in opts.dataframes:
        df = None
        try:
            df = pd.read_csv(fp, sep=opts.delimiter, index_col=0, header=opts.header_index)
        except (pd.errors.EmptyDataError, FileNotFoundError) as e:
            print("[Skipping] {}".format(e), file=sys.stderr)
        if df is not None:
            dataframes.append(df)

        
    df_concat = pd.concat(dataframes, axis=opts.axis)

    # Output
    df_concat.to_csv(opts.output, sep=opts.delimiter)


if __name__ == "__main__":
    main()
    
                

