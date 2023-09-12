#!/usr/bin/env python
import sys, os, glob, argparse, warnings
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.8.30"

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
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")
    parser.add_argument("-a","--axis", type=int,  default=0, help = "axis to merge dataframes [Default: 0]")
    parser.add_argument("-i", "--index_column", type=str,  default="0", help = "Comma-separeted list for index (e.g., 0,1) [Default: 0]")
    parser.add_argument("-c", "--header", type=str,  default="0", help = "Comma-separeted list for header (e.g., 0,1) [Default: 0]")
    parser.add_argument("-n", "--no_header", action="store_true", help = "No header")
    parser.add_argument("-d","--delimiter", type=str,  default="\t", help = "delimter [Default: <tab>]")
    parser.add_argument("-e", "--allow_empty_or_missing_files", action="store_true", help = "Allow empty or missing files")
    parser.add_argument("--prepend_index_levels", type=str,  help = "Comma-separeted list to prepend to index of each dataframe.")
    parser.add_argument("--prepend_column_levels", type=str,  help = "Comma-separeted list to prepend to header of each dataframe.  Must match the number of dataframes exactly.")
    parser.add_argument("--sort_by", type=str,  help = "Sort by values.  This happens before (to remove duplicates) and after concatenation.")
    parser.add_argument("--ascending", action="store_true", help = "Ascending order instead of descending")

    # Options
    opts = parser.parse_args(argv)
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert opts.axis in {0,1}

    # Open output file
    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Format header
    if opts.no_header:
        opts.header = None
    else:
        opts.header = list(map(int, opts.header.split(",")))

    # Format index
    opts.index_column = list(map(int, opts.index_column.split(",")))
    assert all(list(map(lambda x: x >= 0, opts.index_column))), "--index_column must contain values ≥ 0"

    # Format column levels
    if opts.prepend_column_levels:
        assert not opts.no_header, "Cannot use --prepend_column_levels with --no_header"
        opts.prepend_column_levels = list(map(lambda x:x.strip(), opts.prepend_column_levels.split(",")))
        assert len(opts.prepend_column_levels) == len(opts.dataframes), "Number of levels in --prepend_column_levels must match number of dataframes"
    
    # Format index levels
    if opts.prepend_index_levels:
        opts.prepend_index_levels = list(map(lambda x:x.strip(), opts.prepend_index_levels.split(",")))    

    if opts.allow_empty_or_missing_files:
        dataframes = list()
        for i, fp in enumerate(opts.dataframes):
            df = None
            try:
                df = pd.read_csv(fp, sep=opts.delimiter, index_col=opts.index_column, header=opts.header)

                # Handle duplicates
                if opts.sort_by:
                    df = df.sort_values(by=opts.sort_by, ascending=bool(opts.ascending))
                number_of_duplicates = len(df.index.value_counts()[lambda x: x > 1].index)
                if number_of_duplicates > 0:
                    warnings.warn(f"{fp} contains {number_of_duplicates} duplicates.  Removing them now and keeping first instance.")
                    df = df[~df.index.duplicated(keep='first')]

                if opts.prepend_column_levels is not None:
                    df.columns = df.columns.map(lambda x: (opts.prepend_column_levels[i], x))
                if opts.prepend_index_levels is not None:
                    df.index = df.index.map(lambda x: (*opts.prepend_index_levels, x))
            except (pd.errors.EmptyDataError, FileNotFoundError) as e:
                print("[Skipping] {}".format(e), file=sys.stderr)
            if df is not None:
                dataframes.append(df)
    else:
        dataframes = list()
        for i, fp in enumerate(opts.dataframes):
            df = pd.read_csv(fp, sep=opts.delimiter, index_col=opts.index_column, header=opts.header)

            # Handle duplicates
            if opts.sort_by:
                df = df.sort_values(by=opts.sort_by, ascending=bool(opts.ascending))
            number_of_duplicates = len(df.index.value_counts()[lambda x: x > 1].index)
            if number_of_duplicates > 0:
                warnings.warn(f"{fp} contains {number_of_duplicates} duplicates.  Removing them now and keeping first instance.")
                df = df[~df.index.duplicated(keep='first')]

            if opts.prepend_column_levels is not None:
                df.columns = df.columns.map(lambda x: (opts.prepend_column_levels[i], x))
            if opts.prepend_index_levels is not None:
                df.index = df.index.map(lambda x: (*opts.prepend_index_levels, x))
            dataframes.append(df)
        
    df_concat = pd.concat(dataframes, axis=opts.axis)

    # Output
    df_concat.to_csv(opts.output, sep=opts.delimiter, header=not opts.no_header)


if __name__ == "__main__":
    main()
    
                

