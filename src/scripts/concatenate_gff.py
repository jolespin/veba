#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse
import pandas as pd

# Soothsayer Ecosystem
from soothsayer_utils import get_file_object, pv, assert_acceptable_arguments

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.02.17"


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__

    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <record_table> -o [id_sample]/concatenated.gff".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", default="stdin", type=str, help = "[id_sample]<tab>[path/to/reference.gff] [Default: stdin]")
    parser.add_argument("-o","--output_directory", default="concatenated_output", type=str, help = "Concatenated output directory [Default: concatenated_output]")
    parser.add_argument("-x", "--extension", type=str, default="gff", help="Concatenated fasta output extension [Default: gff] ")
    parser.add_argument("-b", "--basename", type=str, default="concatenated", help="Concatenated fasta output extension [Default: concatenated] ")
    parser.add_argument("--no_sort", action="store_true", help = "Don't sort the grouped filepaths")
    parser.add_argument("--no_subdirectory", action="store_true", help = "Don't create a nested directory structure")
    parser.add_argument("-M", "--mode", type=str, default="infer", help="Concatenate all references with global and build index or build index for each reference {global, local, infer}")

    # parser.add_argument("--no_ignore_comments", action="store_true", help = "Don't ignore comments")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Input 
    f_in = opts.input
    if opts.input == "stdin":
        f_in = sys.stdin 

    # Output
    os.makedirs(opts.output_directory, exist_ok=True)

    # Read table
    df = pd.read_csv(f_in, sep="\t", index_col=0, header=None)
    m = df.shape[1]

    if opts.mode == "infer":
        assert_acceptable_arguments(m, {0,1})
        opts.mode = {0:"global", 1:"local"}[m]
    
    if opts.mode == "global":
        assert m == 0, "There should be only one column if mode='global'"
    if opts.mode == "local":
        assert m == 1, "There should be two columns if mode='local'"
    
    assert_acceptable_arguments(opts.mode, {"local", "global"})

    if opts.mode == "local":

        # GroupBy and parse
        for id_sample, filepaths in df.groupby(df.index):
            # Get filepaths
            filepaths = list(filepaths.values.ravel())
            if not opts.no_sort:
                filepaths = sorted(filepaths)
            # Create output files (and directories)
            if opts.no_subdirectory:
                f_out = get_file_object(
                    path=os.path.join(opts.output_directory, "{}.{}".format(id_sample, opts.extension)),
                    mode="write", 
                    safe_mode=False, 
                    verbose=False,
                )
            else:
                os.makedirs(os.path.join(opts.output_directory, id_sample), exist_ok=True)

                f_out = get_file_object(
                    path=os.path.join(opts.output_directory, id_sample, "{}.{}".format(opts.basename, opts.extension)),
                    mode="write", 
                    safe_mode=False, 
                    verbose=False,
                )
            
            # Read input gff, filter out short sequences, and write to concatenated file
            for fp in pv(filepaths, description=id_sample, unit= " files"):
                f_query = get_file_object(fp, mode="read", verbose=False)
                for line in f_query:
                    line = line.strip()
                    if line:
                        if not line.startswith("#"):
                            print(line, file=f_out)

                f_query.close()

            f_out.close()


    if opts.mode == "global":
        # GroupBy and parse
        filepaths = df.index 
        if not opts.no_sort:
            filepaths = sorted(filepaths)

        # Create output files (and directories)
        f_out = get_file_object(
            path=os.path.join(opts.output_directory, "{}.{}".format(opts.basename, opts.extension)),
            mode="write", 
            safe_mode=False, 
            verbose=False,
        )

        # Read input fasta, filter out short sequences, and write to concatenated file
        for fp in pv(filepaths, unit= " files"):
            f_query = get_file_object(fp, mode="read", verbose=False)
            for line in f_query:
                line = line.strip()
                if line:
                    if not line.startswith("#"):
                        print(line, file=f_out)


            f_query.close()

        f_out.close()

if __name__ == "__main__":
    main()
