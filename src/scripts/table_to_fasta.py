#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse
from collections import OrderedDict
import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.4.19"

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
    parser.add_argument("-i","--input", default="stdin", type=str, help = "path/to/table.tsv [Default: stdin]")
    parser.add_argument("-o","--output", default="stdout", type=str, help = "path/to/output.fasta [Default: stdout]")
    parser.add_argument("-n","--name_index", default=0, type=int, help = "Name index position [Default: 0]")
    parser.add_argument("-d","--description_index", type=int, help = "Description index position")
    parser.add_argument("-s","--sequence_index", default=2, type=int, help = "Sequence index position [Default: 2]")
    parser.add_argument("--no_header",action="store_true", help = "Use if there is no header")
    parser.add_argument("--sep", type=str, default="\t", help = "Separator [Default: <tab>]")
    parser.add_argument("--skiprows", type=int, help = "Skiprows")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename


    # I/O
    if opts.input == "stdin":
        opts.input = sys.stdin 

    if opts.output == "stdout":
        opts.output = sys.stdout
    else:
        opts.output = open(opts.output, "w")

    # Read Table
    if opts.no_header:
        df = pd.read_csv(opts.input, sep=opts.sep, index_col=None, skiprows=opts.skiprows, header=None)
    else:
        df = pd.read_csv(opts.input, sep=opts.sep, index_col=None, skiprows=opts.skiprows)#, header=bool(opts.no_header))

    if opts.description_index is not None:
        assert max(opts.name_index, opts.description_index, opts.sequence_index) < df.shape[1]
        for i, (id, description, seq) in tqdm(df.iloc[:,[opts.name_index, opts.description_index, opts.sequence_index]].iterrows(), "Writing sequence records", total=df.shape[0]):
            print(">{} {}\n{}".format(id, description, seq), file=opts.output)
    else:
        assert max(opts.name_index,  opts.sequence_index) < df.shape[1]
        for i, (id, seq) in tqdm(df.iloc[:,[opts.name_index, opts.sequence_index]].iterrows(), "Writing sequence records", total=df.shape[0]):
            print(">{}\n{}".format(id, seq), file=opts.output)

    if opts.output != sys.stdout:
        opts.output.close()
if __name__ == "__main__":
    main()
