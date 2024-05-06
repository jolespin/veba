#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, gzip
from collections import defaultdict
import pandas as pd

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.3.13"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <preprocess_directory> -b cleaned > <output_table>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline

    parser_io = parser.add_argument_group('I/O')
    parser_io.add_argument("-i","--input",  type=str, default="stdin", help = "STAR log file (i.e., Log.final.out[.gz]) [Default: stdin]")
    parser_io.add_argument("-o","--output", default="stdout", type=str, help = "Output filepath [Default: stdout]")
    parser_io.add_argument("-n","--name", type=str, help = "Name of sample")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Input
    if opts.input == "stdin":
        f_in = sys.stdin 
    else:
        if opts.input.endswith(".gz"):
            f_in = gzip.open(opts.input, "rt")
        else:
            f_in = open(opts.input, "r")

    # Output
    if opts.output == "stdout":
        opts.output = sys.stdout

    """
    STAR: https://github.com/alexdobin/STAR
    FILE: star_Log.final.out
    Designed using STAR v2.7.3a
    """
    data = dict()
    # category = "UTILITY"
    for line in f_in:
        line = line.strip()
        if line:
            if "|" not in line:
                category = line[:-1]
            else:
                field, info = map(lambda x:x.strip(), line.split("|"))
            if isinstance(info, str):
                try:
                    info = eval(info.replace("%",""))
                except (SyntaxError, NameError):
                    pass
            data[field] = info
    output = pd.Series(data)
    df = output.to_frame(opts.name)
    df.index.name = "id_metric"
    df.to_csv(opts.output, sep="\t", header=bool(opts.name))

if __name__ == "__main__":
    main()
