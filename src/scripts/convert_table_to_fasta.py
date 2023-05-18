#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, warnings, gzip 
from tqdm import tqdm 

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.5.17"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <table.tsv[.gz]> -o <output.fasta[.gz]>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", default="stdin", type=str, help = "path/to/input.tsv[.gz] 3 columns [id]<tab>[description]<tab>[sequence]")
    parser.add_argument("-o","--output",  default="stdout", type=str, help = "path/to/output.fasta[.gz]")
    parser.add_argument("--header",  action="store_true",  help = "Contains header")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Input
    if opts.input == "stdin":
        opts.input = sys.stdin
    else:
        if opts.input.endswith(".gz"):
            opts.input = gzip.open(opts.input, "rt")
        else:
            opts.input = open(opts.input, "r")
    if opts.header:
        next(opts.input)

    # Output
    if opts.output == "stdout":
        opts.output = sys.stdout
    else:
        if opts.output.endswith(".gz"):
            opts.output = gzip.open(opts.output, "wt")
        else:
            opts.output = open(opts.output, "w")

    # Convert
    for line in tqdm(opts.input, "Converting table to fasta"):
        line = line.strip()
        id, description, sequence = line.split("\t")
        print(">{} {}\n{}".format(id, description, sequence), file=opts.output)
    
    # Close
    if opts.input != sys.stdin:
        opts.input.close()
    if opts.output != sys.stdout:
        opts.output.close()

if __name__ == "__main__":
    main()
