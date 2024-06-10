#!/usr/bin/env python
import sys, os, argparse, gzip, warnings
# from collections import OrderedDict
# import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.5.21"

    
def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <input.txt[.gz]> -o <output.tsv[.gz]".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog)
    # Pipeline
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-i","--input", type=str, default="stdin", help = "path/to/minpath_report.txt [Default: stdin]")
    parser_io.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv[.gz] formatted output [Default: stdout]")
    parser_io.add_argument("--no_header", action="store_true", help="Specify if header should be in output")

  
    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Input
    if opts.input == "stdin":
        f_input = sys.stdin
    else:
        if opts.input.endswith(".gz"):
            f_input = gzip.open(opts.input, "rt")
        else:
            f_input = open(opts.input, "r")
    
    # Output
    if opts.output == "stdout":
        f_output = sys.stdout
    else:
        if opts.output.endswith(".gz"):
            f_output = gzip.open(opts.output, "wt")
        else:
            f_output = open(opts.output, "w")

    # Parse
    if not opts.no_header:
        header = ["id_minpath", "mode", "reconstruction_available", "naive_reconstructed", "minpath_passed", "number_of_families_in_reference_pathway", "number_of_families_annotated", "name"]
        print(*header, sep="\t", file=f_output)

    for line in tqdm(f_input):
        line = line.strip()
        if line:
            left, right = line.split("  name  ")
            fields = left.split(" ")
            fields = list(filter(bool, fields))
            reconstruction_available = fields[3]
            
            row = [
                fields[1],
                fields[2],
                reconstruction_available if reconstruction_available != "n/a" else False,
                bool(eval(fields[5])),
                bool(eval(fields[7])),
                fields[9],
                fields[11],
                right,
            ]
            print(*row, sep="\t", file=f_output)

    # Closing
    if f_input != sys.stdin:
        f_input.close()
    if f_output != sys.stdout:
        f_output.close()

if __name__ == "__main__":
    main()
    
                

