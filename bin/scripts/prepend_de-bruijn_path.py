#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, gzip
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.12.11"

        
def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "cat <scaffolds.path> | {} --prefix <prefix>  > <scaffolds.prefixed.path>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jol.espinoz@gmail.com)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", type=str, default="stdin", help = "path/to/input.fasta")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.fasta")
    parser.add_argument("-p","--prefix", type=str, required=True, help = "prefix for header")
    # parser.add_argument("-s","--suffix", type=str,  help = "prefix for header")
    parser.add_argument("-P","--program", type=str,  choices={"spades", "flye"}, required=True,help = "Assembly program used to build paths")


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
        f_out = sys.stdout 
    else:
        if opts.output.endswith(".gz"):
            f_out = gzip.open(opts.output, "wt")
        else:
            f_out = open(opts.output, "w")

    # Defaults
    if opts.prefix is None:
        opts.prefix = ""
    # if opts.suffix is None:
    #     opts.suffix = ""


    # Main loop
    if opts.program == "spades":
        for line in tqdm(f_in, desc=f"Prepending {opts.prefix} to contig headers in de Bruijn graph paths"):
            line = line.strip()
            if line:
                if line.startswith("NODE_"):
                    print(f"{opts.prefix}{line}", file=f_out)
                else:
                    print(line, file=f_out)
                    
    if opts.program == "flye":
        for line in tqdm(f_in, desc=f"Prepending {opts.prefix} to contig headers in de Bruijn graph paths"):
            line = line.strip()
            if line:
                if "contig_" in line:
                    line = line.replace("contig_", f"{opts.prefix}contig_")
                    print(line, file=f_out)
                else:
                    print(line, file=f_out)

    # Close files
    if f_in is not sys.stdin:
        f_in.close()
    if f_out is not sys.stdout: 
        f_out.close()


if __name__ == "__main__":
    main()
