#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse


# import pandas as pd
# import numpy as np

from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.07.31"

# .............................................................................
# Notes
# .............................................................................
# * Make batch version that takes in a manifest file
# .............................................................................
# Primordial
# .............................................................................
# # Check filename
# def check_filename(filename, acceptable_characters={".","-","_"}):
#     """
#     Source: https://github.com/jolespin/genopype
#     """
#     status_ok = True
#     for character in str(filename).lower():
#         conditions = [
#         character.isalnum(),
#         character in acceptable_characters,
#         # character != " ",
#         ]
#         if not any(conditions):
#             status_ok = False
#             break
#     return status_ok

# # Format filename
# def format_filename(name, replacement_character="_", acceptable_characters={".","-","_"}):
#     """
#     Source: https://github.com/jolespin/genopype
#     """
#     listed_string = list(name)
#     idx_nonalnum = list()
#     for i, character in enumerate(listed_string):
#         if not check_filename(character):
#             idx_nonalnum.append(i)
#     for i in idx_nonalnum:
#         listed_string[i] = replacement_character
#     return "".join(listed_string)


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "cat <input.fasta> | {} --prefix <prefix> --suffix <suffix> --minimum_length <int> > <output.fasta>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", type=str, default="stdin", help = "path/to/input.fasta")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.fasta")
    parser.add_argument("-p","--prefix", type=str,  help = "prefix for header")
    parser.add_argument("-s","--suffix", type=str,  help = "prefix for header")
    parser.add_argument("-m","--minimum_length", type=int,  help = "Include only sequences that are greater than or equal to this length")

    # parser.add_argument("--simplify", type=str,  help = "Simplify header")
    # parser.add_argument("--lengths", type=str,  help = "Simplify header")
    # parser.add_argument("--gc_content", type=str,  help = "Simplify header")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Input
    if opts.input == "stdin":
        f_in = sys.stdin 
    else:
        f_in = open(opts.input, "r")

    # Output
    if opts.output == "stdout":
        f_out = sys.stdout 
    else:
        f_out = open(opts.output, "w")

    # Defaults
    if opts.prefix is None:
        opts.prefix = ""
    if opts.suffix is None:
        opts.suffix = ""
    if opts.minimum_length is None:
        opts.minimum_length = 1

    # Read and write fasta
    for id, seq in tqdm(SimpleFastaParser(f_in), desc="Processing fasta file", unit=" records"):
        id = "{}{}{}".format(opts.prefix, id, opts.suffix).strip()
        if len(seq) >= opts.minimum_length:
            print(">{}\n{}".format(id, seq), file=f_out)

    # Close files
    if f_in is not sys.stdin:
        f_in.close()
    if f_out is not sys.stdout: 
        f_out.close()


if __name__ == "__main__":
    main()
