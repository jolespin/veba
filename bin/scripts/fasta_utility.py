#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm
from typing import TextIO

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.11.9"

# Genomics
def fasta_writer(header:str, seq:str, file:TextIO, wrap:int=1000):
    # Write the FASTA header
    print(f">{header}", file=file)
    
    if wrap:
        # Write the sequence with lines of length 'wrap'
        for i in range(0, len(seq), wrap):
            # Get a chunk of the sequence with a max length of 'wrap'
            line = seq[i:i+wrap]
            print(line, file=file)
    else:
        print(seq, file=file)
        
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
    parser.add_argument("-w", "--wrap", type=int, default=1000, help = f"Line width for wrapping fasta lines.  If 0, then no wrapping is used. [Default: 1000]")

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
    if opts.suffix is None:
        opts.suffix = ""
    if opts.minimum_length is None:
        opts.minimum_length = 1
    if opts.wrap == 0:
        opts.wrap = None

    # Read and write fasta
    for header, seq in tqdm(SimpleFastaParser(f_in), desc="Processing fasta file", unit=" records"):
        id, *description = header.split(" ")
        description = " ".join(description)
        id = "{}{}{}".format(opts.prefix, id, opts.suffix).strip()
        if len(seq) >= opts.minimum_length:
            fasta_writer(
                header="{} {}".format(id, description),
                seq=seq,
                file=f_out,
                wrap=opts.wrap,
            )

    # Close files
    if f_in is not sys.stdin:
        f_in.close()
    if f_out is not sys.stdout: 
        f_out.close()


if __name__ == "__main__":
    main()
