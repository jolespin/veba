#!/usr/bin/env python
import sys, os, glob, argparse 
from collections import OrderedDict
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.08.03"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <tblout> -f <extension> -o <output.tsv>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--tblout", required=True, type=str, help = "path/to/hmmer/tblout")
    parser.add_argument("-f","--fasta",  type=str, default="stdin", help = "path/to/protein/fasta/file [Default: stdin]")
    parser.add_argument("-o","--output",  type=str, default="stdout", help = "path/to/output.faa [Default: stdout]")
    parser.add_argument("--identifiers_only", action="store_true", help="Output identifiers only")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Input
    if opts.fasta == "stdin":
        f_in = sys.stdin 
    else:
        f_in = open(opts.fasta, "r")

    # Output
    if opts.output == "stdout":
        f_out = sys.stdout 
    else:
        f_out = open(opts.output, "w")

    if not opts.identifiers_only:
        assert opts.fasta is not None, "If not --identifiers_only then must specify --fasta"
    
    identifiers = set()
    with open(opts.tblout, "r") as f:
        for line in f.readlines():
            line = line.strip()
            if not line.startswith("#"):
                fields = line.split(" ")
                id = fields[0]
                identifiers.add(id)
    identifiers = sorted(identifiers)

    if opts.identifiers_only:
        for id in identifiers:
            print(id, file=f_out)
    else:
        from Bio.SeqIO.FastaIO import SimpleFastaParser
        for (header, seq) in SimpleFastaParser(f_in):
            id = header.split(" ")[0]
            if id in identifiers:
                print(">{}\n{}".format(header, seq), file=f_out)

    if f_out is not sys.stdout: 
        f_out.close()



if __name__ == "__main__":
    main()
    
                

