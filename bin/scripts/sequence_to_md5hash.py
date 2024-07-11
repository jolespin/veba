#!/usr/bin/env python
import sys, os, argparse, warnings
import hashlib
# from collections import OrderedDict
# import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.6.11"
    
def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "grep v '^>' sequences.fasta | {} -f fasta > sequences-md5hash.fasta".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-f","--output_format", type=str, choices={"fasta", "list"}, default="fasta", help = "[Default: fasta]")
    parser.add_argument("-c","--case", type=str, choices={"upper", "lower", "insensitive"}, default="upper", help = "[Default: upper]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    if opts.output_format == "list":
        # Insensitive
        if opts.case == "insensitive":
            for line in tqdm(sys.stdin):
                line = line.strip()
                if line:
                    seq = line
                    id_md5 = hashlib.md5(seq.encode('utf-8')).hexdigest()
                    print(id_md5, file=sys.stdout)
        # Upper
        if opts.case == "upper":
            for line in tqdm(sys.stdin):
                line = line.strip()
                if line:
                    seq = line.upper()
                    id_md5 = hashlib.md5(seq.encode('utf-8')).hexdigest()
                    print(id_md5, file=sys.stdout)
             
        # Lower       
        if opts.case == "lower":
            for line in tqdm(sys.stdin):
                line = line.strip()
                if line:
                    seq = line.lower()
                    id_md5 = hashlib.md5(seq.encode('utf-8')).hexdigest()
                    print(id_md5, file=sys.stdout)
                    
    if opts.output_format == "fasta":
        # Insensitive
        if opts.case == "insensitive":
            for line in tqdm(sys.stdin):
                line = line.strip()
                if line:
                    seq = line
                    id_md5 = hashlib.md5(seq.encode('utf-8')).hexdigest()
                    print(">{}\n{}".format(id_md5, seq), file=sys.stdout)
        # Upper
        if opts.case == "upper":
            for line in tqdm(sys.stdin):
                line = line.strip()
                if line:
                    seq = line.upper()
                    id_md5 = hashlib.md5(seq.encode('utf-8')).hexdigest()
                    print(">{}\n{}".format(id_md5, seq), file=sys.stdout)
             
        # Lower       
        if opts.case == "lower":
            for line in tqdm(sys.stdin):
                line = line.strip()
                if line:
                    seq = line.lower()
                    id_md5 = hashlib.md5(seq.encode('utf-8')).hexdigest()
                    print(">{}\n{}".format(id_md5, seq), file=sys.stdout)
                
if __name__ == "__main__":
    main()
    
                

