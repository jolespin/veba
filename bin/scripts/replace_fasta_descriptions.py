#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, gzip
from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.11.05"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -f <input.fasta> -d <descriptions.tsv -o <output.fasta>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-f","--fasta", required=True, type=str, help = "path/to/fasta[.gz]. Use --fasta stdin for piping")
    parser.add_argument("-d","--descriptions", required=True, type=str, help = "Either 1) a string; or 2) path/to/descriptions. Tab separated table with header [id, field_1, ..., field_2]")
    parser.add_argument("-o","--output", default="stdout", type=str, help = "Output filepath [Default: stdout]")
    parser.add_argument("-a","--allow_missing_descriptions", action="store_true", help = "Allow for descriptions for identifiers that are not in fasta and vice versa")
    parser.add_argument("-n","--no_header", action="store_true", help = "No header in description")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename



    # Input
    if opts.fasta == "stdin":
        f_in = sys.stdin
    else:
        if opts.fasta.endswith(".gz"):
            f_in = gzip.open(opts.fasta, "rt")
        else:
            f_in = open(opts.fasta, "r")
    
    # Output
    if opts.output == "stdout":
        f_out = sys.stdout 
    else:
        f_out = open(opts.output, "w")


    # Table descriptions
    if os.path.exists(opts.descriptions):
        if opts.no_header:
            df_descriptions = pd.read_csv(opts.descriptions, sep="\t", index_col=0, header=None)
            id_to_description = df_descriptions.apply(lambda x: ";".join(x))
        else:
            df_descriptions = pd.read_csv(opts.descriptions, sep="\t", index_col=0)
            id_to_description = dict()
            for id, descriptions in df_descriptions.iterrows():
                fields = list()
                for k,v in descriptions.items():
                    fields.append("[{}={}]".format(k,v))
                id_to_description[id] = "|".join(fields)
        id_to_description = pd.Series(id_to_description)

        # Append descriptions to fasta
        if opts.allow_missing_descriptions:
            for id, seq in tqdm(SimpleFastaParser(f_in), "Reading fasta file: {}".format(opts.fasta)):
                id = id.split(" ")[0]
                if id in id_to_description:
                    description = id_to_description[id]
                    print(">{} {}\n{}".format(id, description, seq), file=f_out)
                else:
                    print(">{}\n{}".format(id, seq), file=f_out)
        else:
            for id, seq in tqdm(SimpleFastaParser(f_in), "Reading fasta file: {}".format(opts.fasta)):
                id = id.split(" ")[0]
                assert id in id_to_description, "Fasta record identifier not in --descriptions table"
                description = id_to_description[id]
                print(">{} {}\n{}".format(id, description, seq), file=f_out)


    # String descriptions
    else:
        for id, seq in tqdm(SimpleFastaParser(f_in), "Reading fasta file: {}".format(opts.fasta)):
            id = id.split(" ")[0]
            description = opts.descriptions
            print(">{} {}\n{}".format(id, description, seq), file=f_out)
    
    # Closing files
    if f_in != sys.stdin:
        f_in.close()
    if f_out != sys.stdout:
        f_out.close()

if __name__ == "__main__":
    main()
