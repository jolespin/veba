#!/usr/bin/env python
import sys, os, glob, argparse, gzip
from collections import OrderedDict
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.09.30"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <hmmsearch_tblout> -n <name> -a <proteins> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--hmmsearch_tblout", type=str, help = "path/to/hmmsearch_tblout.tsv", required=True)
    parser.add_argument("-a","--proteins", type=str, help = "path/to/proteins.faa",  required=True)
    parser.add_argument("-o","--output_directory", type=str,  default="input", help = "Output directory [Default: Directory of input hmmsearch tblout]")
    parser.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    # parser.add_argument("-p", "--identifier_prefix", type=str,  help="Identifier prefix [Default: '")
    # parser.add_argument("-s", "--identifier_suffix", type=str,  help="Identifier suffix [Default: '")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    if opts.output_directory == "input":
        opts.output_directory = os.path.split(opts.hmmsearch_tblout)[0]
    
    os.makedirs(opts.output_directory, exist_ok=True)

    # Get marker to protein
    marker_to_evalue = dict()
    marker_to_protein = dict()
    with open(opts.hmmsearch_tblout, "r") as f:
        for line in f.readlines():
            line = line.strip()
            if not line.startswith("#"):
                fields = list(filter(bool, line.split(" ")))
                id_query = fields[0]
                id_target = fields[3]
                evalue = eval(fields[4])
                if id_target not in marker_to_protein:
                    print(opts.name, "|", id_target, "-->", id_query, "| E-Value = {:.1e}".format(evalue), file=sys.stderr)
                    marker_to_protein[id_target] = id_query
                    marker_to_evalue[id_target] = evalue
                else:
                    if evalue < marker_to_evalue[id_target]:
                        print(opts.name, "|", id_target, "-->", id_query, "| E-Value = {:.1e}".format(evalue), "[Update] [Overriding: {}]".format(marker_to_protein[id_target]), file=sys.stderr)
                        marker_to_protein[id_target] = id_query
                        marker_to_evalue[id_target] = evalue
    
    # Protein sequences 
    protein_to_seq = dict() 
    marker_proteins = set(marker_to_protein.values())
    with open(opts.proteins, "r") as f:
        for header, seq in SimpleFastaParser(f):
            id_protein = header.split(" ")[0]
            if id_protein in marker_proteins:
                protein_to_seq[id_protein] = seq 

    # Write sequences 
    for id_marker, id_protein in marker_to_protein.items():
        with gzip.open(os.path.join(opts.output_directory, "{}.faa.gz".format(id_marker)), "wt") as f:
            seq = protein_to_seq[id_protein]
            # header = id_protein
            # if opts.identifier_prefix:
            #     header = "{}{}".format(opts.identifier_prefix, header)
            # if opts.identifier_suffix:
            #     header = "{}{}".format(header, opts.identifier_suffix) 
            # header = "{} {}".format(header, opts.hmmsearch_tblout)
            header = "{}  {} {} {}".format(opts.name, id_marker, id_protein, opts.hmmsearch_tblout)
            print(">{}\n{}".format(header,seq), file=f)

if __name__ == "__main__":
    main()
    
                

