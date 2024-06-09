#!/usr/bin/env python
import sys, os, glob, argparse 
from collections import OrderedDict
# import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.6.29"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <input.gff> -a gene_id -o <output.gff>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", type=str, default="stdin", help = "path/to/gene_models.gff [Default: stdin]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/gene_models.updated.gff [Default: stdout]")
    parser.add_argument("-a","--attribute", type=str, default="gene_id",  help = "Attribute to add to GFF [Default: gene_id]")


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
 
    # Update GFF
    for line in f_in:
        line = line.strip()
        if line.startswith("#"):
            print(line, file=f_out)
        else:
            id_contig = line.split("\t")[0]
            id = line.split("ID=")[1].split(";")[0]
            print(
                line,
                "contig_id",
                "=",
                id_contig,
                ";",
                opts.attribute, 
                "=",
                id_contig,
                "_",
                id.split("_")[-1],
                ";",
                "gene_biotype=protein_coding",
                ";",
                sep="",
            file=f_out)

    # Close files
    if f_in is not sys.stdin:
        f_in.close()
    if f_out is not sys.stdout: 
        f_out.close()


if __name__ == "__main__":
    main()
    
                
#  dict(zip(faa.index.map(lambda x: x.split("ID=")[1].split(";")[0]), faa.index.map(lambda x: x.split(" ")[0])))
