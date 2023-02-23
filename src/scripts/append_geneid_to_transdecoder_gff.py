#!/usr/bin/env python
import sys, os, glob, argparse 
from collections import OrderedDict
# import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.2.22"

transcript_to_gene = dict()
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
    # parser.add_argument("-a","--attribute", type=str, default="gene_id",  help = "Attribute to add to GFF [Default: gene_id]")


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
        if line:
            if line.startswith("#"):
                print(line, file=f_out)
            else:
                if "Name" in line:
                    line = line.split(";Name")[0]
                

                fields = line.split("\t")

                id_contig = fields[0]

                if fields[2] not in {"gene", "exon", "CDS", "mRNA"}:
                    print(line + ";", file=f_out)
                else:

                    # Update the descriptors
                    if fields[2] != "exon":

                        # Get transcript to gene mapping
                        if fields[2] == "gene":
                            g,t = fields[-1][3:].split(";")[0].split("~~")
                            transcript_to_gene[t] = g

                        updated_descriptors = list()
                        if fields[2] in {"gene", "mRNA"}:
                            for x in fields[-1].split(";"):
                                if "~~" in x:
                                    x = x.split("~~")[0]
                                updated_descriptors.append(x)
                        if fields[2] == "CDS":
                            for x in fields[-1].split(";"):
                                if x.startswith("ID="):
                                    t = x[7:]
                                    x = "ID={}".format(t)
                                updated_descriptors.append(x)
                        fields[-1] = ";".join(updated_descriptors)

                    
                    # Get gene identifier
                    if fields[2] == "gene":
                        id_gene = fields[-1].split("ID=")[1].split(";")[0]
                        print(
                            "\t".join(fields),
                            # contig
                            ";",
                            "contig_id",
                            "=",
                            id_contig,

                            # gene
                            ";",
                            "gene_id",
                            "=",
                            id_gene,
                            ";",

                            sep="",
                        file=f_out)

                    if fields[2] in {"mRNA", "CDS"}:
                        # print(fields)
                        id_transcript = fields[-1].split("ID=")[1].split(";")[0]
                        id_gene = transcript_to_gene[id_transcript]

                        print(
                            "\t".join(fields),
                            # contig
                            ";",
                            "contig_id",
                            "=",
                            id_contig,

                            # gene
                            ";",
                            "gene_id",
                            "=",
                            id_gene,

                            # transcript
                            ";",
                            "transcript_id",
                            "=",
                            id_transcript,
                            ";",

                            sep="",
                        file=f_out)

                    if fields[2] == "exon":
                        id_exon = fields[-1].split(";")[0][3:]
                        id_transcript = id_exon.split(".exon")[0]
                        id_gene = transcript_to_gene[id_transcript]

                        print(
                            "\t".join(fields),
                            # contig
                            ";",
                            "contig_id",
                            "=",
                            id_contig,

                            # gene
                            ";",
                            "gene_id",
                            "=",
                            id_gene,

                            # transcript
                            ";",
                            "transcript_id",
                            "=",
                            id_transcript,

                            # exon
                            ";",
                            "exon_id",
                            "=",
                            id_exon,
                            ";",

                            sep="",
                        file=f_out)                        


    # # Close files
    if f_in is not sys.stdin:
        f_in.close()
    if f_out is not sys.stdout: 
        f_out.close()


if __name__ == "__main__":
    main()
    
                
