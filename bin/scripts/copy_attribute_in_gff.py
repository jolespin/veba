#!/usr/bin/env python
import sys, os, glob, argparse 
from collections import OrderedDict
# import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.12.23"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <input.gff> -s ID -dgene_id -o <output.gff>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", type=str, default="stdin", help = "path/to/gene_models.gff [Default: stdin]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/gene_models.updated.gff [Default: stdout]")
    parser.add_argument("-s","--source_attribute", type=str, default="ID",  help = "Attribute to add to GFF [Default: gene_id]")
    parser.add_argument("-d","--destination_attribute", type=str, default="gene_id",  help = "Attribute to add to GFF [Default: gene_id]")
    parser.add_argument("-f","--feature", type=str,help = "Only include this feature type. Case sensitive.  Default behavior is to update all features.")
    parser.add_argument("-c", "--add_contig_attribute", action="store_true",  help = "Add contig_id to GFF")
    parser.add_argument("-e", "--exclude_fasta", action="store_true",  help = "If fasta is included in GFF then exclude it")

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
    # Feature filter
    if opts.feature:
        if opts.add_contig_attribute:
            for line in f_in:
                line = line.strip()
                if line:
                    if opts.exclude_fasta:
                        if line.startswith(("##FASTA", ">")):
                            break
                    if line.startswith("#") or "\t" not in line:
                        print(line, file=f_out)
                    else:
                        if not line.endswith(";"):
                            line += ";"
                        try:
                            id_contig, tool, feature, *_ = line.split("\t")
                        except ValueError:
                            raise ValueError(f"Error parsing line: {line}")
                        if all([
                            f"{opts.source_attribute}=" in line,
                            feature == opts.feature,
                            ]):
                            value = line.split(f"{opts.source_attribute}=")[1].split(";")[0].strip()
                            print(
                                line,
                                "contig_id",
                                "=",
                                id_contig,
                                ";",
                                opts.destination_attribute, 
                                "=",
                                value,
                                ";",
                                sep="",
                            file=f_out)
                        else:
                            print(line, file=f_out)
                    
                    
        else:
            for line in f_in:
                line = line.strip()
                if line:
                    if opts.exclude_fasta:
                        if line.startswith(("##FASTA", ">")):
                            break
                    if line.startswith("#") or "\t" not in line:
                        print(line, file=f_out)
                    else:
                        if not line.endswith(";"):
                            line += ";"
                        try:
                            id_contig, tool, feature, *_ = line.split("\t")
                        except ValueError:
                            raise ValueError(f"Error parsing line: {line}")
                        if all([
                            f"{opts.source_attribute}=" in line,
                            feature == opts.feature,
                            ]):
                            value = line.split(f"{opts.source_attribute}=")[1].split(";")[0].strip()
                            print(
                                line,
                                opts.destination_attribute, 
                                "=",
                                value,
                                ";",
                                sep="",
                            file=f_out)
                        else:
                            print(line, file=f_out)
                    
    # No feature filter   
    else:
        if opts.add_contig_attribute:
            for line in f_in:
                line = line.strip()
                if line:
                    if opts.exclude_fasta:
                        if line.startswith(("##FASTA", ">")):
                            break
                    if line.startswith("#") or "\t" not in line:
                        print(line, file=f_out)
                    else:
                        if not line.endswith(";"):
                            line += ";"
                        try:
                            id_contig, tool, feature, *_ = line.split("\t")
                        except ValueError:
                            raise ValueError(f"Error parsing line: {line}")
                        if f"{opts.source_attribute}=" in line:
                            value = line.split(f"{opts.source_attribute}=")[1].split(";")[0].strip()
                            print(
                                line,
                                "contig_id",
                                "=",
                                id_contig,
                                ";",
                                opts.destination_attribute, 
                                "=",
                                value,
                                ";",
                                sep="",
                            file=f_out)
                        else:
                            print(line, file=f_out)
                    

        for line in f_in:
            line = line.strip()
            if line:
                if opts.exclude_fasta:
                    if line.startswith(("##FASTA", ">")):
                        break
                if line.startswith("#") or "\t" not in line:
                    print(line, file=f_out)
                else:
                    if not line.endswith(";"):
                        line += ";"
                    try:
                        id_contig, tool, feature, *_ = line.split("\t")
                    except ValueError:
                        raise ValueError(f"Error parsing line: {line}")
                    if f"{opts.source_attribute}=" in line:
                        value = line.split(f"{opts.source_attribute}=")[1].split(";")[0].strip()
                        print(
                            line,
                            opts.destination_attribute, 
                            "=",
                            value,
                            ";",
                            sep="",
                        file=f_out)
                    else:
                        print(line, file=f_out)


    # Close files
    if f_in is not sys.stdin:
        f_in.close()
    if f_out is not sys.stdout: 
        f_out.close()


if __name__ == "__main__":
    main()
    