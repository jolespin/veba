#!/usr/bin/env python
import sys, os, glob, argparse, gzip
from collections import defaultdict
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.11.07"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i  <scaffolds_to_bins.tsv> -g <gene_models.gff> -d <gene_models.ffn> -a <gene_models.faa> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--scaffolds_to_bins", type=str, required=True,  help = "path/to/scaffolds_to_bins.tsv, Format: [id_scaffold]<tab>[id_bin], No header")
    parser.add_argument("-f","--fasta", type=str, required=False,  help = "path/to/scaffolds.fasta")
    parser.add_argument("-g","--gff", type=str, required=True, help = "path/to/gene_models.gff")
    parser.add_argument("-d","--cds", type=str,  required=False, help = "path/to/gene_models.ffn")
    parser.add_argument("-a","--protein", type=str, required=False,  help = "path/to/gene_models.faa")
    parser.add_argument("-o","--output_directory", type=str, default="gene_models_partitioned",  help = "Output directory for gene models [Default: gene_models_partitioned/]")
    parser.add_argument("--separate", action="store_true", help = "Separate bins into their own subdirectories")
    parser.add_argument("--include_unbinned", action="store_true", help = "Include unbinned contigs from identifier_mapping.tsv")
    parser.add_argument("--include_header", action="store_true", help = "Include headers on output tables")
    parser.add_argument("-M", "--use_mag_as_description", action="store_true", help = "Include MAG identifier for each contig in fasta header description")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Import fasta parser
    if opts.cds or opts.protein:
        from Bio.SeqIO.FastaIO import SimpleFastaParser

    # Parse scaffolds to bins
    scaffold_to_bin = dict()
    bins = set()
    print("Parsing scaffolds to bins file: {}".format(opts.scaffolds_to_bins), file=sys.stderr)
    with open(opts.scaffolds_to_bins) as f_scaffoldstobins:
        for line in f_scaffoldstobins.readlines():
            line = line.strip()
            if line:
                id_scaffold, id_bin = line.split("\t")
                scaffold_to_bin[id_scaffold] = id_bin
                bins.add(id_bin)

    # Make directories and file objects
    os.makedirs(opts.output_directory, exist_ok=True)

    file_objects = defaultdict(dict)
    if opts.separate:
        for id_bin in bins:
            os.makedirs(os.path.join(opts.output_directory,id_bin), exist_ok=True)
            file_objects[id_bin]["gff"] = open(os.path.join(opts.output_directory,id_bin, "gene_models.gff"),"w")

            if opts.cds:
                file_objects[id_bin]["cds"] = open(os.path.join(opts.output_directory,id_bin, "gene_models.ffn"),"w")
            if opts.protein:
                file_objects[id_bin]["protein"] = open(os.path.join(opts.output_directory,id_bin, "gene_models.faa"),"w")
            if opts.fasta:
                file_objects[id_bin]["genome"] = open(os.path.join(opts.output_directory,id_bin, "genome.fa"),"w")
    else:
        for id_bin in bins:
            file_objects[id_bin]["gff"] = open(os.path.join(opts.output_directory, "{}.gff".format(id_bin)),"w")
            if opts.cds:
                file_objects[id_bin]["cds"] = open(os.path.join(opts.output_directory, "{}.ffn".format(id_bin)),"w")
            if opts.protein:
                file_objects[id_bin]["protein"] = open(os.path.join(opts.output_directory, "{}.faa".format(id_bin)),"w")
            if opts.fasta:
                file_objects[id_bin]["genome"] = open(os.path.join(opts.output_directory, "{}.fa".format(id_bin)),"w")

    # ===
    # GFF
    # ===
    print("Parsing GFF file: {}".format(opts.gff), file=sys.stderr)
    # Parse
    f_comments = open(os.path.join(opts.output_directory,"gff_comments.txt"), "w")
    identifier_mapping = defaultdict(dict)
    with open(opts.gff, "r") as f_gff:
        for line in f_gff:
            line = line.strip()
            if  line.startswith("#"):
                print(line, file=f_comments)
            else:
                id_scaffold = line.split("\t")[0]
                # id_gene = line.split("gene_id=")[1].replace(";","")
                id_gene = "_".join([ 
                    id_scaffold,
                    line.split("ID=")[1].split(";")[0].split("_")[-1],
                ])
                identifier_mapping[id_gene]["id_contig"] = id_scaffold
                if id_scaffold in scaffold_to_bin:
                    id_bin = scaffold_to_bin[id_scaffold]
                    identifier_mapping[id_gene]["id_mag"] = id_bin
                    print(line, file=file_objects[id_bin]["gff"])
    df_identifiers = pd.DataFrame(identifier_mapping).T
    df_identifiers.index.name = "id_orf"

    if not opts.include_unbinned:
        df_identifiers = df_identifiers.dropna(how="any", axis=0)

    df_identifiers.to_csv(os.path.join(opts.output_directory, "identifier_mapping.tsv"), sep="\t", header=bool(opts.include_header))

    # Close
    f_comments.close()
    for id_bin, files in file_objects.items():
        files["gff"].close()
        
    # ===
    # CDS
    # ===
    if opts.cds:
        print("Parsing CDS file: {}".format(opts.cds), file=sys.stderr)
        # Parse
        with open(opts.cds, "r") as f_cds:
            for header, seq in SimpleFastaParser(f_cds):
                id_gene = header.split(" ")[0]
                id_scaffold = "_".join(id_gene.split("_")[:-1])
                if id_scaffold in scaffold_to_bin:
                    id_bin = scaffold_to_bin[id_scaffold]
                    print(">{}\n{}".format(header, seq), file=file_objects[id_bin]["cds"])
        # Close
        for id_bin, files in file_objects.items():
            files["cds"].close()

    # ===
    # Protein
    # ===
    if opts.protein:
        print("Parsing protein file: {}".format(opts.protein), file=sys.stderr)
        # Parse
        with open(opts.protein, "r") as f_protein:
            for header, seq in SimpleFastaParser(f_protein):
                id_gene = header.split(" ")[0]
                id_scaffold = "_".join(id_gene.split("_")[:-1])
                if id_scaffold in scaffold_to_bin:
                    id_bin = scaffold_to_bin[id_scaffold]
                    print(">{}\n{}".format(header, seq), file=file_objects[id_bin]["protein"])
        # Close
        for id_bin, files in file_objects.items():
            files["protein"].close()

    # ===
    # Fasta
    # ===
    if opts.fasta:
        print("Parsing assembly file: {}".format(opts.fasta), file=sys.stderr)
        # Parse
        with open(opts.fasta, "r") as f_assembly:
            if opts.use_mag_as_description:
                for header, seq in SimpleFastaParser(f_assembly):
                    id_scaffold = header.split(" ")[0]
                    if id_scaffold in scaffold_to_bin:
                        id_bin = scaffold_to_bin[id_scaffold]
                        print(">{} {}\n{}".format(id_scaffold, id_bin, seq), file=file_objects[id_bin]["genome"])
            else:
                for header, seq in SimpleFastaParser(f_assembly):
                    id_scaffold = header.split(" ")[0]
                    if id_scaffold in scaffold_to_bin:
                        id_bin = scaffold_to_bin[id_scaffold]
                        print(">{}\n{}".format(header, seq), file=file_objects[id_bin]["genome"])
        # Close
        for id_bin, files in file_objects.items():
            files["genome"].close()


if __name__ == "__main__":
    main()
    
                
