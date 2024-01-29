#!/usr/bin/env python
import sys, os, glob, argparse, gzip, warnings
from collections import OrderedDict
import pandas as pd
from tqdm import tqdm
from Bio.SeqIO.FastaIO import SimpleFastaParser

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.7.7"

def gc_content(seq):
    seq = seq.upper()
    number_of_gc = seq.count("G") + seq.count("C")
    return number_of_gc/len(seq)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__

    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -n <genome_id> -f <genome_assembly.fasta> -o <output.gff[.gz]> -c <gff_cds(Pyrodigal|MetaEuk)> -r <gff_rRNA(BARRNAP)> -t <gff_tRNA(tRNAscan-SE) -d <organism_type>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-f", "--fasta", type=str, default="stdin", help = "path/to/assembly.fasta.  Adds contig regions to output GFF. [Default: stdin]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.gff[.gz] [Default: stdout]")
    parser.add_argument("-n", "--name", type=str, required=True, help="Genome ID")
    parser.add_argument("-c","--gff_cds", type=str, required=True, help = "GFF for CDS genes (Pyrodigal|MetaEuk - modified)")
    parser.add_argument("-r","--gff_rRNA", type=str, help = "GFF for rRNA genes (BARRNAP - modified) [Optional]")
    parser.add_argument("-t","--gff_tRNA", type=str, help = "GFF for tRNA genes (tRNAscan-SE - modified) [Optional]")
    parser.add_argument("-d","--organism_type", type=str,  help = "Organism type e.g., {Prokaryote, Eukaryote, Bacteria, Archaea, Virus} [Optional]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Output 
    if opts.output == "stdout":
        f_out = sys.stdout
    else:
        if opts.output.endswith(".gz"):
            f_out = gzip.open(opts.output, "wt")
        else:
            f_out = open(opts.output, "w")
    
    # GFF (Output)
    print("##gff-version 3", file=f_out)
    print("##Program: VEBA (github.com/jolespin/veba)", file=f_out)
    print("##ID: {}".format(opts.name), file=f_out)
    print("##Source: Metagenome-Assembled Genome", file=f_out)
    if opts.organism_type:
        print("##Organism-type: {}".format(opts.organism_type), file=f_out)

    # Contigs
    if opts.fasta == "stdin":
        f_fa = sys.stdin
    else:
        f_fa = open(opts.fasta, "r")

    for header, seq in tqdm(SimpleFastaParser(f_fa), "Reading fasta file: {}".format(opts.fasta)):
        id = header.split(" ")[0]
        attributes = "ID={};genome_id={};gc_cont={:.3f};".format(id, opts.name, gc_content(seq))
        gff_fields = [id, "VEBA", "region", 1, len(seq), ".", "+", ".", attributes]
        print(*gff_fields, sep="\t", file=f_out)

    if f_fa != sys.stdin:
        f_fa.close()

    # GFF (CDS)
    with open(opts.gff_cds, "r") as f_gff:
        for line in tqdm(f_gff, "Reading GFF file (CDS): {}".format(opts.gff_cds)):
            line = line.strip()
            if not line.startswith("#"):
                print(line, file=f_out)

    # GFF (rRNA)
    if opts.gff_rRNA:
        with open(opts.gff_rRNA, "r") as f_gff:
            for line in tqdm(f_gff, "Reading GFF file (rRNA): {}".format(opts.gff_rRNA)):
                line = line.strip()
                if not line.startswith("#"):
                    print(line, file=f_out)
    else:
        warnings.warn("No --gff_rRNA was provided")

    # GFF (tRNA)
    if opts.gff_tRNA:
        with open(opts.gff_tRNA, "r") as f_gff:
            for line in tqdm(f_gff, "Reading GFF file (tRNA): {}".format(opts.gff_tRNA)):
                line = line.strip()
                if not line.startswith("#"):
                    if not line.endswith(";"):
                        line += ";"
                    print(line, file=f_out)
    else:
        warnings.warn("No --gff_tRNA was provided")

    if f_out != sys.stdout:
        f_out.close()

if __name__ == "__main__":
    main()
    
                