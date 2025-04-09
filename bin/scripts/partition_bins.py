#!/usr/bin/env python
import sys, os, glob, argparse 
from collections import defaultdict
# from Bio.SeqIO.FastaIO import SimpleFastaParser
from pyexeggutor import (
    open_file_writer,
    fasta_writer,
)
import pyfastx

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2025.4.8"

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
    parser.add_argument("-f","--fasta", type=str, required=True, help = "path/to/fasta")
    parser.add_argument("-o","--output_directory", type=str, required=True,  help = "Output directory for bins")
    parser.add_argument("-x","--extension", type=str, default="fa",  help = "File extension for bin fasta.  Supports gzip [Default: fa]")
    parser.add_argument("-w", "--wrap", type=int, default=1000, help = f"Line width for wrapping fasta lines.  If 0, then no wrapping is used. [Default: 1000]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Output
    os.makedirs(opts.output_directory, exist_ok=True)

    # Scaffolds
    print("Reading scaffolds to bins file: {}".format(opts.scaffolds_to_bins), file=sys.stderr)
    scaffold_to_bin = dict()
    with open(opts.scaffolds_to_bins, "r") as f:
        for line in f.readlines():
            line = line.strip()
            if line:
                id_scaffold, id_bin = line.split("\t")
                scaffold_to_bin[id_scaffold] = id_bin
    # Files
    bin_to_file = dict()
    for id_bin in set(scaffold_to_bin.values()):
        filepath = os.path.join(opts.output_directory, "{}.{}".format(id_bin, opts.extension))
        bin_to_file[id_bin] = open_file_writer(filepath)
        
    # Fasta
    print("Reading fasta file: {}".format(opts.fasta), file=sys.stderr)
    for header, seq in pyfastx.Fasta(opts.fasta, build_index=False, full_name=True):
        id = header.split(" ")[0]
        id_bin = scaffold_to_bin[id]
        f_out = bin_to_file[id_bin]
        fasta_writer(
            header=header,
            seq=seq,
            file=f_out,
            wrap=opts.wrap,
        )

    # Close
    for id_bin in bin_to_file:
        bin_to_file[id_bin].close()

if __name__ == "__main__":
    main()
    
                
