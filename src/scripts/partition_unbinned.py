#!/usr/bin/env python
import sys, os, glob, argparse 
from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.12.18"

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
    parser.add_argument("-b","--bins", type=str, required=True, help = "path/to/bins.list, No header")
    parser.add_argument("-f","--fasta", type=str, required=True, help = "path/to/fasta")
    parser.add_argument("-o","--output", type=str, default="stdout",  help = "Output fasta file [Default: stdout]")
    parser.add_argument("-m", "--minimum_contig_length", type=int, default=1, help="Minimum contig length.  [Default: 1] ")
    parser.add_argument("--mode", type=str, default="unbinned", help="Get 'unbinned' or 'binned' contigs [Default: 'unbinned'] ")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    assert opts.mode in {"unbinned", "binned"}, "--mode must be 'unbinned' or 'binned'"

    # Output
    if opts.output == "stdout":
        f_out = sys.stdout 
    else:
        f_out = open(opts.output, "w")

    # Bins
    print("Reading bins file: {}".format(opts.bins), file=sys.stderr)
    bins = set()
    with open(opts.bins, "r") as f_bins:
        for line in f_bins.readlines():
            line = line.strip()
            if line:
                id_bin = line
                bins.add(id_bin)

    # Scaffolds
    print("Reading scaffolds to bins file: {}".format(opts.scaffolds_to_bins), file=sys.stderr)
    scaffolds = set()
    with open(opts.scaffolds_to_bins, "r") as f_s2b:
        for line in f_s2b.readlines():
            line = line.strip()
            if line:
                id_scaffold, id_bin = line.split("\t")
                if id_bin in bins:
                    scaffolds.add(id_scaffold)

    # Fasta
    print("Reading fasta file: {}".format(opts.fasta), file=sys.stderr)
    func_filter = {
        "unbinned": lambda id: id not in scaffolds,
        "binned": lambda id: id in scaffolds,
    }[opts.mode]
    with open(opts.fasta, "r") as f_fasta:
        for header, seq in SimpleFastaParser(f_fasta):
            if len(seq) >= opts.minimum_contig_length:
                id = header.split(" ")[0]
                if func_filter(id):
                    print(">{}\n{}".format(header, seq), file=f_out)


    
    if f_out is not sys.stdout: 
        f_out.close()

if __name__ == "__main__":
    main()
    
                
