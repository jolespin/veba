#!/usr/bin/env python
import sys
import os
import argparse
import gzip
from collections import defaultdict
# import glob
from tqdm import tqdm
import pysam

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2025.5.22"

def main():
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i  <bam> -o <counts.tsv>".format(__program__)
    epilog = "Copyright 2025 Josh L. Espinoza"


    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", type=str,  help = "Input BAM")
    parser.add_argument("-o","--output", type=str, default="stdout",  help = "Output tsv [Default: stdout]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    # Output
    if opts.output == "stdout":
        f_out = sys.stdout
    else:
        if opts.output.endswith(".gz"):
            f_out = gzip.open(opts.output, "rt")
        else:
            f_out = open(opts.output, "w")

    # Run
    counts = defaultdict(int)
    with pysam.AlignmentFile(opts.input, "rb") as bamfile:
        for alignment in tqdm(bamfile, f"Reading {opts.input}", unit=" alignments"):
            # Filter for primary alignments:
            # A primary alignment is NOT secondary (flag 0x100)
            # AND NOT supplementary (flag 0x800).
            if not alignment.is_secondary and not alignment.is_supplementary:
                # Also, ensure the alignment is actually mapped to a reference
                if not alignment.is_unmapped:
                    # id_read = alignment.query_name      # QNAME field (Read ID)
                    id_contig = alignment.reference_name # RNAME field (Reference ID)
                    counts[id_contig] += 1
    # Write
    for id_contig, count in tqdm(counts.items(), f"Writing counts: {f_out}", unit= " contigs"):
        print(id_contig, count, sep="\t", file=f_out)
        
    if f_out != sys.stdout:
        f_out.close()
        
if __name__ == "__main__":
    main()
    
