#!/usr/bin/env python
import sys, os, argparse, gzip 
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.8.29"

# Reverse Complement
def reverse_complement(seq):
    """
    Inputs a sequence of DNA and returns the reverse complement se$
    Outputs a reversed string of the input where A>T, T>A, G>C, an$
    """

    conversion = str.maketrans('ACGTacgt','TGCAtgca')

    comp = seq.translate(conversion)
    rev_comp = comp[::-1]
    return rev_comp


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__

    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <input.fasta> -o <output.fasta>)".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-L","--input", default="stdin", type=str, help = "Input fasta file [Default: stdin]")
    parser.add_argument("-1","--forward_reads",  type=str, help = "Output forwad fastq file")
    parser.add_argument("-2","--reverse_reads",  type=str, help = "Output reverse fastq file")
    parser.add_argument("-I","--interleaved_reads",  default="stdout", type=str, help = "Output interleaved file [Default: stdout]")
    parser.add_argument("-r","--read_length",  type=int, default=300, help = "Length of output reads [Default: 300]")
    
    # parser.add_argument("-m","--minimum_read_length",  type=int, default=150, help = "Minimum length of output reads (end of long read) [Default: 150]")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    paired = False
    assert bool(opts.forward_reads) == bool(opts.reverse_reads), "If --forward_reads are provided then --reverse_reads must also be provided"
    if any([opts.forward_reads, opts.reverse_reads]):
        assert opts.interleaved_reads == "stdout", "Cannot provide --interleaved_reads if --forward_reads and --reverse_reads are provided"
        paired = True
    
    if opts.input == "stdin":
        f_in = sys.stdin
    else:
        if opts.input.endswith(".gz"):
            f_in = gzip.open(opts.input, "rt")
        else:
            f_in = open(opts.input, "r")
  
            
    # Paired end reads
    # @SRR17458614.3972019
    # AATGCGATACCACTAATACTAATACTTCGATAGGAAAAGAGCGACATCCAGAAAACGTTTGGAACAGCAACACAAGGAGAGAGAATGACTTTCCAACACGCAAGGAAGCACGACTAATTTCAAAATCGAAAATGAAATTTTCATTTTCAAT
    # +
    # AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<AEEEEAAEE
    if paired:
        if opts.forward_reads.endswith(".gz"):
            f_forward = gzip.open(opts.forward_reads, "wt")
        else:
            f_forward = open(opts.forward_reads, "w")
        if opts.reverse_reads.endswith(".gz"):
            f_reverse = gzip.open(opts.reverse_reads, "wt")
        else:
            f_reverse = open(opts.reverse_reads, "w") 
            
        for header, seq_long, qual_long in tqdm(FastqGeneralIterator(f_in), f"Coverting fastq input: {opts.input}"):
            id_long, *description = header.split(" ")
            description = " ".join(description)
            
            for i in range(0, len(seq_long), opts.read_length):
                start = i
                end = i + opts.read_length
                id_short = f"{id_long}_{start}-{end}"
                seq_short = seq_long[start:end]
                # if len(seq_short) >= opts.minimum_read_length:
                qual_short = qual_long[start:end]
                
                print(f"@{id_short} {description}\n{seq_short}\n+\n{qual_short}", file=f_forward)
                print(f"@{id_short} {description}\n{reverse_complement(seq_short)}\n+\n{qual_short[::-1]}", file=f_reverse)

        f_forward.close()
        f_reverse.close()
        
    # Interleaved
    else:        
        if opts.interleaved == "stdout":
            f_interleaved = sys.stdout
        else:
            if opts.interleaved.endswith(".gz"):
                f_interleaved = gzip.open(opts.interleaved, "wt")
            else:
                f_reverse = open(opts.interleaved, "w") 
            
        for header, seq_long, qual_long in tqdm(FastqGeneralIterator(f_in), f"Coverting fastq input: {opts.input}"):
            id_long, *description = header.split(" ")
            description = " ".join(description)
            
            for i in range(0, len(seq_long), opts.read_length):
                start = i
                end = i + opts.read_length
                id_short = f"{id_long}_{start}-{end}"
                seq_short = seq_long[start:end]
                # if len(seq_short) >= opts.minimum_read_length:
                qual_short = qual_long[start:end]
                
                print(f"@{id_short}/1 {description}\n{seq_short}\n+\n{qual_short}", file=f_interleaved)
                print(f"@{id_short}/2 {description}\n{reverse_complement(seq_short)}\n+\n{qual_short[::-1]}", file=f_interleaved)
                
        if f_interleaved != sys.stdout:
            f_interleaved.close()
        
            
    if f_in != sys.stdin:
        f_in.close()

if __name__ == "__main__":
    main()
    
                

