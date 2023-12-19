#!/usr/bin/env python
import sys, os, argparse, gzip 
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.12.18"

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
    parser.add_argument("-i","--input", default="stdin", type=str, help = "Input fasta file")
    parser.add_argument("-o","--output", default="stdout", type=str, help = "Output fasta file")
    parser.add_argument("-n", "--name", type=str, required=True, help = "Name to use for pseudo-scaffold")
    parser.add_argument("-N", "--pad", type=int, default=100, help = "Number of N to use for joining contigs")
    parser.add_argument("-d", "--description", type=str, help = "Description to use [Default: Input filepath]")
    parser.add_argument("-m","--minimum_sequence_length", default=1, type=int, help = "Minimum sequence length accepted [Default: 1]")
    parser.add_argument("-w","--wrap", default=1000, type=int, help = "Wrap fasta. Use 0 for no wrapping [Default: 1000]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert opts.minimum_sequence_length > 0
    assert opts.pad >= 0

    # Input
    f_in = None
    if opts.input == "stdin":
        f_in = sys.stdin 
    else:
        if opts.input.endswith(".gz"):
            f_in = gzip.open(opts.input, "rt")
        else:
            f_in = open(opts.input, "r")
    assert f_in is not None

    # Output
    f_out = None
    if opts.output == "stdout":
        f_out = sys.stdout 
    else:
        if opts.output.endswith(".gz"):
            f_out = gzip.open(opts.output, "wt")
        else:
            f_out = open(opts.output, "w")
    assert f_out is not None
    
    # Concatenated assembly

    if not opts.description:
        opts.description = "assembly_filepath: {}".format(opts.input)
    else:
        if opts.description == "NONE":
            opts.description = ""
    pseudoscaffold_header = "{} {}".format(opts.name, opts.description)

    print(">{}".format(pseudoscaffold_header), file=f_out)
    sequences = list()
    for header, seq in tqdm(SimpleFastaParser(f_in), "Reading fasta input"):
        if len(seq) >= opts.minimum_sequence_length:
            sequences.append(seq)
    number_of_sequences = len(sequences)
    sequences = ("N"*opts.pad).join(sequences)
    
    # Open output file
    if opts.wrap > 0:
        for i in range(0, len(sequences), opts.wrap):
            wrapped_sequence = sequences[i:i+opts.wrap]
            # Write header and wrapped sequence
            print(wrapped_sequence, file=f_out)
    else:
        print(sequences, file=f_out)


    # Close
    if f_in != sys.stdin:
        f_in.close()
    if f_out != sys.stdout:
        f_out.close()

if __name__ == "__main__":
    main()
    
                

