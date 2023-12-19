#!/usr/bin/env python
import sys, os, argparse, gzip 
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.11.10"

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
    parser.add_argument("-r","--retain_description", action="store_true", help = "Retain description")
    parser.add_argument("-s","--retain_stop_codon", action="store_true", help = "Retain stop codon character (if one exists)")
    parser.add_argument("-m","--minimum_sequence_length", default=1, type=int, help = "Minimum sequence length accepted [Default: 1]")
    parser.add_argument("--stop_codon_character", default="*", type=str, help = "Stop codon character [Default: *] ")
    # parser.add_argument("-t","--molecule_type",  help = "Comma-separated list of names for the --scaffolds_to_bins")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert opts.minimum_sequence_length > 0

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
    
    # retain_description=True
    # retain_stop_codon=True
    if all([
        opts.retain_description,
        opts.retain_stop_codon,
        ]):
        for header, seq in tqdm(SimpleFastaParser(f_in), "Reading fasta input"):
            header = header.strip()
            if len(seq) >= opts.minimum_sequence_length:
                assert ">" not in seq, "`{}` has a '>' character in the sequence which will cause an error.  This can arise from concatenating fasta files where a record is missing a final linebreak".format(header)
                print(">{}\n{}".format(header,seq), file=f_out)

    # retain_description=False
    # retain_stop_codon=True
    if all([
        not opts.retain_description,
        opts.retain_stop_codon,
        ]):
        for header, seq in tqdm(SimpleFastaParser(f_in), "Reading fasta input"):
            id = header.split(" ")[0].strip()
            if len(seq) >= opts.minimum_sequence_length:
                assert ">" not in seq, "`{}` has a '>' character in the sequence which will cause an error.  This can arise from concatenating fasta files where a record is missing a final linebreak".format(header)
                print(">{}\n{}".format(id,seq), file=f_out)

    # retain_description=True
    # retain_stop_codon=False
    if all([
        opts.retain_description,
        not opts.retain_stop_codon,
        ]):
        for header, seq in tqdm(SimpleFastaParser(f_in), "Reading fasta input"):
            header = header.strip()
            if seq.endswith(opts.stop_codon_character):
                seq = seq[:-1]
            if len(seq) >= opts.minimum_sequence_length:
                assert ">" not in seq, "`{}` has a '>' character in the sequence which will cause an error.  This can arise from concatenating fasta files where a record is missing a final linebreak".format(header)
                print(">{}\n{}".format(header,seq), file=f_out)

    # retain_description=False
    # retain_stop_codon=False
    if all([
        not opts.retain_description,
        not opts.retain_stop_codon,
        ]):
        for header, seq in tqdm(SimpleFastaParser(f_in), "Reading fasta input"):
            id = header.split(" ")[0].strip()
            if seq.endswith(opts.stop_codon_character):
                seq = seq[:-1]
            if len(seq) >= opts.minimum_sequence_length:
                assert ">" not in seq, "`{}` has a '>' character in the sequence which will cause an error.  This can arise from concatenating fasta files where a record is missing a final linebreak".format(header)
                print(">{}\n{}".format(id,seq), file=f_out)

    # Close
    if f_in != sys.stdin:
        f_in.close()
    if f_out != sys.stdout:
        f_out.close()

if __name__ == "__main__":
    main()
    
                

