#!/usr/bin/env python
import sys, os, argparse, gzip 
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.12.5"

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
    parser.add_argument("-c","--minimum_coverage", default=0, type=int, help = "Minimum coverage accepted [Default: 0.0]")
    parser.add_argument("-m","--minimum_sequence_length", default=1, type=int, help = "Minimum sequence length accepted [Default: 1]")

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
    
    if opts.retain_description:
        for header, seq in tqdm(SimpleFastaParser(f_in), "Reading fasta input"):
            id = header.split(" ")[0].strip()
            fields = id.split("_")
            try:
                length_index = fields.index("length")
                coverage_index = fields.index("cov")
            except ValueError:
                raise "Your fastq identifiers do not look like they are from SPAdes: {}".format(id)
                sys.exit(1)
            assert ">" not in seq, "`{}` has a '>' character in the sequence which will cause an error.  This can arise from concatenating fasta files where a record is missing a final linebreak".format(header)
            coverage = float(fields[coverage_index + 1])
            length = int(fields[length_index + 1])
            if all([coverage >= opts.minimum_coverage, length >= opts.minimum_sequence_length]):
                print(">{}\n{}".format(header,seq), file=f_out)
    else:
        for header, seq in tqdm(SimpleFastaParser(f_in), "Reading fasta input"):
            id = header.split(" ")[0].strip()
            fields = id.split("_")
            try:
                length_index = fields.index("length")
                coverage_index = fields.index("cov")
            except ValueError:
                raise "Your fastq identifiers do not look like they are from SPAdes: {}".format(id)
                sys.exit(1)
            assert ">" not in seq, "`{}` has a '>' character in the sequence which will cause an error.  This can arise from concatenating fasta files where a record is missing a final linebreak".format(header)
            coverage = float(fields[coverage_index + 1])
            length = int(fields[length_index + 1])
            if all([coverage >= opts.minimum_coverage, length >= opts.minimum_sequence_length]):
                print(">{}\n{}".format(id,seq), file=f_out)

    # Close
    if f_in != sys.stdin:
        f_in.close()
    if f_out != sys.stdout:
        f_out.close()

if __name__ == "__main__":
    main()
    
                

