#!/usr/bin/env python
import sys, os, argparse, gzip, hashlib, warnings
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.3.12"

# Thanks to @shenwei356 for hash function:
# https://github.com/shenwei356/seqkit/issues/444#issuecomment-1987356376

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__

    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <input.fasta> -o <output.fasta> -s <source_organism>)".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jol.espinoz@gmail.com)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-i","--input", default="stdin", type=str, help = "Input fasta file")
    parser.add_argument("-o","--output", default="stdout", type=str, help = "Output fasta file")
    parser.add_argument("-s","--source_organism", required=True, type=str, help = "Source organism")
    parser.add_argument("-m","--mode", required=True, choices={"numeric", "hash"}, type=str, help = "`hash` is useful for dereplicating sequences for building databases but will remove duplicate proteins in a genome")
    parser.add_argument("-t","--tsv", type=str, help = "Tabular output for identifier mapping [id_original]<tab>[id_new]")
    parser.add_argument("-d","--delimiter", default="|", type=str, help = "Delimiter between source organism and md5hash. Duplicate sequences are removed. [id_source_organism][delimiter][md5hash] [id_original] [Default: | ]")
    parser.add_argument("-p", "--prefix", default="g", type=str, help = "Prefix for protein if --mode numeric [Default: g]")
    parser.add_argument("--minimum_sequence_length", default=1, type=int, help = "Minimum sequence length accepted [Default: 1]")
    parser.add_argument("--stop_codon_character", default="*", type=str, help = "Stop codon character [Default: *] ")


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

    # Table output
    f_tsv = None
    if opts.tsv:
        if opts.tsv.endswith(".gz"):
            f_tsv = gzip.open(opts.tsv, "wt")
        else:
            f_tsv = open(opts.tsv, "w")
    
    if opts.mode == "hash":
        hashes = set()
        # tsv_output=True
        if f_tsv:
            for header, seq in tqdm(SimpleFastaParser(f_in), "Reading fasta input: {}".format(opts.input)):
                header = header.strip()
                id_original = header.split(" ")[0]
            
                if seq.endswith(opts.stop_codon_character):
                    seq = seq[:-1]
                
                if len(seq) >= opts.minimum_sequence_length:
                    assert ">" not in seq, "`{}` has a '>' character in the sequence which will cause an error.  This can arise from concatenating fasta files where a record is missing a final linebreak".format(header)
                    seq = seq.upper()
                    id_hash = hashlib.md5(seq.encode("utf-8")).hexdigest()
                    if id_hash not in hashes:
                        id_reformatted = "{}{}{}".format(opts.source_organism, opts.delimiter, id_hash)
                        header_reformatted = "{} {}".format(id_reformatted, id_original)
                        print(">{}\n{}".format(header_reformatted,seq), file=f_out)
                        print(id_reformatted, id_original, sep="\t", file=f_tsv)
                        hashes.add(id_hash)
                    else:
                        print("[Duplicate Removed] {} {}".format(id_original, id_hash), file=sys.stderr)

        # tsv_output=False
        else:
            for header, seq in tqdm(SimpleFastaParser(f_in), "Reading fasta input: {}".format(opts.input)):
                header = header.strip()
                id_original = header.split(" ")[0]
            
                if seq.endswith(opts.stop_codon_character):
                    seq = seq[:-1]
                
                if len(seq) >= opts.minimum_sequence_length:
                    assert ">" not in seq, "`{}` has a '>' character in the sequence which will cause an error.  This can arise from concatenating fasta files where a record is missing a final linebreak".format(header)
                    seq = seq.upper()
                    id_hash = hashlib.md5(seq.encode("utf-8")).hexdigest()
                    if id_hash not in hashes:
                        id_reformatted = "{}{}{}".format(opts.source_organism, opts.delimiter, id_hash)
                        header_reformatted = "{} {}".format(id_reformatted, id_original)
                        print(">{}\n{}".format(header_reformatted,seq), file=f_out)
                        hashes.add(id_hash)
                    else:
                        print("[Duplicate Removed] {} {}".format(id_original, id_hash), file=sys.stderr)

    if opts.mode == "numeric":
        # tsv_output=True
        if f_tsv:
            for i, (header, seq) in tqdm(enumerate(SimpleFastaParser(f_in), start=1), "Reading fasta input: {}".format(opts.input)):
                header = header.strip()
                id_original = header.split(" ")[0]
            
                if seq.endswith(opts.stop_codon_character):
                    seq = seq[:-1]
                
                if len(seq) >= opts.minimum_sequence_length:
                    assert ">" not in seq, "`{}` has a '>' character in the sequence which will cause an error.  This can arise from concatenating fasta files where a record is missing a final linebreak".format(header)
                    seq = seq.upper()

                    id_reformatted = "{}{}{}{}".format(opts.source_organism, opts.delimiter, opts.prefix, i)
                    header_reformatted = "{} {}".format(id_reformatted, id_original)
                    print(">{}\n{}".format(header_reformatted,seq), file=f_out)
                    print(id_reformatted, id_original, sep="\t", file=f_tsv)


        # tsv_output=False
        else:
            for i, (header, seq) in tqdm(enumerate(SimpleFastaParser(f_in), start=1), "Reading fasta input: {}".format(opts.input)):
                header = header.strip()
                id_original = header.split(" ")[0]
            
                if seq.endswith(opts.stop_codon_character):
                    seq = seq[:-1]
                
                if len(seq) >= opts.minimum_sequence_length:
                    assert ">" not in seq, "`{}` has a '>' character in the sequence which will cause an error.  This can arise from concatenating fasta files where a record is missing a final linebreak".format(header)
                    seq = seq.upper()

                    id_reformatted = "{}{}{}{}".format(opts.source_organism, opts.delimiter, opts.prefix, i)
                    header_reformatted = "{} {}".format(id_reformatted, id_original)
                    print(">{}\n{}".format(header_reformatted,seq), file=f_out)

    # Close
    if f_in != sys.stdin:
        f_in.close()
    if f_out != sys.stdout:
        f_out.close()
    if f_tsv is not None:
        f_tsv.close()

if __name__ == "__main__":
    main()
    
                

