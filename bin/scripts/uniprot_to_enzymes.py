#!/usr/bin/env python
import sys, os, argparse, gzip, warnings, re, pickle
# from collections import OrderedDict
# import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.4.23"

# Define the regex pattern
pattern = r'\(EC [\d.-]+\)'
# Compile the pattern
regex = re.compile(pattern)
def parse_enzyme_commissions(description, regex=regex):
    # Find all matches in the input string
    matches = regex.findall(description)
    if matches:
        # Extract only the alphanumeric sequences from the matches
        return {re.sub(r'[^\w.-]', '', match)[2:] for match in matches}
    

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <input.tsv[.gz]> -o <output.tsv[.gz] -d <uniprot_to_ecs.dict.pkl".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-i","--input", type=str, help = "path/to/input.tsv[.gz] UniProt output from https://www.uniprot.org/uniprotkb?query=ec%3A* [Default: stdin]")
    parser_io.add_argument("-o","--output", type=str, help = "path/to/output.tsv[.gz] formatted output containing 2 columns [id_uniprot]<tab>[EC:[id_ec],[id_ec]] [Default: stdout]")
    parser_io.add_argument("--header", action="store_true", help="Specify if header should be in output")

    parser_serialization = parser.add_argument_group('Serialization arguments')
    parser_serialization.add_argument("-d", "--export_dict", type=str,   help = "prefix/to/dict pickled output file: {id_uniprot: {id_ec1, id_ec2, ...}} suggested prefix is .dict.pkl")

    parser_fasta = parser.add_argument_group('Fasta arguments')
    parser_fasta.add_argument("--fasta_input", type=str, help = "prefix/to/proteins.fasta[.gz]")
    parser_fasta.add_argument("--fasta_output", type=str, help = "prefix/to/proteins.reformatted.fasta[.gz]")
    # parser_fasta.add_argument("--identifier_index", type=int, default=1, help = "Identifier index separated by | characters. (e.g., sp|A0A009IHW8|ABTIR_ACIB9 selecting 1 would yield A0A009IHW8) [Default: 1]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert bool(opts.fasta_input) == bool(opts.fasta_output), "If --fasta_input is provided then --fasta_output must be provided and vice versa"

    # Input
    if opts.input == "stdin":
        f_input = sys.stdin
    else:
        if opts.input.endswith(".gz"):
            f_input = gzip.open(opts.input, "rt")
        else:
            f_input = open(opts.input, "r")
    
    # Output
    if opts.output == "stdout":
        f_output = sys.stdout
    else:
        if opts.output.endswith(".gz"):
            f_output = gzip.open(opts.output, "wb")
        else:
            f_output = open(opts.output, "w")
    if opts.header:
        print("id_uniprot", "enzymes", sep="\t", file=f_output)

    # Parse
    uniprot_to_reformattedlabel = dict()

    # Output table and export dictionary pickle
    if opts.export_dict:
        if opts.export_dict.endswith((".gz",".pgz")):
            f_pickle = gzip.open(opts.export_dict, "wb")
        else:
            f_pickle = open(opts.export_dict, "wb")

        uniprot_to_ecs = dict()
        for line in tqdm(f_input, desc="Reading from file: {}".format(f_input)):
            line = line.strip()
            if line:
                fields = line.split("\t")
                id_uniprot = fields[0]
                description = fields[3]
                enzymes = parse_enzyme_commissions(description)
                if enzymes:
                    uniprot_to_ecs[id_uniprot] = enzymes
                    label = "EC:{}".format(",".join(sorted(enzymes)))
                    uniprot_to_reformattedlabel[id_uniprot] = label
                    print(id_uniprot, label, sep="\t", file=f_output)
        pickle.dump(uniprot_to_ecs, f_pickle)
        f_pickle.close()

    # Output table only
    else:
        for line in tqdm(f_input, desc="Reading from file: {}".format(f_input)):
            line = line.strip()
            if line:
                fields = line.split("\t")
                id_uniprot = fields[0]
                description = fields[3]
                enzymes = parse_enzyme_commissions(description)
                if enzymes:
                    label = "EC:{}".format(",".join(sorted(enzymes)))
                    uniprot_to_reformattedlabel[id_uniprot] = label
                    print(id_uniprot, label, sep="\t", file=f_output)

    
    # Fasta
    if opts.fasta_input:
        from Bio.SeqIO.FastaIO import SimpleFastaParser
        # Fasta Input
        if opts.fasta_input.endswith(".gz"):
            f_fasta_input = gzip.open(opts.fasta_input, "rt")
        else:
            f_fasta_input = open(opts.fasta_input, "r")
        # Fasta Output
        if opts.fasta_output.endswith(".gz"):
            f_fasta_output = gzip.open(opts.fasta_output, "wb")
        else:
            f_fasta_output = open(opts.fasta_output, "w") 

        # Reformat labels
        for header, seq in tqdm(SimpleFastaParser(f_fasta_input), desc="Reading from file: {}".format(f_fasta_input)):
            id_long = header.split(" ")[0]
            id = id_long.split("|")[1]
            if id in uniprot_to_reformattedlabel:
                print(">{}_{}\n{}".format(id, uniprot_to_reformattedlabel[id], seq), file=f_fasta_output)

        f_fasta_input.close()
        f_fasta_output.close()
    # Closing
    if f_input != sys.stdin:
        f_input.close()
    if f_output != sys.stdout:
        f_output.close()

if __name__ == "__main__":
    main()
    
                

