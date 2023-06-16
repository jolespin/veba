#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, gzip
from collections import OrderedDict
import pandas as pd
from tqdm import tqdm
from Bio.SeqIO.FastaIO import SimpleFastaParser

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.6.14"

#
def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <index.list> -t <table.tsv> -o <output.tsv>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", default="stdin", type=str, help = "path/to/table.tsv [Default: stdin]")
    parser.add_argument("-o","--output", default="stdout", type=str, help = "path/to/output.fasta [Default: stdout]")
    parser.add_argument("-n","--name_index", default=0, type=int, help = "Name index position [Default: 0]")
    parser.add_argument("-d","--description_index", type=int, help = "Description index position")
    parser.add_argument("-s","--sequence_index", default=2, type=int, help = "Sequence index position.  Not used if --fasta is provided. [Default: 2]")
    parser.add_argument("-f","--fasta",  type=str, help = "Separate fasta file to use instead of column")
    parser.add_argument("-x","--fasta_id_index",  type=int,  help = "Index to use for fasta ID in table.  Most likely will be --name_index or --description_index.  Required if --fasta is used.")
    parser.add_argument("--no_header",action="store_true", help = "Use if there is no header")
    parser.add_argument("--sep", type=str, default="\t", help = "Separator [Default: <tab>]")
    parser.add_argument("--skiprows", type=int, help = "Skiprows")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # I/O
    if opts.input == "stdin":
        opts.input = sys.stdin 

    if opts.output == "stdout":
        opts.output = sys.stdout
    else:
        if opts.output.endswith(".gz"):
            opts.output = gzip.open(opts.output, "wt")
        else:
            opts.output = open(opts.output, "w")

    # Read Table
    if opts.no_header:
        df_input = pd.read_csv(opts.input, sep=opts.sep, index_col=None, skiprows=opts.skiprows, header=None)
    else:
        df_input = pd.read_csv(opts.input, sep=opts.sep, index_col=None, skiprows=opts.skiprows)#, header=bool(opts.no_header))

    if opts.description_index is not None:
        if opts.fasta:
            assert opts.fasta_id_index is not None, "If --fasta is provided then --fasta_id_index must also be provided"
            assert max(opts.name_index, opts.description_index, opts.fasta_id_index) < df_input.shape[1]

            index = df_input.iloc[:,opts.fasta_id_index].values

            
            assert len(set(index)) == len(index), "--fasta_id_index contains duplicates which is not allowed for fasta files"

            df = df_input.iloc[:,[opts.name_index, opts.description_index]]
            file_open_function = {True:gzip.open, False:open}[opts.fasta.endswith(".gz")]
            
            id_to_seq = dict()
            with file_open_function(opts.fasta, "rt") as f_fasta:
                for header, seq in tqdm(SimpleFastaParser(f_fasta), "Reading fasta sequences"):
                    id = header.split(" ")[0]
                    if id in index:
                        # assert id in index, f"{id} not in --input table under index {opts.fasta_id_index}"
                        id_to_seq[id] = seq
            id_to_seq = pd.Series(id_to_seq).loc[index]
            assert id_to_seq.isnull().sum() == 0, f"There are ids in --input table under column index {opts.fasta_id_index} that do not have sequences in --fasta"
            df[2] = id_to_seq.values
        else:
            assert max(opts.name_index, opts.description_index, opts.sequence_index) < df_input.shape[1]
            df = df_input.iloc[:,[opts.name_index, opts.description_index, opts.sequence_index]]

        for i, (id, description, seq) in tqdm(df.iterrows(), "Writing sequence records", total=df.shape[0]):
            print(">{} {}\n{}".format(id, description, seq), file=opts.output)
    else:
        if opts.fasta:
            assert opts.fasta_id_index is not None
            assert max(opts.name_index,  opts.sequence_index, opts.fasta_id_index) < df_input.shape[1]

            index = df_input.iloc[:,opts.fasta_id_index].values

            assert len(set(index)) == len(index), "--fasta_id_index contains duplicates which is not allowed for fasta files"

            df = df_input.iloc[:,[opts.name_index]]
            file_open_function = {True:gzip.open, False:open}[opts.fasta.endswith(".gz")]

            id_to_seq = dict()
            with file_open_function(opts.fasta, "rt") as f_fasta:
                for header, seq in tqdm(SimpleFastaParser(f_fasta), "Reading fasta sequences"):
                    id = header.split(" ")[0]
                    if id in index:
                        # assert id in index, f"{id} not in --input table under index {opts.fasta_id_index}"
                        id_to_seq[id] = seq
            id_to_seq = pd.Series(id_to_seq).loc[index]
            assert id_to_seq.isnull().sum() == 0, f"There are ids in --input table under column index {opts.fasta_id_index} that do not have sequences in --fasta"
            df[1] = id_to_seq.values

        else:
            assert max(opts.name_index,  opts.sequence_index) < df_input.shape[1]
            df = df_input.iloc[:,[opts.name_index, opts.sequence_index]]
            
        for i, (id, seq) in tqdm(df.iterrows(), "Writing sequence records", total=df.shape[0]):
            print(">{}\n{}".format(id, seq), file=opts.output)

    if opts.output != sys.stdout:
        opts.output.close()
if __name__ == "__main__":
    main()
