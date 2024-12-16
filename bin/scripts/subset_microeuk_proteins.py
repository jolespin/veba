#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, gzip, pickle
from tqdm import tqdm
from pandas import read_pickle
import pyfastx

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.10.2"

    
#
def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <identifiers.list> -o <output.fasta> -d <database.pkl.gz>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    
    # Pipeline
    parser_io = parser.add_argument_group('I/O')
    parser.add_argument("-i","--identifiers", default="stdin", type=str, help = "path/to/identifiers.list|.tsv [Default: stdin]")
    parser.add_argument("-o","--output", default="stdout", type=str, help = "path/to/output.fasta [Default: stdout]")
    parser.add_argument("--input_type", choices={"source", "target", "batch-source", "batch-target"}, default="target", type=str, help = "Input type for identifiers where batch refers to tabular format: [id_genome]<tab>[id_source/id_target]")
    parser.add_argument("--output_type", choices={"fasta", "table"}, default="fasta", type=str, help = "Output format")

    parser_database = parser.add_argument_group('Database')
    parser_database.add_argument("-s","--sequences", type=str, help = "path/to/MicroEuk100.faa.gz which is a lower memory footprint but takes much longer") # Lower memory footprint
    parser_database.add_argument("-d","--sequence_database",  type=str, help = "path/to/MicroEuk100.source_to_target_to_sequence.sources_from_eukaryota_odb10.dict.pkl.gz") 
    parser_database.add_argument("-t","--target_database", type=str, help = "path/to/target_to_source.dict.pkl.gz")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    # Checks
    if opts.input_type in {"batch-source", "batch-target"}:
        raise Exception(f"--input_type {opts.input_type} not yet supported")
    if opts.input_type == "target":
        if not opts.target_database:
            raise Exception("If --input_type 'target' then --target_databas must be provided")
        if not os.path.exists(opts.target_database):
            raise FileNotFoundError(f"If --input_type 'target' then --target_databas must be provided. The following file does not exist: {opts.target_database}")
    assert bool(opts.sequences) != bool(opts.sequence_database), "Must provide either --sequences or --sequence_database but not both"
    assert any([bool(opts.sequences), bool(opts.sequence_database)]), "Must provide either --sequences or --sequence_database but not both"

    # I/O
    if opts.identifiers == "stdin":
        f_in = sys.stdin 
    else:
        if opts.identifiers.endswith(".gz"):
            f_in = gzip.open(opts.identifiers, "rt")
        else:
            f_in = open(opts.identifiers, "r")

    if opts.output == "stdout":
        f_out = sys.stdout 
    else:
        if opts.output.endswith(".gz"):
            f_out = gzip.open(opts.output, "wt")
        else:
            f_out = open(opts.output, "w")       
    
    if opts.input_type == "target":
        print(f"Loading target database: {opts.target_database}", file=sys.stderr)
        target_to_source = read_pickle(opts.target_database)
        
        source_identifiers = set()
        for line in tqdm(f_in, desc=f"Reading identifiers [type={opts.input_type}]: {f_in}"):
            line = line.strip()
            if line:
                id_source = target_to_source[line]
                source_identifiers.add(id_source)
        del target_to_source
    elif opts.input_type == "source":
        source_identifiers = set()
        for line in tqdm(f_in, desc=f"Reading identifiers [type={opts.input_type}]: {f_in}"):
            line = line.strip()
            if line:
                id_source = str(line)
                source_identifiers.add(id_source)
                
    # Extract from sequence database
    if opts.sequence_database:
        print(f"Loading sequence database: {opts.sequence_database}", file=sys.stderr)
        source_to_target_to_sequence = read_pickle(opts.sequence_database)
        if opts.output_type == "fasta":
            for id_source in tqdm(source_identifiers, desc=f"Writing output [{opts.output_type}]: {f_out}"):
                for id_target, seq in source_to_target_to_sequence[id_source].items():
                    print(f">{id_target} {id_source}\n{seq}", file=f_out)
        elif opts.output_type == "table":
            for id_source in tqdm(source_identifiers, desc=f"Writing output [{opts.output_type}]: {f_out}"):
                for id_target, seq in source_to_target_to_sequence[id_source].items():
                    print(id_target, id_source, seq, sep="\t", file=f_out)
                    
    # Extract fasta sequences    
    elif opts.sequences:
        print(f"Reading sequences: {opts.sequences}", file=sys.stderr)
        if opts.output_type == "fasta":
            for id_target, seq in tqdm(pyfastx.Fasta(opts.sequences, build_index=False), desc=f"Writing output [{opts.output_type}]: {f_out}"):
                id_source = target_to_source.pop(id_target)
                if id_source in source_identifiers:
                    print(f">{id_target} {id_source}\n{seq}", file=f_out)

        elif opts.output_type == "table":
            for id_target, seq in tqdm(pyfastx.Fasta(opts.sequences, build_index=False), desc=f"Writing output [{opts.output_type}]: {f_out}"):
                id_source = target_to_source.pop(id_target)
                if id_source in source_identifiers:
                    print(id_target, id_source, seq, sep="\t", file=f_out)
    
    if f_in != sys.stdin:
        f_in.close()
    if f_out != sys.stdout:
        f_out.close()

if __name__ == "__main__":
    main()
