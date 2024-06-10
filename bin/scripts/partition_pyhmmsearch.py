#!/usr/bin/env python
import sys, os, glob, argparse, gzip
from collections import defaultdict
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.6.7"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <pyhmmsearch_results.tsv> -a <proteins> -o <output_directory> -s '|--|'".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--pyhmmsearch_results", type=str, help = "path/to/pyhmmsearch_output.tsv[.gz] assumes that `-s/--scores_cutoff` was used with `pyhmmsearch.py`", required=True)
    parser.add_argument("-a","--proteins", type=str, help = "path/to/proteins.faa[.gz]",  required=True)
    parser.add_argument("-o","--output_directory", type=str, help = "Output directory [Default: Directory of --pyhmmsearch_results")
    parser.add_argument("-d", "--sep", default="|--|", type=str, help="[id_organism]<sep>[id_protein] [Default: |--|]")
    parser.add_argument("-g", "--gzip", action="store_true", help="Gzip the protein fasta files")
    parser.add_argument("--no_header", action="store_true", help = "No header on --pyhmmsearch_results")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    if not opts.output_directory:
        opts.output_directory = os.path.split(opts.pyhmmsearch_results)[0]
    
    os.makedirs(opts.output_directory, exist_ok=True)


    # Get marker to query protein
    organism_to_marker_to_evalue = defaultdict(dict)
    organism_to_marker_to_query = defaultdict(dict)
    markers = list()
    
    print(" * Compiling organism query markers from PyHMMSearch results: {}".format(opts.pyhmmsearch_results), file=sys.stderr)
    if opts.pyhmmsearch_results.endswith(".gz"):
        f_in = gzip.open(opts.pyhmmsearch_results, "rt")
    else:
        f_in = open(opts.pyhmmsearch_results, "r")

    if not opts.no_header:
        next(f_in)
    
    for line in tqdm(f_in):
        line = line.strip()
        if line:
            fields = line.split("\t")
            id_organism, id_query = fields[0].split(opts.sep)
            id_marker = fields[1]
            evalue = eval(fields[7])

            if id_marker not in organism_to_marker_to_query[id_organism]:
                print(id_organism, "|", id_marker, "-->", id_query, "| E-Value = {:.1e}".format(evalue), file=sys.stdout)
                organism_to_marker_to_query[id_organism][id_marker] = id_query
                organism_to_marker_to_evalue[id_organism][id_marker] = evalue
            else:
                if evalue < organism_to_marker_to_evalue[id_organism][id_marker]:
                    print(id_organism, "|", id_marker, "-->", id_query, "| E-Value = {:.1e}".format(evalue), "[Update] [Overriding: {}]".format(organism_to_marker_to_query[id_organism][id_marker]), file=sys.stdout)
                    organism_to_marker_to_query[id_organism][id_marker] = id_query
                    organism_to_marker_to_evalue[id_organism][id_marker] = evalue
            markers.append(id_marker)
    markers = set(markers)

    f_in.close()
   

    print(" * Writing marker list: {}".format(os.path.join(opts.output_directory, "markers.tsv")), file=sys.stderr)
    marker_to_organisms = defaultdict(list)
    marker_to_queries = defaultdict(list)
    for id_organism, marker_to_query in organism_to_marker_to_query.items():
        for id_marker, id_query in marker_to_query.items():
            marker_to_organisms[id_marker].append(id_organism)
            marker_to_queries[id_marker].append(id_query)
    df_markers = pd.concat([ 
        pd.Series(marker_to_organisms).to_frame("organisms"),
        pd.Series(marker_to_queries).to_frame("proteins"),
    ], axis=1)
    df_markers.to_csv(os.path.join(opts.output_directory, "markers.tsv"), sep="\t", header=None)

    # Create file objects
    print(" * Creating fasta file objects for N = {} markers".format(len(markers)), file=sys.stderr)
    files = dict()
    if opts.gzip:
        for id_marker in markers:
            files[id_marker] = gzip.open(os.path.join(opts.output_directory, "{}.faa.gz".format(id_marker)), "wt")
    else:
        for id_marker in markers:
            files[id_marker] = open(os.path.join(opts.output_directory, "{}.faa".format(id_marker)), "w")

    # Protein sequences 
    print(" * Getting query protein sequences: {}".format(opts.proteins), file=sys.stderr)
    # organism_to_query_to_seq = defaultdict(dict) 
    sequences = dict()
    if opts.proteins.endswith(".gz"):
        f_proteins = gzip.open(opts.proteins, "rt")
    else:
        f_proteins = open(opts.proteins, "r")
    for header, seq in SimpleFastaParser(f_proteins):
        header = header.split(" ")[0]
        id_organism, id_query = header.split(opts.sep)
        if id_query in organism_to_marker_to_query[id_organism].values():
            sequences[(id_organism, id_query)] = seq
            # organism_to_query_to_seq[id_organism][id_query] = seq 
    f_proteins.close()

    # Write sequences 
    print(" * Writing protein sequences for each marker: {}".format(opts.output_directory), file=sys.stderr)
    for id_organism, marker_to_query in organism_to_marker_to_query.items():
        for id_marker, id_query in marker_to_query.items():
            # seq = organism_to_query_to_seq[id_organism][id_query]
            seq = sequences[(id_organism, id_query)]
            header = "{}  {} {}".format(id_organism, id_marker, id_query)
            print(">{}\n{}".format(header,seq), file=files[id_marker])

    # Close files
    for id_marker, f in files.items():
        f.close()

if __name__ == "__main__":
    main()
    
                

