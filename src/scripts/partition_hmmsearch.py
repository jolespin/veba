#!/usr/bin/env python
import sys, os, glob, argparse, gzip
from collections import defaultdict
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.10.01"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <hmmsearch_tblout> -a <proteins> -o <output_directory> -s '|--|'".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--hmmsearch_tblout", type=str, help = "path/to/hmmsearch_tblout.tsv", required=True)
    parser.add_argument("-a","--proteins", type=str, help = "path/to/proteins.faa",  required=True)
    parser.add_argument("-o","--output_directory", type=str,  default="input", help = "Output directory [Default: Directory of --hmmsearch_tblout")
    parser.add_argument("-s", "--sep", default="|--|", type=str, help="[id_organism]<sep>[id_protein] [Default: |--|]")
    parser.add_argument("-f", "--hmm_marker_field", default="accession", type=str, help="HMM reference type (accession, name) [Default: accession")



    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    if opts.output_directory == "input":
        opts.output_directory = os.path.split(opts.hmmsearch_tblout)[0]
    
    os.makedirs(opts.output_directory, exist_ok=True)
    assert opts.hmm_marker_field in {"accession", "name"}, "--hmm_marker_field must be either accession or name"
    hmm_indexer = None
    if opts.hmm_marker_field == "accession":
        hmm_indexer = 3
    if opts.hmm_marker_field == "name":
        hmm_indexer = 2
    assert hmm_indexer is not None, "Please use --hmm_marker_field to specify HMM indexer"

    # Get marker to query protein
    print(" * Parsing HMMSearch tblout to get mapping of organisms and query markers: {}".format(opts.hmmsearch_tblout), file=sys.stderr)
    organism_to_marker_to_evalue = defaultdict(dict)
    organism_to_marker_to_query = defaultdict(dict)
    markers = list()
    with open(opts.hmmsearch_tblout, "r") as f:
        for line in f.readlines():
            line = line.strip()
            if not line.startswith("#"):
                fields = list(filter(bool, line.split(" ")))
                id_organism, id_query = fields[0].split(opts.sep)

                id_marker = fields[hmm_indexer]
                assert id_marker != "-", "Cannot use --hmm_marker_field {} with current HMM output because missing identifiers".format(opts.hmmer_marker_field)
                evalue = eval(fields[4])

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
    for id_marker in markers:
        files[id_marker] = gzip.open(os.path.join(opts.output_directory, "{}.faa.gz".format(id_marker)), "wt")

    
    # Protein sequences 
    print(" * Getting query protein sequences: {}".format(opts.proteins), file=sys.stderr)
    # organism_to_query_to_seq = defaultdict(dict) 
    sequences = dict()
    with open(opts.proteins, "r") as f:
        for header, seq in SimpleFastaParser(f):
            header = header.split(" ")[0]
            id_organism, id_query = header.split(opts.sep)
            if id_query in organism_to_marker_to_query[id_organism].values():
                sequences[(id_organism, id_query)] = seq
                # organism_to_query_to_seq[id_organism][id_query] = seq 

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
    
                

