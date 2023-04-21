#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, shutil, time
from multiprocessing import cpu_count
from collections import OrderedDict, defaultdict
import pandas as pd

# Soothsayer Ecosystem
from soothsayer_utils import *

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.3.17"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <sequences.fasta> -c <clusters.tsv> -o <output>".format(__program__)

    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i", "--fasta", type=str, required=True, help = "path/to/representative_sequences.fasta (e.g., mmseqs2_rep_seqs.fasta)")
    parser.add_argument("-c", "--clusters", type=str, required=True, help = "path/to/clusters.tsv, Format: [id_sequence]<tab>[id_cluster] (No header)")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output [Default: stdout]")
    parser.add_argument("-f","--output_format", type=str, default="table", help = "Format of output: {table, fasta} [Default: table]")
    parser.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
   
    # Options
    opts = parser.parse_args()

    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    if opts.output == "stdout":
        opts.output = sys.stdout
    
    assert_acceptable_arguments(opts.output_format, {"table", "fasta"})

    id_to_cluster = pd.read_csv(opts.clusters, sep="\t", index_col=0, header=None).iloc[:,0]
    representative_sequences = read_fasta(opts.fasta, description=False, verbose=True)

    assert set(representative_sequences.index) <= set(id_to_cluster.index), "Not all representative sequences are in the clustering table"

    if opts.output_format == "table":
        df = id_to_cluster[representative_sequences.index].to_frame("id_cluster")
        df["representative_sequence"] = representative_sequences
        df.index.name = "id_representative_sequence"
        df = df.reset_index(drop=False).set_index("id_cluster")[["id_representative_sequence", "representative_sequence"]]
        df.to_csv(opts.output, sep="\t")

    if opts.output_format == "fasta":
        representative_sequences.index = representative_sequences.index.map(lambda id_sequence: "{} {}".format(id_to_cluster[id_sequence], id_sequence))
        write_fasta(representative_sequences, opts.output)


if __name__ == "__main__":
    main(sys.argv[1:])
