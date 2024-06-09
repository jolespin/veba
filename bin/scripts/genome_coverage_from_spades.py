#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse
import numpy as np
import pandas as pd
from collections import defaultdict

# Soothsayer Ecosystem
from soothsayer_utils import get_file_object, pv

# Bioinformatics
from Bio.SeqIO.FastaIO import SimpleFastaParser

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.7.14"



def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <genomes>  -o <output>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-i","--genomes", required=True,nargs="+", type=str, help = "path/to/fasta.fa[.gz]")
    # parser.add_argument("-c","--coverage", default="stdin", type=str, help = "path/to/contig_coverage.tsv [id_contig]<tab>[coverage_value]")
    parser.add_argument("-x","--extension", type=str, help = "Remove extension")
    parser.add_argument("-o","--output", default="stdout", type=str, help = "path/to/output.tsv. Calculated via sum of (contig_length * contig_coverage) divided by sum of all contig lengths that belong to a genome")
    parser.add_argument("-f","--float", action="store_true",  help = "Output as float. Default is to round and add X")
    parser.add_argument("--include_full_path", action="store_true",  help = "Include full path instead of just basename")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # if opts.coverage == "stdin":
    #     opts.coverage = sys.stdin
    if opts.output == "stdout":
        opts.output = sys.stdout

    # Read fasta
    contig_to_coverage = dict()
    contig_to_size = dict()
    contig_to_genome = dict() 
    genome_to_size = defaultdict(int)
    for fp in pv(opts.genomes, "Parsing genome assemblies"):
        if opts.include_full_path:
            id_genome = fp 
        else:
            id_genome = fp.split("/")[-1]
            if opts.extension:
                id_genome = id_genome.split(opts.extension)[0]
        with get_file_object(fp, mode="read", verbose=False) as f:        
            for id, seq in SimpleFastaParser(f):
                id = id.split(" ")[0]
                assert id not in contig_to_genome, "Cannot have duplicate contig identifiers.  If you have duplicate contig ids, then run separately. First duplicate encountered [{}]".format(id)
                contig_size = len(seq)
                contig_to_size[id] = contig_size
                contig_to_genome[id] = id_genome 
                genome_to_size[id_genome] += contig_size
                fields = id.split("_")
                contig_to_coverage[id] = float(fields[fields.index("cov") + 1])
    contig_to_coverage = pd.Series(contig_to_coverage)
    contig_to_size = pd.Series(contig_to_size)
    contig_to_genome = pd.Series(contig_to_genome)
    genome_to_size = pd.Series(genome_to_size)

    # Coverage
    # Calculated via sum of (contig_length * contig_coverage) divided by sum of all contig lengths that belong to a genome
    # contig_to_coverage = pd.read_csv(opts.coverage, sep="\t", index_col=0, header=None).squeeze("columns")
    # assert set(contig_to_genome.index) <= set(contig_to_coverage.index), "Not all contigs from genomes have coverage values"
    genome_to_coverage = dict()
    for id_genome, group in pv(contig_to_genome.groupby(contig_to_genome), "Calculating coverage for genomes"):
        contigs = group.index
        genome_size = genome_to_size[id_genome]
        aggregated_coverage = np.sum(contig_to_size[contigs] * contig_to_coverage[contigs])
        genome_to_coverage[id_genome] = aggregated_coverage/genome_size
    genome_to_coverage = pd.Series(genome_to_coverage)

    if not opts.float:
        genome_to_coverage = genome_to_coverage.map(lambda x: "{}X".format(round(x)))
    genome_to_coverage.to_frame().to_csv(opts.output, sep="\t", header=None)
        
if __name__ == "__main__":
    main()
