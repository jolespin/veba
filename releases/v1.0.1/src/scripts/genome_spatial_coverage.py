#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, gzip
from collections import OrderedDict, defaultdict
import pandas as pd
from tqdm import tqdm
from Bio.SeqIO.FastaIO import SimpleFastaParser

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.08.17"


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -f <reference.fasta[.gz]> -i <scaffolds_to_bins.tsv> -o <output> coverage_1.tsv[.gz] [coverage_2, ..., coverage_n]".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("coverage", type=str, nargs="+", help = "path/to/coverage[s]. One or more samtools coverage files of reads mapped to contigs.  Assumes mapped to same reference. [Required]")
    parser.add_argument("-i","--scaffolds_to_bins", type=str, required=True, help = "path/to/scaffolds_to_bins.tsv [id_scaffold]<tab>[id_bin] [Required]")
    parser.add_argument("-f","--fasta", type=str, required=True, help = "path/to/reference.fasta.  Must contain all contigs from --scaffolds_to_bins [Required]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output [Default: stdout]")
    parser.add_argument("-b","--basename",  action="store_true", help = "Output basename for multiple bam files. Equivalent to --index_split_position -1.  Cannot use with --index_split_position [Default: path]")
    parser.add_argument("-s", "--index_split_position",  type=int, help = "Filename index split (e.g., output_directory/[id_name]/mapped.sorted.bam.cov you would choose either 1 or -2 for [id_name].   Cannot use with basename.")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Assertion
    if opts.basename:
        assert opts.index_split_position is None, "Cannot use --basename and --index_split_position."
        opts.index_split_position = -1

    # Output
    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Read scaffolds to bins
    contig_to_mag = pd.read_csv(opts.scaffolds_to_bins, sep="\t", index_col=0, header=None).iloc[:,0]

    # Read fasta file
    if opts.fasta.endswith(".gz"):
        f_fasta = gzip.open(opts.fasta, "rt")
    else:
        f_fasta = open(opts.fasta, "r")

    contig_to_length = dict()
    for id, seq in tqdm(SimpleFastaParser(f_fasta), "Reading fasta file: {}".format(opts.fasta)):
        id = id.split(" ")[0]
        contig_to_length[id] = len(seq)
    f_fasta.close()

    # Contig lengths
    contig_to_length = pd.Series(contig_to_length)
    assert set(contig_to_mag.index) <= set(contig_to_length.index), "Please ensure all contigs from --scaffolds_to_bins are available in --fasta"

    # MAG lengths
    mag_to_length = contig_to_length.groupby(contig_to_mag).sum()

    # Read coverage files
    if len(opts.coverage) == 1:
        fp = opts.coverage[0]
        # Load in covered bases
        try:
            contig_to_coverage = pd.read_csv(fp, sep="\t", index_col=0)["covbases"]
        except UnicodeDecodeError as e:
            print("\nUnicodeDecodeError: {}\nDid you accidentally use a bam file ({}) instead of samtools coverage table?".format(e, fp), file=sys.stderr)
            sys.exit(1)
        # Aggregate w.r.t MAG
        mag_to_coverage = contig_to_coverage.groupby(contig_to_mag).sum()
        # Get ratio of covered bases w.r.t MAG
        spatial_coverage = mag_to_coverage/mag_to_length[mag_to_coverage.index]
        # Output
        spatial_coverage.to_csv(opts.output, sep="\t", header=None)
    else:
        # Load in covered bases
        id_to_cov = dict()
        for fp in tqdm(opts.coverage, "Reading coverage files"):
            if  opts.index_split_position is not None:
                id = fp.split("/")[opts.index_split_position]
            else:
                id = fp
            assert id not in id_to_cov, "{} is a duplicate.  Please ensure that there are no duplicates filenames.  Try excluding the --basename option.".format(id)
            try:
                id_to_cov[id] = pd.read_csv(fp, sep="\t", index_col=0)["covbases"]
            except UnicodeDecodeError as e:
                print("\nUnicodeDecodeError: {}\n\nDid you accidentally use a bam file instead of samtools coverage table?\n\t{}".format(e, fp), file=sys.stderr)
                sys.exit(1)
        df_contig_coverage = pd.DataFrame(id_to_cov).T

        # Aggregate w.r.t MAG
        df_mag_coverage = df_contig_coverage.groupby(contig_to_mag, axis=1).sum()

        # Get ratio of covered bases w.r.t MAG
        df_spatial_coverage = df_mag_coverage/mag_to_length[df_mag_coverage.columns].values

        # Output
        df_spatial_coverage.to_csv(opts.output, sep="\t")

if __name__ == "__main__":
    main()
