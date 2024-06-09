#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
import pandas as pd



pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.03.08"




#
def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <fasta_file> --compression infer > <saf_file>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input_directory", required=True, type=str, help = "path/to/input_directory")
    parser.add_argument("-o","--output", default="stdout", type=str, help = "path/to/output.tsv")
    parser.add_argument("-m","--mode", type=str, default="paired", help = "Compile mode {paired, separate} [Default: paired]")
    parser.add_argument("-r", "--ratios", action="store_true", help = "Add ratios in addition to counts. This will create a multiindex header.")

    # parser.add_argument("--retain_extension", action="store_true", help = "Retain file extensions (e.g., .fastq.gz)")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Assertions
    assert opts.mode in {"paired", "separate"}, "--mode must be in {paired, separate}"

    # I/O
    if opts.output == "stdout":
        opts.output = sys.stdout

    # Data
    data = dict()
    for fp in glob.glob(os.path.join(opts.input_directory, "*", "output", "seqkit_stats.concatenated.tsv")):
        id_sample = fp.split("/")[-3]
        file_to_depth = pd.read_csv(fp, sep="\t", index_col=0)["num_seqs"]
        file_to_depth.index = ["unprocessed_1.fastq.gz", "unprocessed_2.fastq.gz"] + file_to_depth.index[2:].map(lambda x: "/".join(x.split("/")[-2:])).tolist()
        file_to_depth = file_to_depth[file_to_depth.index.map(lambda x: "singletons" not in x)]
        data[id_sample] = file_to_depth 

    df_counts = pd.DataFrame(data).T

    # Remove extension
    if opts.mode == "paired":
        df_counts = df_counts.loc[:,df_counts.columns.map(lambda x: x.endswith("_1.fastq.gz"))]
        
    # Order the columns and clean up labels
    if opts.mode == "paired":
        df_counts.columns = df_counts.columns.map(lambda x: x[:-11])
        column_order = [
            "unprocessed",
            "1__fastp/trimmed", 
            "2__bowtie2/cleaned", "2__bowtie2/contaminated",
            "3__bbduk/non-kmer_hits", "3__bbduk/kmer_hits", 
        ]

    if opts.mode == "separate":
        df_counts.columns = df_counts.columns.map(lambda x: x[:-9])
        column_order = [
            "unprocessed_1", "unprocessed_2", 
            "1__fastp/trimmed_1", "1__fastp/trimmed_2", 
            "2__bowtie2/cleaned_1", "2__bowtie2/cleaned_2", "2__bowtie2/contaminated_1", "2__bowtie2/contaminated_2",
            "3__bbduk/non-kmer_hits_1", "3__bbduk/non-kmer_hits_2", "3__bbduk/kmer_hits_1", "3__bbduk/kmer_hits_2", 
        ]
    df_counts = df_counts.loc[:,column_order]


    # Output
    if not opts.ratios:
        df_counts.index.name = "SampleID"
        df_counts.to_csv(opts.output, sep="\t")

    if opts.ratios:
        assert opts.mode == "paired", "--ratios only implemented for --mode paired"

        df_ratios = pd.DataFrame()
        df_ratios["1__fastp/trimmed"] = df_counts["1__fastp/trimmed"]/df_counts["unprocessed"]
        df_ratios["2__bowtie2/cleaned"] = df_counts["2__bowtie2/cleaned"]/df_counts["1__fastp/trimmed"]
        df_ratios["2__bowtie2/contaminated"] = df_counts["2__bowtie2/contaminated"]/df_counts["1__fastp/trimmed"]
        df_ratios["3__bbduk/non-kmer_hits"] = df_counts["3__bbduk/non-kmer_hits"]/df_counts["1__fastp/trimmed"]
        df_ratios["3__bbduk/kmer_hits"] = df_counts["3__bbduk/kmer_hits"]/df_counts["1__fastp/trimmed"]


        df_ratios.columns = df_ratios.columns.map(lambda x: ("Ratios", x))
        df_counts.columns = df_counts.columns.map(lambda x: ("Counts", x))

        df_output = pd.concat([df_counts, df_ratios], axis=1)
        df_output.index.name = "SampleID"
        df_output.to_csv(opts.output, sep="\t")


if __name__ == "__main__":
    main()
