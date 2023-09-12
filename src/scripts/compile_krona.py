#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
import numpy as np
import pandas as pd

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.9.5"


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <input.tsv> -m <mode> > <output_table>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline

    parser.add_argument("-i","--input",  type=str, default="stdin", help = "path/to/input.tsv [Default: stdin]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")
    parser.add_argument("-m","--mode", type=str, required=True, help = "{prokaryotic, eukaryotic, biosynthetic-global, biosynthetic-local}")
    parser.add_argument("-R", "--retain_rank_prefix", action="store_true", help = "Retain rank prefixes for modes {prokaryotic, eukaryotic} (e.g., d__, p__, c__)")
    parser.add_argument("-u", "--unclassified_label", default="Unclassified", type=str, help = "Unclassified label [Default: Unclassified]")

    # parser.add_argument("-C", "--remove_incomplete", action="store_true", help = "Remove BGCs on contig edge for biosynthetic mode")
    # parser.add_argument("-G", "--remove_genome_column", action="store_true", help = "Remove the genome ID for biosynthetic mode.  Useful when using lots of genomes.")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert opts.mode in {"prokaryotic", "eukaryotic", "biosynthetic-global", "biosynthetic-local"}

    if opts.input == "stdin":
        opts.input = sys.stdin   

    if opts.output == "stdout":
        opts.output = sys.stdout 

    df_input = pd.read_csv(opts.input, sep="\t", index_col=0)
    if opts.mode == "prokaryotic":
        value_counts = df_input["classification"].value_counts()
        try:
            df_output = pd.DataFrame(np.stack(value_counts.index.map(lambda classification: list(classification.split(";")))))
        except ValueError:
            sizes = value_counts.index.map(lambda classification: len(classification.split(";")))
            index = list()
            for i,s in zip(value_counts.index, sizes):
                if s < 7:
                    i = ";".join(map(lambda x: f"{x}__{opts.unclassified_label}", list("dpcofgs")))
                else:
                    if s > 7:
                        print(f"Not sure why {i} has more than 7 fields.  Please double check this.", file=sys.stderr)
                index.append(i)
            value_counts.index = index
            df_output = pd.DataFrame(np.stack(value_counts.index.map(lambda classification: list(classification.split(";")))))

        if not opts.retain_rank_prefix:
            df_output = df_output.applymap(lambda x: x.split("__")[-1])
        df_output.columns = ["domain", "phylum","class", "order", "family", "genus", "species"]
        df_output.insert(0, "count", value_counts.values)

    if opts.mode == "eukaryotic":
        value_counts = df_input["consensus_classification"].value_counts()
        try:
            df_output = pd.DataFrame(np.stack(value_counts.index.map(lambda classification: list(classification.split(";")))))
        except ValueError:
            sizes = value_counts.index.map(lambda classification: len(classification.split(";")))
            index = list()
            for i,s in zip(value_counts.index, sizes):
                if s < 5:
                    i = ";".join(map(lambda x: f"{x}__{opts.unclassified_label}", list("cofgs")))
                else:
                    if s > 5:
                        print(f"Not sure why {i} has more than 5 fields.  Please double check this.", file=sys.stderr)
                index.append(i)
            value_counts.index = index
            df_output = pd.DataFrame(np.stack(value_counts.index.map(lambda classification: list(classification.split(";")))))
        if not opts.retain_rank_prefix:
            df_output = df_output.applymap(lambda x: x.split("__")[-1])
        df_output.columns = ["class", "order", "family", "genus", "species"]
        df_output.insert(0, "count", value_counts.values)

    if opts.mode == "biosynthetic-global":
        ## From synopsis
        # if opts.remove_incomplete:
        #     df_input = df_input.query("cluster_on_contig_edge == False")
        # gb = df_input.groupby(["sample_of_origin", "genome_id", "bgc_type"])
        # df_output = gb.size().to_frame("count").reset_index()
        # df_output = df_output.loc[:,["count", "sample_of_origin", "genome_id", "bgc_type"]]
        # if opts.remove_genome_column:
        #     df_output = df_output.drop("genome_id", axis=1)

        df_output = df_input.reset_index().loc[:,["number_of_bgcs", "bgc_type","id_genome"]]
        df_output.columns = [ "count", "bgc_type","id_genome"]

    if opts.mode == "biosynthetic-local":
        # if opts.remove_incomplete:
        #     df_input = df_input.query("cluster_on_contig_edge == False")
        # gb = df_input.groupby(["genome_id", "bgc_type"])
        # df_output = gb.size().to_frame("count").reset_index()
        # df_output = df_output.loc[:,["count", "genome_id", "bgc_type"]]
        # if opts.remove_genome_column:
        #     df_output = df_output.drop("genome_id", axis=1)

        df_output = df_input.reset_index().loc[:,[ "number_of_bgcs", "bgc_type"]]
        df_output.columns = ["count", "bgc_type"]

    df_output.sort_values("count", ascending=False).to_csv(opts.output, sep="\t", header=None, index=None)

if __name__ == "__main__":
    main()
