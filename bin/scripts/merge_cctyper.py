#!/usr/bin/env python
import sys, os, glob, argparse, gzip
import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.3.1"

def main(argv=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <cctyper_output> -t <type> -o <output.tsv>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--cctyper_directory", required=True, type=str, help = "path/to/cctyper/ with subdirectories for each genome")
    parser.add_argument("-t","--type", type=str, choices={"cas_operons", "crispr_cas", "spacers"}, required=True, help = "Type of file to merge")
    parser.add_argument("-o","--output", type=str,  default="stdout", help = "Output merged file [Default: stdout]")

    # Options
    opts = parser.parse_args(argv)
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Open output file
    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Merge tables
    if opts.type in {"cas_operons", "crispr_cas"}:
        dataframes = list()

        if opts.type == "cas_operons":
            filename = "cas_operons.tab"
        if opts.type == "crispr_cas":
            filename = "CRISPR_Cas.tab"
        for fp in tqdm(glob.glob(os.path.join(opts.cctyper_directory, "*", filename)), "Reading {} files: {}".format(opts.type, filename)):
            id_genome = fp.split("/")[-2]
            df = pd.read_csv(fp, sep="\t", index_col=0)
            df.index = df.index.map(lambda x: (id_genome, x))
            dataframes.append(df)
        df_output = pd.concat(dataframes, axis=0)
        df_output.index.names = ["id_genome", "id_contig"]

        df_output.to_csv(opts.output, sep="\t")

    if opts.type == "spacers":
        from Bio.SeqIO.FastaIO import SimpleFastaParser
        filename = "*.fa"
        if opts.output != sys.stdout:
            if opts.output.endswith(".gz"):
                opts.output = gzip.open(opts.output, "wt")
            else:
                opts.output = open(opts.output, "w")
        for fp in tqdm(glob.glob(os.path.join(opts.cctyper_directory, "*", "spacers", filename)), "Reading {} files: {}".format(opts.type, filename)):
            id_genome = fp.split("/")[-3]

            with open(fp, "r") as f_in:
                for id, seq in SimpleFastaParser(f_in):
                    header = "{} {}".format(id.split(" ")[0].strip(), id_genome)
                    print(">{}\n{}".format(header, seq), file=opts.output)

        if opts.output != sys.stdout:
            opts.output.close()

if __name__ == "__main__":
    main()