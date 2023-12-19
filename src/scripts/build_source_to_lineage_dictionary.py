#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, gzip, pickle
from tqdm import tqdm
import pandas as pd

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.11.13"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__

    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <source_lineage.tsv[.gz]> -o <output.dict.pkl[.gz]>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", default="stdin", type=str, help = "Path to table [id_source]<tab>[class]<tab>[order]<tab>[family]<tab>[genus]<tab>[species], with header.  Can include more columns but the first column must be `id_source`. [Default: stdin]")
    parser.add_argument("-o","--output", required=True, type=str, help = "Path to dictionary pickle object.  Can be gzipped. (Recommended name: source_to_lineage.dict.pkl.gz)")
    parser.add_argument("--separator", default=";", type=str, help = "Separator field for taxonomy [Default: ; ]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Input 
    if opts.input == "stdin":
        opts.input = sys.stdin 


    print(" * Reading identifier mappings from the following file: {}".format(opts.input), file=sys.stderr)
    source_to_lineage = dict()
    df_input = pd.read_csv(opts.input, sep="\t", index_col=0)
    for id_source, row in tqdm(df_input.loc[:,["class", "order", "family", "genus", "species"]].iterrows(), total=df_input.shape[0]):
        lineage = list()
        for level, taxon in row.items():
            v = level[0] + "__"
            if pd.notnull(taxon):
                v += taxon 
            lineage.append(v)
        source_to_lineage[id_source] = opts.separator.join(lineage)



    print(" * Writing Python dictionary: {}".format(opts.output), file=sys.stderr)
    f_out = None 
    if opts.output.endswith((".gz", ".pgz")):
        f_out = gzip.open(opts.output, "wb")
    else:
        f_out = open(opts.output, "wb")
    assert f_out is not None, "Unrecognized file format: {}".format(opts.output)
    pickle.dump(source_to_lineage, f_out)


   




if __name__ == "__main__":
    main()
