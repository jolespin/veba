#!/usr/bin/env python
import sys, os, argparse 
from collections import OrderedDict
import numpy as np
import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.10.26"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -o <output_directory> [file_1.tsv, file_2.tsv.gz, ..., file_3.tsv]".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("files", nargs="+", help = "Counts tables of 2 columns where 1st column is the feature and second column is the value")
    parser.add_argument("-o","--output", type=str,  help = "path/to/output.tsv [Default: stdout]", default="stdout")
    parser.add_argument("-d","--sep", type=str, default="\t", help = "Delimiter for table [Default: <tab>]")
    parser.add_argument("-c", "--sample_column_label", type=str, default="id_sample", help = "Sample column label [Default: id_sample]")
    parser.add_argument("-r", "--feature_row_label", type=str, default="id_feature", help = "Feature row label [Default: id_feature]")
    parser.add_argument("-i", "--sample_index", type=int,default=-3, help = "Sample index in filepath (e.g., veba_output/mapping/global/[ID]/output/counts.tsv it would be -3) [Default: -3]")
    parser.add_argument("-a", "--allow_missing_values", action="store_true", help = "Allow missing values instead of filling with zeros")
    parser.add_argument("-e", "--remove_empty_features", action="store_true", help = "Remove empty features")
    parser.add_argument("--pickle", type=str, help = "path/to/pandas.pkl output")
    parser.add_argument("--comment", type=str, default="#", help = "Comments prefixed with this character [Default: #]")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    if opts.output == "stdout":
        opts.output = sys.stdout

    # Merge feature counts
    output = dict()
    for fp in tqdm(opts.files, "Reading processed featureCounts outputs"):
        id_sample = fp.split("/")[opts.sample_index]
        counts = pd.read_csv(fp, sep="\t", index_col=0, header=None, comment=opts.comment).iloc[:,-1]
        if opts.remove_empty_features:
            counts = counts[counts > 0]
        output[id_sample] = counts
    X = pd.DataFrame(output).T

    if not opts.allow_missing_values:
        X = X.fillna(0)

    if np.allclose(X, X.astype(int), rtol=1e-05, atol=1e-08, equal_nan=True):
        X = X.astype(int)
        
    X.index.name = opts.sample_column_label
    X.columns.name = opts.feature_row_label
    X.to_csv(opts.output, sep=opts.sep)

    if opts.pickle:
        X.to_pickle(opts.pickle)

    sparsity = (X.values == 0).ravel().mean() * 100
    print("There are n={} samples and m={} features in the concatenated output table ({}% sparse).".format(*X.shape, sparsity), file=sys.stderr)
    
if __name__ == "__main__":
    main()
    
                

