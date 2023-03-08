#!/usr/bin/env python
import sys, os, glob, argparse 
from collections import OrderedDict
import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.3.7"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <binning_directory> -m <mapping_directory> -c <clusters.tsv> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-m","--mapping_directory", type=str, help = "path/to/mapping_directory")
    parser.add_argument("-o","--output", type=str,  help = "path/to/output.tsv [Default: stdout]", default="stdout")
    parser.add_argument("-d","--sep", type=str, default="\t", help = "Delimiter for table [Default: <tab>]")
    parser.add_argument("--sample_column_label", type=str, default="id_sample", help = "Sample column label [Default: id_sample]")
    parser.add_argument("--feature_row_label", type=str, default="id_feature", help = "Feature row label [Default: id_feature]")
    parser.add_argument("--allow_missing_values", action="store_true", help = "Allow missing values instead of filling with zeros")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    if opts.output == "stdout":
        opts.output = sys.stdout

    # Merge feature counts
    counts = dict()
    for fp in tqdm(glob.glob(os.path.join(opts.mapping_directory, "*", "output", "counts.tsv.gz")), "Reading processed featureCounts outputs"):
        id_sample = fp.split("/")[-3]
        counts[id_sample] = pd.read_csv(fp, sep="\t", index_col=0).iloc[:,-1]
    X = pd.DataFrame(counts).T
    if not opts.allow_missing_values:
        X = X.fillna(0).astype(int)
        
    X.index.name = opts.sample_column_label
    X.columns.name = opts.feature_row_label
    X.to_csv(opts.output, sep=opts.sep)

    print("The following directory:{} includes n={} samples and m={} features in the concatenated table.".format(opts.mapping_directory, *X.shape), file=sys.stderr)

if __name__ == "__main__":
    main()
    
                

