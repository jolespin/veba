#!/usr/bin/env python
import sys, os, glob, argparse 
from collections import OrderedDict
import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.12.12"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -m <mapping_directory> -c <clusters.tsv> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-m","--mapping_directory", type=str, help = "path/to/mapping_directory")
    parser.add_argument("-o","--output", type=str, default="stdout",  help = "path/to/output.tsv [Default: stdout]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Merge genome spatial coverage
    coverage = dict()
    for fp in tqdm(glob.glob(os.path.join(opts.mapping_directory, "*", "output", "genome_spatial_coverage.tsv.gz")), "Reading genome spatial coverage tables"):
        id_sample = fp.split("/")[-3]
        coverage[id_sample] = pd.read_csv(fp, sep="\t", index_col=0).iloc[:,-1]

    X_coverage = pd.DataFrame(coverage).T
    X_coverage.index.name = "id_sample"
    X_coverage.columns.name = "id_genome"

    X_coverage.to_csv(opts.output, sep="\t")


if __name__ == "__main__":
    main()
    
                

