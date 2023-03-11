#!/usr/bin/env python
import sys, os, argparse
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.3.9"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <features> -s <synopsis> -d <diamond> -o <output.tsv[.gz]>".format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i","--features", type=str, required=True, help = "path/to/biosynthetic_gene_clusters.features.tsv.gz")
    parser_io.add_argument("-s","--synopsis", type=str, required=True, help = "path/to/biosynthetic_gene_clusters.synopsis.tsv.gz")
    parser_io.add_argument("-d","--diamond", type=str, required=True, help = "path/to/homology.mibig.tsv.gz")
    parser_io.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")

    parser_thresholds = parser.add_argument_group('Optional thresholds arguments')
    parser.add_argument("--pident", type=float, default=0.0, help = "pident lower bound [float:0 ≤ x < 100] [Default: 0]")
    parser.add_argument("--qcovhsp", type=float, default=0.0, help = "qcovhsp lower bound [float:0 ≤ x < 100] [Default: 0]")
    parser.add_argument("--scovhsp", type=float, default=0.0, help = "scovhsp lower bound [float:0 ≤ x < 100] [Default: 0]")
    parser.add_argument("--evalue", type=float, default=1e-3, help = "e-value lower bound [float:0 < x < 1] [Default: 1e-3]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Output
    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Load tables
    df_features = pd.read_csv(opts.features, sep="\t", index_col=None)
    df_synopsis = pd.read_csv(opts.synopsis, sep="\t", index_col=0)
    df_diamond = pd.read_csv(opts.diamond, sep="\t", index_col=0)

    # Gene -> BGC
    gene_to_bgc = dict(zip(df_features["gene_id"], df_features["bgc_id"]))

    # Filter diamond 
    df_diamond = df_diamond.query("pident >= {}".format(opts.pident))
    df_diamond = df_diamond.query("qcovhsp >= {}".format(opts.qcovhsp))
    df_diamond = df_diamond.query("scovhsp >= {}".format(opts.scovhsp))
    df_diamond = df_diamond.query("evalue <= {}".format(opts.evalue))

    # BGC -> Number of hits
    bgc_to_nhits = pd.Series([0]*df_synopsis.shape[0], df_synopsis.index)
    bgc_to_nhits.update(df_diamond.index.map(lambda x: gene_to_bgc[x]).value_counts())
    df_synopsis["number_of_mibig_hits"] = bgc_to_nhits
    df_synopsis["novelty_score"] = 1 - df_synopsis["number_of_mibig_hits"]/df_synopsis["number_of_genes"]

    # Output
    df_synopsis.to_csv(opts.output, sep="\t")

if __name__ == "__main__":
    main()
    
                
