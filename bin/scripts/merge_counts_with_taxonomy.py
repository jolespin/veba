#!/usr/bin/env python
import sys, os, glob, argparse, warnings
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.3.8"

def main(argv=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -a <axis> -d <delimiter> -o <output.tsv> [table_1] [table_2]".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-X","--counts",  help = "path/to/counts.tsv [rows=samples, columns=features, tab-delimited]")
    parser.add_argument("-t","--taxonomy", type=str, default="stdin", help = "path/to/merged_taxonomy.tsv [id_genome/cluster]<tab>[classification] or [id_genome/cluster]<tab>[organism_type]<tab>[classification].  Recommended to use `merge_taxonomy_classifications.py`.  If more than 3 columns are found, only the first 3 will be considered.")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")
    parser.add_argument("-c","--clusters", type=str,  help = "path/to/clusters.tsv [id_genome]<tab>[id_genome_cluster], No header.  Not really applicable when --counts are SLC-level counts.")
    parser.add_argument("-T", "--transpose", action="store_true", help = "Transpose output dataframe [Default: rows=features, columns=samples]")
    parser.add_argument("--no_header_taxonomy", action="store_true", help = "No header on taxonomy table")
    parser.add_argument("--assert_features", action="store_true", help = "Assert that all features are in both tables")
    parser.add_argument("--fillna", type=float, default=0, help = "Fill missing values [Default: 0]")

    # Options
    opts = parser.parse_args(argv)
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Open output file
    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Load data
    X = pd.read_csv(opts.counts, sep="\t", index_col=0)
    if opts.no_header_taxonomy:
        df_taxonomy = pd.read_csv(opts.taxonomy, sep="\t", index_col=0, header=None)
    else:
        df_taxonomy = pd.read_csv(opts.taxonomy, sep="\t", index_col=0)

    A = set(X.columns)
    B = set(df_taxonomy.index)
    if opts.assert_features:
        assert A == B, "--counts columns must have the same values (not order) as --taxonomy rows.\nThe following are unique to --counts: {}\n\nThe following are unique to --taxonomy:{}".format(", ".join(sorted(A - B)), ", ".join(sorted(B - A)))
    else:
        if A != B:
            warnings.warn("--counts columns don't have the same values (not order) as --taxonomy rows.\nThe following are unique to --counts: {}\n\nThe following are unique to --taxonomy:{}".format(", ".join(sorted(A - B)), ", ".join(sorted(B - A))))
        assert B >= A, "All features in --taxonomy must be in --counts.\nMissing the following: {}".format(", ".join(sorted(A - B)))

    # Merge
    df_output = X.T.copy()

    n, m = df_taxonomy.shape
    if m > 2:
        df_taxonomy = df_taxonomy.iloc[:,:2]
        m = 2
    feature_to_taxonomy = df_taxonomy.iloc[:,-1]

    # assert m in {1,2}, "--taxonomy should have the following columns: [id_genome/cluster]<tab>[classification] or [id_genome/cluster]<tab>[organism_type]<tab>[classification]"
    if m == 1:
        df_output.insert(0, "taxonomy", feature_to_taxonomy)
    if m == 2:
        feature_to_domain = df_taxonomy.iloc[:,-2]
        df_output.insert(0, "organism_type", feature_to_domain)
        df_output.insert(1, "taxonomy", feature_to_taxonomy)

    if opts.clusters:
        feature_to_cluster = pd.read_csv(opts.clusters, sep="\t", index_col=0, header=None).iloc[:,0]
        assert set(feature_to_cluster.index) <= set(df_output.index), "All of the genomes in --clusters must be included in --counts"
        df_output.index = df_output.index.map(lambda x: (feature_to_cluster[x], x))

    # Output
    if opts.transpose: 
        df_output = df_output.T
    df_output.to_csv(opts.output, sep="\t")


if __name__ == "__main__":
    main()
    
                

