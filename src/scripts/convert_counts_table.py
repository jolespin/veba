#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, warnings
import pandas as pd
from soothsayer_utils import check_packages, assert_acceptable_arguments

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.5.28"


@check_packages(["anndata"])
def table_to_anndata(X, df_samples, df_features, output, compression):
    import anndata as ad
    if compression:
        compression = "gzip"

    adata = ad.AnnData(
        X=X, 
        obs=df_samples, 
        var=df_features,
    )
    adata.raw = adata
    adata.write(filename=output, compression=compression)

@check_packages(["biom"])
def table_to_biom(X, df_samples, df_features, output, compression):
    from biom.util import biom_open  
    from biom.table import Table

    table = Table(
    # Counts
    data=X.T.values, 
    # Samples
    sample_ids=X.index, 
    sample_metadata=list(map(lambda x: x[1].to_dict(), df_samples.iterrows())) if df_samples is not None else None,
    # Features
    observation_ids=X.columns, 
    observation_metadata=list(map(lambda x: x[1].to_dict(), df_features.iterrows())) if df_features is not None else None,
    )
    
    with biom_open(output, "w") as f:  
        table.to_hdf5(h5grp=f, generated_by=X.index.name, compress=bool(compression))

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <table.tsv> -f <format> -o <output>|-n sample_metadata.tsv -m feature_metadata.tsv".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", required=True, default="stdin", type=str, help = "path/to/input.tsv[.gz] [Default: stdin]")
    parser.add_argument("-f","--format", required=True, type=str, help = "{pickle, anndata, biom}")
    parser.add_argument("-o","--output", required=True, type=str, help = "Recommended extensions: {pickle:pkl, anndata:h5ad, biom:biom}")
    parser.add_argument("-n", "--sample_metadata", type=str, help = "Sample metadata")
    parser.add_argument("-m", "--feature_metadata", type=str, help = "Feature metadata")
    parser.add_argument("-d", "--sep", type=str, default="\t", help = "Input delimiter [Default: <tab>]")
    parser.add_argument("-c", "--compress", action="store_true", help = "Gzip compress output")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert_acceptable_arguments(opts.format, {"pickle", "anndata", "biom"})

    # Input
    index_name = opts.input
    if opts.input == "stdin":
        opts.input = sys.stdin

    df_counts = pd.read_csv(opts.input, sep=opts.sep, index_col=0, comment="#")
    df_counts.index.name = index_name
    df_samples = None
    df_features = None 

    if opts.sample_metadata is not None:
        df_samples = pd.read_csv(opts.sample_metadata, sep=opts.sep, index_col=0, comment="#")
        assert set(df_samples.index) >= set(df_counts.index), "--sample_metadata is missing the following samples that are in --input:\n{}".format(*sorted(set(df_counts.index)- set(df_samples.index)), sep="\n")
        df_samples = df_samples.loc[df_counts.index]

    if opts.feature_metadata is not None:
        df_features = pd.read_csv(opts.feature_metadata, sep=opts.sep, index_col=0, comment="#")
        assert set(df_features.index) >= set(df_counts.columns), "--feature_metadata is missing the following features that are in --input:\n{}".format(*sorted(set(df_counts.columns)- set(df_features.index)), sep="\n")
        df_features = df_features.loc[df_counts.columns]

    print("Converting --input {} to --output {} [Format: {}]".format(opts.input, opts.output, opts.format), file=sys.stderr)

    if opts.format not in {"pickle"}:
        if opts.sample_metadata is None:
            warnings.warn("--sample_metadata is not provided")
        if opts.feature_metadata is  None:
            warnings.warn("--feature_metadata is not provided")

    # Output
    if opts.format == "pickle":
    
        if opts.sample_metadata is not None:
            warnings.warn("--sample_metadata is not used with pickle output format")
        if opts.feature_metadata is not None:
            warnings.warn("--feature_metadata is not used with pickle output format")
        df_counts.to_pickle(opts.output, compression="gzip" if opts.compress else "infer")

    if opts.format == "anndata":
        table_to_anndata(X=df_counts, df_samples=df_samples, df_features=df_features, output=opts.output, compression=opts.compress)
    
    if opts.format == "biom":
        table_to_biom(X=df_counts, df_samples=df_samples, df_features=df_features, output=opts.output, compression=opts.compress)





if __name__ == "__main__":
    main()
