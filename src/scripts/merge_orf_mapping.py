#!/usr/bin/env python
import sys, os, glob, argparse 
from collections import OrderedDict
import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.03.26"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <binning_directory> -m <mapping_directory> -c <orthogroups.tsv> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--binning_directory", type=str, help = "path/to/binning_directory")
    parser.add_argument("-c","--orthogroups", type=str, help = "path/to/orthogroups")
    parser.add_argument("-m","--mapping_directory", type=str, help = "path/to/mapping_directory")
    parser.add_argument("-o","--output_directory", type=str, help = "path/to/output_directory [Default: counts_output]", default="counts_output")
    parser.add_argument("-f","--format", type=str, default="tsv", help = "Output format: {tsv, csv, pickle} [Future will support feather, parquet] [Default: tsv]")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    os.makedirs(opts.output_directory, exist_ok=True)


    # ORF -> ORTHOGROUP
    orf_to_orthogroup = pd.read_csv(opts.orthogroups, sep="\t", index_col=0, header=None).iloc[:,0]

    # Merge orf counts
    counts = dict()
    for fp in tqdm(glob.glob(os.path.join(opts.mapping_directory, "*", "output", "featurecounts.tsv.gz")), "Reading featureCounts output"):
        id_sample = fp.split("/")[-3]
        counts[id_sample] = pd.read_csv(fp, sep="\t", index_col=0, comment="#").iloc[:,-1]

    X_orfs = pd.DataFrame(counts).T
    X_orthogroups = X_orfs.groupby(orf_to_orthogroup, axis=1).sum() 

    if opts.format == "tsv":
        fp = os.path.join(opts.output_directory, "X_orfs.tsv.gz")
        print(f"Writing ORFs counts table: {fp}")
        X_orfs.to_csv(fp, sep="\t")

        fp = os.path.join(opts.output_directory, "X_orthogroups.tsv.gz")
        print(f"Writing orthogroups counts table: {fp}")
        X_orthogroups.to_csv(fp, sep="\t")



    if opts.format == "csv":
        fp = os.path.join(opts.output_directory, "X_orfs.csv.gz")
        print(f"Writing ORFs counts table: {fp}")
        X_orfs.to_csv(fp, sep=",")

        fp = os.path.join(opts.output_directory, "X_orthogroups.csv.gz")
        print(f"Writing orthogroups counts table: {fp}")
        X_orthogroups.to_csv(fp, sep=",")


    if opts.format == "pickle":
        fp = os.path.join(opts.output_directory, "X_orfs.pkl")
        print(f"Writing ORFs counts table: {fp}")
        X_orfs.to_pickle(fp)

        fp = os.path.join(opts.output_directory, "X_orthogroups.pkl")
        print(f"Writing orthogroups counts table: {fp}")
        X_maX_orthogroupsgs.to_pickle(fp)



    # if opts.format == "feather":
    #     fp = os.path.join(opts.output_directory, "X_contigs.feather")
    #     print(f"Writing contigs counts table: {fp}")
    #     X_contigs.to_feather(fp)

    #     fp = os.path.join(opts.output_directory, "X_mags.feather")
    #     print(f"Writing MAGs counts table: {fp}")
    #     X_mags.to_feather(fp)

    #     fp = os.path.join(opts.output_directory, "X_slcs.feather")
    #     print(f"Writing SLCs counts table: {fp}")
    #     X_slcs.to_feather(fp)

    # if opts.format == "parquet":
    #     fp = os.path.join(opts.output_directory, "X_contigs.parquet")
    #     print(f"Writing contigs counts table: {fp}")
    #     X_contigs.to_parquet(fp)

    #     fp = os.path.join(opts.output_directory, "X_mags.parquet")
    #     print(f"Writing MAGs counts table: {fp}")
    #     X_mags.to_parquet(fp)

    #     fp = os.path.join(opts.output_directory, "X_slcs.parquet")
    #     print(f"Writing SLCs counts table: {fp}")
    #     X_slcs.to_parquet(fp)


if __name__ == "__main__":
    main()
    
                

