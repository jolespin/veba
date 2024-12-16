#!/usr/bin/env python
import sys, os, argparse, gzip
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.11.15"


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <diamond_output.tsv[.gz]> -l <label> -a <annotations.tsv[.gz] -o <output.tsv[.gz]>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", type=str, default="stdin", help = "path/to/diamond_output.tsv[.gz] with header [Default: stdin]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/annotations.MODIFIED.tsv[.gz] [Default: stdout]")
    parser.add_argument("-a","--annotations", type=str, required=True, help = "path/to/annotations.tsv[.gz]")
    parser.add_argument("-l","--label", type=str, required=True,  help = "Label to add for level 0 (first row) in VEBA multi-header (e.g., PlasticDB)")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Input
    if opts.input == "stdin":
        opts.input = sys.stdin 


    # Output
    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Annotations
    df_annotation = pd.read_csv(opts.annotations, sep="\t", index_col=0, header=[0,1])
    index_name = df_annotation.index.name
    
    # Input
    df_input = pd.read_csv(opts.input, sep="\t", index_col=0)
    df_input.columns = df_input.columns.map(lambda x: (opts.label, x))
    if not set(df_input.index) & set(df_annotation.index):
        raise ValueError("No overlap in identifiers between input and annotation files.")
    df_input = df_input.reindex(df_annotation.index)
    
    # Output
    df_annotation = pd.concat([df_annotation, df_input], axis=1)
    df_annotation.index.name = index_name
    df_annotation.to_csv(opts.output, sep="\t")


if __name__ == "__main__":
    main()
    
            