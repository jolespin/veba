#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
from collections import defaultdict
# import numpy as np
import pandas as pd
from tqdm import tqdm 

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.6.7"


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <annotation_results.tsv[.gz]> -l genome -o <output_table>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline

    parser.add_argument("-i","--annotation_results",  type=str, default = "stdin", help = "path/to/annotation.tsv from annotate.py [Default: stdin]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/genomes_table.tsv [Default: stdout]")
    parser.add_argument("-l","--level", type=str, default="genome", help = "level {genome, genome_cluster} [Default: genome]")
    parser.add_argument("-p", "--include_protein_identifiers", action="store_true", help = "Write protein identifiers")
    parser.add_argument("--header", action="store_true", help = "Write header")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert opts.level in {"genome", "genome_cluster"}, "--level must be either {genome, genome_cluster}"

    if opts.level == "genome":
        level_field = ("Identifiers", "id_genome")
    if opts.level == "genome_cluster":
        level_field = ("Identifiers", "id_genome_cluster")

    if opts.annotation_results == "stdin":
        opts.annotation_results = sys.stdin 
    if opts.output == "stdout":
        opts.output = sys.stdout 
    
    df_annotations = pd.read_csv(opts.annotation_results, sep="\t", index_col=0, header=[0,1])

    output = list()
    for id_protein, (id_organism, ko_ids) in tqdm(df_annotations.loc[:,[level_field, ("KOfam", "ids")]].iterrows(), "Compiling KO identifiers", total=df_annotations.shape[0]):
        ko_ids = eval(ko_ids)
        if len(ko_ids):
            for id_ko in ko_ids:
                output.append([id_protein, id_organism, id_ko])
    df_output = pd.DataFrame(output, columns=["id_protein", level_field[1], "id_kegg-ortholog"])

    if not opts.include_protein_identifiers:
        df_output = df_output.drop("id_protein", axis=1)

    df_output.to_csv(opts.output, sep="\t", header=bool(opts.header), index=False)

if __name__ == "__main__":
    main()
