#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
import numpy as np
import pandas as pd
from tqdm import tqdm 

pd.options.display.max_colwidth = 100
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.5.18"


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}

    Create input using the following commands: 
    ```
    cut -f1,3 metaeuk_output/genomes/identifier_mapping.tsv > metaeuk_output/genomes/protein_to_genome.tsv
    concatenate_dataframes.py -n -a 1 metaeuk_output/genomes/protein_to_genome.tsv mmseqs2_output/output/clusters.tsv | cut_table_by_column_index.py -f2,1,3 > input.tsv
    ```
    """.format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <input.tsv> -o <output.tsv>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline

    parser.add_argument("-i","--input",  type=str, default="stdin", help = "path/to/input.tsv [id_genome]<tab>[id_protein]<tab>[id_protein-cluster](No header) [Default: stdin]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")
    parser.add_argument("-b","--boolean", action="store_true", help = "Return True/False instead of integer counts")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # I/O
    if opts.input == "stdin":
        opts.input = sys.stdin 
    if opts.output == "stdout":
        opts.output = sys.stdout

    # Read Input
    df_input = pd.read_csv(opts.input, sep="\t", index_col=None, header=None)
    genomes = sorted(df_input.iloc[:,0].unique())
    clusters = sorted(df_input.iloc[:,2].unique())

    # Create array
    A = np.zeros((len(genomes), len(clusters)), dtype=int)

    for _, (id_genome, id_protein, id_cluster) in tqdm(df_input.iterrows(), total=df_input.shape[0]):
        i = genomes.index(id_genome)
        j = clusters.index(id_cluster)
        A[i,j] += 1

    # Create output
    df_output = pd.DataFrame(A, index=genomes, columns=clusters)
    df_output.index.name = "id_genome"
    df_output.columns.name = "id_protein-cluster"

    if opts.boolean:
        df_output = df_output > 0

    df_output.to_csv(opts.output, sep="\t")






  

if __name__ == "__main__":
    main()
