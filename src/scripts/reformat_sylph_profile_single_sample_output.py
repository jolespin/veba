#!/usr/bin/env python
import sys, os, argparse, gzip 
import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.11.10"

def filepath_to_genome(fp, extension):
    assert fp.endswith(extension)
    fn = os.path.split(fp)[1]
    return fn[:-(len(extension) + 1)]

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__

    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <input.fasta> -o <output.fasta>)".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-i","--input", default="stdin", type=str, help = "Input fasta file")
    parser.add_argument("-o","--output_directory", required=True, type=str, help = "Output directory to write output files")
    # parser.add_argument("-n", "--name", type=str, required=False, help="Name of sample")
    parser.add_argument("-c","--genome_clusters", type=str, help = "path/to/mags_to_slcs.tsv. [id_genome]<tab>[id_genome-cluster], No header.")
    parser.add_argument("-f","--field", type=str, default="Taxonomic_abundance", help = "Field to use for reformating [Default: Taxonomic_abundance]")
    parser.add_argument("-x","--extension", type=str, default="fa", help = "Fasta file extension for bins [Default: fa]")
    parser.add_argument("--header", action="store_true",  help = "Do not include header.  Doesn't apply to unstacked dataframe.")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Input
    if opts.input == "stdin":
        opts.input = sys.stdin 
    
    # Output
    os.makedirs(opts.output_directory, exist_ok=True)

    # Process
    df_sylph = pd.read_csv(opts.input, sep="\t")
    assert opts.field in df_sylph.columns, "--field {} not in --input columns: {}".format(opts.field, ", ".join(df_sylph.columns))

    genome_to_value = df_sylph.set_index("Genome_file")[opts.field]
    genome_to_value.index = genome_to_value.index.map(lambda fp: filepath_to_genome(fp, opts.extension))

    # Output genome values
    genome_to_value.to_frame(opts.field.lower()).to_csv(os.path.join(opts.output_directory, "{}.tsv.gz".format(opts.field.lower())), sep="\t", header=bool(opts.header))

    if opts.genome_clusters:
        genome_to_slc = pd.read_csv(opts.genome_clusters, sep="\t", index_col=0).iloc[:,0]
        slc_to_value = genome_to_value.groupby(genome_to_slc).sum()
        slc_to_value.to_frame(opts.field.lower()).to_csv(os.path.join(opts.output_directory, "{}.clusters.tsv.gz".format(opts.field.lower())), sep="\t", header=bool(opts.header))

if __name__ == "__main__":
    main()
    
                

