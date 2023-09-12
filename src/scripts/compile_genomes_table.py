#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
from collections import defaultdict
import numpy as np
import pandas as pd
from tqdm import tqdm 

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.8.28"


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <preprocess_directory> -b cleaned > <output_table>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline

    parser_io = parser.add_argument_group('[Mode 1] Preprocess Directory arguments')
    parser_io.add_argument("-i","--binning_directory",  type=str, help = "path/to/binning_directory (e.g., veba_output/binning)")
    parser_io.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")
    # parser_io.add_argument("-t","--organism_type", type=str, default="infer", help = "organism type [Default: infer]")
    # parser_io.add_argument("--genome_fasta_extension", default="fa", type=str, help = "File extension. Don't include period/fullstop/. [Default: fa]")
    # parser_io.add_argument("--protein_fasta_extension", default="faa", type=str, help = "File extension. Include the period/fullstop/. [Default: faa]")
    parser_io.add_argument("-a", "--absolute", action="store_true", help = "Use absolute paths instead of relative paths")
    parser_io.add_argument("-e", "--allow_missing_files", action="store_true", help = "Allow missing files")
    parser_io.add_argument("--volume_prefix", type=str, help = "Docker container prefix to volume path")
    parser_io.add_argument("--header", action="store_true", help = "Write header")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    output = defaultdict(dict)
    # Build table from binning directory
    if opts.binning_directory:
        for fp in tqdm(glob.glob(os.path.join(opts.binning_directory, "*", "*", "output", "genomes", "*")), "Reading files in {}".format(opts.binning_directory)):
            organism_type = fp.split("/")[-5]
            id_sample = fp.split("/")[-4]
            id_mag = ".".join(fp.split("/")[-1].split(".")[:-1])

            if fp.endswith(".fa"):
                output[(organism_type, id_sample, id_mag)]["genome"] = fp
            if fp.endswith(".faa"):
                output[(organism_type, id_sample, id_mag)]["proteins"] = fp
            if fp.endswith(".ffn"):
                output[(organism_type, id_sample, id_mag)]["cds"] = fp
            if fp.endswith(".gff"):
                output[(organism_type, id_sample, id_mag)]["gene_models"] = fp

    df_output = pd.DataFrame(output).reindex([ "genome", "proteins", "cds", "gene_models"]).T
    assert not df_output.empty, "Did not find any matches in the following directory: {}".format(opts.binning_directory)
    df_output.index.names = ["organism_type",  "id_sample", "id_mag"]

    # Check missing values
    if not opts.allow_missing_files:
        assert not np.any(df_output.isnull()), "Missing files detected.  Check error or allow with --allow_missing_files: \n{}".format(df_output.loc[df_output.isnull().sum(axis=1) > 0].to_string())
    
    # Absolute paths
    if opts.absolute:
        df_output = df_output.applymap(lambda fp: os.path.abspath(fp))

    # Docker volume prefix
    if opts.volume_prefix:
        df_output = df_output.applymap(lambda fp: os.path.join(opts.volume_prefix, fp) if pd.notnull(fp) else fp)

    if opts.output == "stdout":
        opts.output = sys.stdout 

    df_output.to_csv(opts.output, sep="\t", header=bool(opts.header))

if __name__ == "__main__":
    main()
