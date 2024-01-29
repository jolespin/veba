#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
from collections import defaultdict
import pandas as pd

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.7.18"



#
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
    parser.add_argument("-i","--preprocess_directory", required=True, type=str, help = "path/to/preprocess directory")
    parser.add_argument("-b","--basename", default="cleaned", type=str, help = "File basename to search for in [preprocess_directory]/[id_sample]/[output]/[basename]_1/2.fastq.gz [Default: cleaned]")
    parser.add_argument("-o","--output", default="stdout", type=str, help = "Output filepath [Default: stdout]")
    parser.add_argument("-a", "--absolute", action="store_true", help = "Use absolute paths instead of relative paths")
    parser.add_argument("--header", action="store_true", help = "Write header")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Build table
    output = defaultdict(dict)
    for fp in glob.glob(os.path.join(opts.preprocess_directory, "*", "output", "{}_1.fastq.gz".format(opts.basename))):
        id_sample = fp.split("/")[-3]
        output[id_sample]["Forward_Reads"] = fp
    for fp in glob.glob(os.path.join(opts.preprocess_directory, "*", "output", "{}_2.fastq.gz".format(opts.basename))):
        id_sample = fp.split("/")[-3]
        output[id_sample]["Reverse_Reads"] = fp
    df_output = pd.DataFrame(output).T 
    if opts.absolute:
        df_output = df_output.applymap(lambda fp: os.path.abspath(fp))
    df_output.index.name = "SampleID"

    if opts.output == "stdout":
        opts.output = sys.stdout 
    df_output.to_csv(opts.output, sep="\t", header=bool(opts.header))

if __name__ == "__main__":
    main()
