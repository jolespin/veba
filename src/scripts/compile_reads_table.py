#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, re
from collections import defaultdict
import pandas as pd

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.10.24"

def parse_basename(query: str, naming_scheme: str):
    """
    Adapted from the following source: 
    * @s-b
    * https://stackoverflow.com/questions/73296040/how-to-parse-out-two-fields-out-of-string-using-regex-in-python
    """
    def repl(match_object):
        inside_bracket = match_object.group(1)
        if inside_bracket == "DIRECTION":
            return r"(?P<DIRECTION>[12])"
        if inside_bracket == "ID":
            return r"(?P<ID>[-.\w]+)"
    pattern = re.sub(r"\[(.*?)\]", repl, naming_scheme)
    match = re.match(pattern, query)
    return match["ID"], match["DIRECTION"]


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

    parser_preprocess_directory = parser.add_argument_group('[Mode 1] Preprocess Directory arguments')
    parser_preprocess_directory.add_argument("-i","--preprocess_directory",  type=str, help = "path/to/preprocess directory (e.g., veba_output/preprocess) [Cannot be used with --fastq_directory]")
    parser_preprocess_directory.add_argument("-b","--basename", default="cleaned", type=str, help = "File basename to search VEBA preprocess directory [preprocess_directory]/[id_sample]/[output]/[basename]_1/2.fastq.gz [Default: cleaned]")

    parser_fastq_directory = parser.add_argument_group('[Mode 2] Fastq Directory arguments')
    parser_fastq_directory.add_argument("-f","--fastq_directory",  type=str, help = "path/to/fastq_directory [Cannot be used with --preprocess_directory]")
    parser_fastq_directory.add_argument("-n","--naming_scheme", default="[ID]_R[DIRECTION]_001.fastq.gz", type=str, help = "Naming scheme.  Use [ID] for identifier name and [DIRECTION] for 1 or 2. [Default: [ID]_R[DIRECTION]_001.fastq.gz]")
    parser_fastq_directory.add_argument("-x","--extension", default=".fastq.gz", type=str, help = "File extension. Include the . [Default: .fastq.gz]")

    parser_output = parser.add_argument_group('Output arguments')
    parser_output.add_argument("-o","--output", default="stdout", type=str, help = "Output filepath [Default: stdout]")
    parser_output.add_argument("-r", "--relative", action="store_true", help = "Use relative paths instead of absolute paths")
    parser_output.add_argument("-0", "--sample_label", default="sample-id", type=str, help = "Sample ID column label [Reverse: sample-id]")
    parser_output.add_argument("-1", "--forward_label", default="forward-absolute-filepath", type=str, help = "Forward filepath column label [Default: forward-absolute-filepath]")
    parser_output.add_argument("-2", "--reverse_label", default="reverse-absolute-filepath", type=str, help = "Reverse filepath column label [Default: reverse-absolute-filepath]")
    parser_output.add_argument("--header", action="store_true", help = "Write header")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Checks
    assert bool(opts.preprocess_directory) != bool(opts.fastq_directory), "Please provide either --preprocess_directory or --fastq_directory (not both)"

    output = defaultdict(dict)
    # Build table from preprocess directory
    if opts.preprocess_directory:
        for fp in glob.glob(os.path.join(opts.preprocess_directory, "*", "output", "{}_1.fastq.gz".format(opts.basename))):
            id_sample = fp.split("/")[-3]
            output[id_sample][opts.forward_label] = fp
        for fp in glob.glob(os.path.join(opts.preprocess_directory, "*", "output", "{}_2.fastq.gz".format(opts.basename))):
            id_sample = fp.split("/")[-3]
            output[id_sample][opts.reverse_label] = fp
    # Build table from fastq directory
    if opts.fastq_directory:
        for fp in glob.glob(os.path.join(opts.fastq_directory, "*{}".format(opts.extension))):
            basename = fp.split("/")[-1]
            id_sample, direction = parse_basename(basename, naming_scheme=opts.naming_scheme)
            # id_sample = "_R".join(basename.split("_R")[:-1])
            if direction == "1":
                output[id_sample][opts.forward_label] = fp
            if direction == "2":
                output[id_sample][opts.reverse_label] = fp
    df_output = pd.DataFrame(output).T.sort_index().loc[:,[opts.forward_label, opts.reverse_label]]
    df_output.index.name = opts.sample_label
    
    # Check missing values
    missing_values = df_output.notnull().sum(axis=1)[lambda x: x < 2].index
    assert missing_values.size == 0, "Missing fastq for the following samples: {}".format(missing_values.index)
    
    # Absolute paths
    if not opts.relative:
        df_output = df_output.applymap(lambda fp: os.path.abspath(fp))
    else:
        if opts.header:
            if "absolute" in opts.forward_label.lower():
                print("You've selected --relative and may want to either not use a header or remove 'absolute' from the --forward_label: {}".format(opts.forward_label), file=sys.stderr)
            if "absolute" in opts.reverse_label.lower():
                print("You've selected --relative and may want to either not use a header or remove 'absolute' from the --reverse_label: {}".format(opts.reverse_label), file=sys.stderr)

    if opts.output == "stdout":
        opts.output = sys.stdout 
    df_output.to_csv(opts.output, sep="\t", header=bool(opts.header))

if __name__ == "__main__":
    main()
