#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, gzip, pickle

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.11.15"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__

    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <identifier_mapping.proteins.tsv[.gz]> -o <output.dict.pkl[.gz]>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", default="stdin", type=str, help = "Path to identifier mapping table [id_database]<tab>[id_source]<tab>[id_protein]<tab>[id_hash], No header. [Default: stdin]")
    parser.add_argument("-o","--output", required=True, type=str, help = "Path to dictionary pickle object.  Can be gzipped. (Recommended name: target_to_source.dict.pkl.gz)")
    parser.add_argument("-n","--number_of_sequences",  type=int, help = "Number of sequences.  If used, the tqdm is required.")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Input 
    f_in = None
    if opts.input == "stdin":
        f_in = sys.stdin 
    else:
        if opts.input.endswith(".gz"):
            f_in = gzip.open(opts.input, "rt")
        else:
            f_in = open(opts.input, "r")
    assert f_in is not None, "Unrecognized file format: {}".format(opts.input)

    if opts.number_of_sequences is not None:
        from tqdm import tqdm
        input_iterable = tqdm(f_in, total=opts.number_of_sequences, unit=" sequences")
    else:
        input_iterable = f_in

    print(" * Reading identifier mappings from the following file: {}".format(f_in), file=sys.stderr)
    target_to_source = dict()
    for line in input_iterable:
        line = line.strip()
        if line:
            fields = line.split("\t")
            id_hash = fields[3]
            id_source = fields[1]
            target_to_source[id_hash] = id_source
    if f_in != sys.stdin:
        f_in.close()

    print(" * Writing Python dictionary: {}".format(opts.output), file=sys.stderr)
    f_out = None 
    if opts.output.endswith((".gz", ".pgz")):
        f_out = gzip.open(opts.output, "wb")
    else:
        f_out = open(opts.output, "wb")
    assert f_out is not None, "Unrecognized file format: {}".format(opts.output)
    pickle.dump(target_to_source, f_out)


   




if __name__ == "__main__":
    main()
