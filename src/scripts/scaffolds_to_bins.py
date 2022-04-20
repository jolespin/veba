#!/usr/bin/env python
import sys, os, glob, argparse, gzip
from collections import OrderedDict
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.03.26"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <binning_directory> -x <extension> --sep '\t' > <output.tsv>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--binning_directory", type=str, help = "path/to/binning_directory")
    parser.add_argument("-x","--extension", default="fa", type=str, help = "Binning file extension [Default: fa]")
    parser.add_argument("--sep", type=str, default="\t",  help = "Seperator [Default: '\t'")
    parser.add_argument("--scaffold_column_name", type=str, default="Scaffold", help="Scaffold column name [Default: Scaffold]")
    parser.add_argument("--bin_column_name", type=str, default="Bin", help="Bin column name [Default: Bin]")
    parser.add_argument("--column_order", type=str, default="scaffold,bin", help="Column order.  Specify either 'scaffold,bin' or 'bin,scaffold' [Default:scaffold,bin]")
    parser.add_argument("--bin_prefix", type=str,  help="Bin prefix. Default is to not have a prefix.")
    parser.add_argument("--header", action="store_true", help="Specify if header should be in output")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Parse
    assert opts.column_order in {"scaffold,bin", "bin,scaffold"}, "Must choose either 'scaffold,bin' or 'bin,scaffold' for --column_order"

    assert os.path.exists(opts.binning_directory), "{} does not exist".format(opts.binning_directory)

    if not opts.bin_prefix:
        opts.bin_prefix = ""
        
    scaffold_to_bin = OrderedDict()
    for filepath in glob.glob(os.path.join(opts.binning_directory, "*.{}".format(opts.extension))):
        id_bin = filepath.split("/")[-1][:-1*(len(opts.extension)+1)]
        id_bin = "{}{}".format(opts.bin_prefix, id_bin)

        if opts.extension.endswith(".gz"):
            f = gzip.open(filepath, "rt")
        else:
            f = open(filepath, "r")
        for line in f:
            if opts.header:
                if line.startswith(">"):
                    id_scaffold = line.strip()[1:]
                    scaffold_to_bin[id_scaffold] = id_bin 
            else:
                if line.startswith(">"):
                    id_scaffold = line.strip()[1:].split(" ")[0]
                    scaffold_to_bin[id_scaffold] = id_bin 
        f.close()
        
    df = pd.Series(scaffold_to_bin).to_frame(opts.bin_column_name)
    df.index.name = opts.scaffold_column_name 
    df = df.reset_index()
    if opts.column_order == "bin,scaffold":
        df = df.iloc[:,[1,0]]
    df.to_csv(sys.stdout, sep=opts.sep, header=opts.header, index=None)

if __name__ == "__main__":
    main()
    
                

