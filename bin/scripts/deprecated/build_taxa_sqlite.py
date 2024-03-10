#!/usr/bin/env python
import sys, os, glob, argparse
from ete3 import NCBITaxa

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.04.18"


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <taxdump_file> -o <output_directory".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--taxdump_file", type=str, required=True, help = "path/to/taxdump.tar.gz (i.e., https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")
    parser.add_argument("-o","--output_directory", type=str, help = "Output directory: [Default: path/of/taxdump.tar.gz]")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # I/O
    if opts.output_directory is None:
        opts.output_directory = os.path.split(opts.taxdump_file)[0]
    os.makedirs(opts.output_directory, exist_ok=True)

    # Create
    NCBITaxa(taxdump_file=opts.taxdump_file, dbfile=os.path.join(opts.output_directory, "taxa.sqlite"))
    
if __name__ == "__main__":
    main()
    
                

