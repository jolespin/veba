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
__version__ = "2023.2.9"


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <veba_directory> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline

    parser.add_argument("-i","--veba_directory", type=str, help = "path/to/hmmsearch_tblout.tsv", required=True)
    parser.add_argument("-o","--output_directory", type=str,  default="veba_output/synopsis", help = "Output directory [Default: veba_output/synopsis]")
    parser.add_argument("-a", "--include_assembly", action="store_true", help="Output query identifiers only")
    parser.add_argument("-b", "--include_assembly", action="store_true", help="Output query identifiers only")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

 

if __name__ == "__main__":
    main()
