#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse
from soothsayer_utils import format_header, read_script_as_module

script_directory  =  os.path.dirname(os.path.abspath( __file__ ))

try:
    from fastq_preprocessor import fastq_preprocessor_long
except ImportError:
    fastq_preprocessor_long = read_script_as_module("fastq_preprocessor_long", os.path.join(script_directory, "fastq_preprocessor_long.py"))

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.11.29"

if __name__ == "__main__":
    print(format_header("VEBA Preprocessing Wrapper (fastq_preprocessor v{})".format(fastq_preprocessor_long.__version__)), file=sys.stderr)
    label = "Mode: Long Nanopore and PacBio reads"
    print(label, file=sys.stderr)
    print(len(label)*"-", file=sys.stderr)
    fastq_preprocessor_long.main(sys.argv[1:])
