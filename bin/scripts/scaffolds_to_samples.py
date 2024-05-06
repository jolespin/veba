#!/usr/bin/env python
import sys, os, glob, argparse, gzip
from collections import OrderedDict
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm 

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.2.6"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <assembly_directory> --sep '\t' > <output.tsv>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--assembly_directory", type=str, help = "path/to/assembly_directory")
    parser.add_argument("-m", "--minimum_contig_length", type=int, default=1, help="Minimum contig length.  [Default: 1] ")
    parser.add_argument("--sep", type=str, default="\t",  help = "Seperator [Default: '\t'")
    parser.add_argument("--scaffold_column_name", type=str, default="Scaffold", help="Scaffold column name [Default: Scaffold]")
    parser.add_argument("--sample_column_name", type=str, default="Sample", help="Bin column name [Default: Sample]")
    parser.add_argument("--column_order", type=str, default="scaffold,sample", help="Column order.  Specify either 'scaffold,sample' or 'sample,scaffold' [Default:scaffold,bin]")
    parser.add_argument("--header", action="store_true", help="Specify if header should be in output")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Parse
    assert opts.column_order in {"scaffold,sample", "sample,scaffold"}, "Must choose either 'scaffold,sample' or 'sample,scaffold' for --column_order"

    assert os.path.exists(opts.assembly_directory), "{} does not exist".format(opts.assembly_directory)

    scaffold_to_sample = OrderedDict()
    for filepath in glob.glob(os.path.join(opts.assembly_directory, "*", "output", "scaffolds.fasta")):
        id_sample = filepath.split("/")[-3]
        with open(filepath, "r") as f_in:
            for id_scaffold, seq in tqdm(SimpleFastaParser(f_in), desc="Getting scaffold identifers for each sample: {}".format(id_sample)):
                if len(seq) >= opts.minimum_contig_length:
                    if id_scaffold in scaffold_to_sample:
                        print("[Duplicate] {} from {}".format(id_scaffold, id_sample), file=sys.stdout)
                    scaffold_to_sample[id_scaffold] = id_sample

    df = pd.Series(scaffold_to_sample).to_frame(opts.sample_column_name)
    df.index.name = opts.scaffold_column_name 
    df = df.reset_index()
    if opts.column_order == "sample,scaffold":
        df = df.iloc[:,[1,0]]
    df.to_csv(sys.stdout, sep=opts.sep, header=opts.header, index=None)

if __name__ == "__main__":
    main()
    
                

