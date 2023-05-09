#!/usr/bin/env python
import sys, os, glob, argparse, gzip
from collections import OrderedDict
import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.4.25"

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
    parser.add_argument("-i","--binning_directory", type=str, help = "path/to/binning_directory [Cannot use with --genomes]")
    parser.add_argument("-g","--genomes", type=str, default="stdin", help = "path/to/genomes as either .list with each line as a path to genome.fasta or a table [id_genome]<tab>[path/to/fasta] [Cannot use with --binning_directory] [Default: stdin]")
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
    

    if opts.binning_directory:
        if opts.genomes == "stdin":
            assert sys.stdin.isatty(), "--genomes cannot have any stdin if --binning_directory is selected"
        else:
            assert not opts.genomes,  "--genomes cannot be used if --binning_directory is selected"
    else:
        if opts.genomes != "stdin":
            assert os.path.exists(opts.genomes)
        else:
            assert not sys.stdin.isatty(), "If --binning_directory is not selected and --genomes is not provided explicit file then stdin is expected"

    # assert bool(opts.binning_directory) != bool(opts.genomes), "Must choose either --binning_directory or --genomes, not both."


    if not opts.bin_prefix:
        opts.bin_prefix = ""

    scaffold_to_bin = OrderedDict()

    if opts.binning_directory:
        assert os.path.exists(opts.binning_directory), "{} does not exist".format(opts.binning_directory)
        for filepath in glob.glob(os.path.join(opts.binning_directory, "*.{}".format(opts.extension))):
            id_bin = filepath.split("/")[-1][:-1*(len(opts.extension)+1)]
            id_bin = "{}{}".format(opts.bin_prefix, id_bin)

            if opts.extension.endswith(".gz"):
                f = gzip.open(filepath, "rt")
            else:
                f = open(filepath, "r")
            for line in tqdm(f, filepath):
                if opts.header:
                    if line.startswith(">"):
                        id_scaffold = line.strip()[1:]
                        scaffold_to_bin[id_scaffold] = id_bin 
                else:
                    if line.startswith(">"):
                        id_scaffold = line.strip()[1:].split(" ")[0]
                        scaffold_to_bin[id_scaffold] = id_bin 
            f.close()

    else:
        if opts.genomes:
            if opts.genomes == "stdin":
                opts.genomes = sys.stdin 

            # Load table
            df_genomes = pd.read_csv(opts.genomes, sep="\t", index_col=0, header=None)

            # Convert list to table
            if df_genomes.shape[1] == 0:
                bin_to_filepath = dict()
                for line in df_genomes.index:
                    filepath = line.strip()
                    id_bin = filepath.split("/")[-1][:-1*(len(opts.extension)+1)]
                    id_bin = "{}{}".format(opts.bin_prefix, id_bin)
                    bin_to_filepath[id_bin] = filepath
                df_genomes = pd.Series(bin_to_filepath).to_frame()

            # Read fasta files
            for id_bin, filepath in df_genomes.iloc[:,0].items():
                if opts.extension.endswith(".gz"):
                    f = gzip.open(filepath, "rt")
                else:
                    f = open(filepath, "r")
                for line in tqdm(f, filepath):
                    if opts.header:
                        if line.startswith(">"):
                            id_scaffold = line.strip()[1:]
                            scaffold_to_bin[id_scaffold] = id_bin 
                    else:
                        if line.startswith(">"):
                            id_scaffold = line.strip()[1:].split(" ")[0]
                            scaffold_to_bin[id_scaffold] = id_bin 
                f.close()
        
    df_output = pd.Series(scaffold_to_bin).to_frame(opts.bin_column_name)
    df_output.index.name = opts.scaffold_column_name 
    df_output = df_output.reset_index()
    if opts.column_order == "bin,scaffold":
        df_output = df_output.iloc[:,[1,0]]
    df_output.to_csv(sys.stdout, sep=opts.sep, header=opts.header, index=None)

if __name__ == "__main__":
    main()
    
                

