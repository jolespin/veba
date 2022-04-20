#!/usr/bin/env python
import sys, os, glob, argparse 
from collections import defaultdict
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.04.08"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i  <scaffolds_to_samples.tsv> -b <binning_directory> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-i","--scaffolds_to_samples", type=str, required=True,  help = "path/to/scaffolds_to_bins.tsv, Format: [id_scaffold]<tab>[id_sample], No header")
    parser.add_argument("-b","--binning_directory", type=str, required=True, help = "path/to/binning_directory")
    parser.add_argument("-o","--output_directory", type=str, help = "path/to/output_directory. Default is binning_directory")
    parser.add_argument("-x","--extension", default="fa", type=str, help = "Binning file extension [Default: fa]")
    # parser.add_argument("--bin_prefix", type=str,  default="%s__%s", help="Bin prefix for each partitioned bin. First %s refers to sample and second %s refers to original bin name [Default: %s__%s which is [id_sample]__[id_bin]")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    if not opts.output_directory:
        opts.output_directory = opts.binning_directory

    os.makedirs(opts.output_directory, exist_ok=True)

    #Scaffolds to samples
    scaffolds_to_samples = pd.read_csv(opts.scaffolds_to_samples, sep="\t", index_col=0, header=None).iloc[:,0]

    # Create sample subdirectories
    samples = sorted(scaffolds_to_samples.unique())
    for id_sample in samples:
        os.makedirs(os.path.join(opts.output_directory, id_sample), exist_ok=True)

    # Files
    bin_to_file = dict()
    for fp in glob.glob(os.path.join(opts.binning_directory, "*.{}".format(opts.extension))):
        id_coassembly = os.path.split(fp)[-1][:-(len(opts.extension)+1)]
        
        # Get all the samples that have said bin
        samples_in_bin = set()
        with open(fp, "r") as f_in:
            for id_scaffold, seq in tqdm(SimpleFastaParser(f_in), desc="Identifying samples that contributed to bin {}".format(id_coassembly)):
                id_scaffold = id_scaffold.split(" ")[0]
                id_sample = scaffolds_to_samples[id_scaffold]
                samples_in_bin.add(id_sample)
        for id_sample in samples_in_bin:
            id_bin = "{}__{}".format(id_sample, id_coassembly)
            bin_to_file[id_bin] = open(os.path.join(opts.output_directory, id_sample, "{}.{}".format(id_bin, opts.extension)), "w")

        with open(fp, "r") as f_in:
            for id_scaffold, seq in tqdm(SimpleFastaParser(f_in), "Partitioning {} sample-specific bins from {}".format(len(samples_in_bin), id_coassembly)):
                id_scaffold = id_scaffold.split(" ")[0]
                id_sample = scaffolds_to_samples[id_scaffold]
                id_bin = "{}__{}".format(id_sample, id_coassembly)
                print(">{}\n{}".format(id_scaffold, seq), file=bin_to_file[id_bin])
        for f in bin_to_file.values():
            f.close()

        
if __name__ == "__main__":
    main()
    
                
