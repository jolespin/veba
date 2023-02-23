#!/usr/bin/env python
import sys, os, glob, argparse
from collections import defaultdict
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
    usage = "{} -i  <scaffolds_to_samples.tsv> -b <binning_directory> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-i","--scaffolds_to_samples", type=str, required=True,  help = "path/to/scaffolds_to_bins.tsv, Format: [id_scaffold]<tab>[id_sample], No header")
    parser.add_argument("-b","--binning_directory", type=str, required=True, help = "path/to/binning_directory")
    parser.add_argument("-o","--output_directory", type=str, help = "path/to/output_directory. Default is [binning_directory]/multisplit")
    parser.add_argument("-x","--extension", default="fa", type=str, help = "Binning file extension [Default: fa]")
    parser.add_argument("-d", "--delimiter", type=str,  default="__", help="Delimiter between [id_sample]<delimiter>[id_bin] [Default: __ which is [id_sample]__[id_bin]")
    # parser.add_argument("-r", "--remove_composite", action="store_true", help="Remove the original composite bin files")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    if not opts.output_directory:
        opts.output_directory = os.path.join(opts.binning_directory, "multisplit")

    os.makedirs(opts.output_directory, exist_ok=True)

    #Scaffolds to samples
    scaffold_to_sample = pd.read_csv(opts.scaffolds_to_samples, sep="\t", index_col=0, header=None).iloc[:,0]
    sample_to_scaffolds = defaultdict(set)
    for id_scaffold, id_sample in scaffold_to_sample.items():
        sample_to_scaffolds[id_sample].add(id_scaffold)

    # Create sample subdirectories
    samples = sorted(scaffold_to_sample.unique())
    for id_sample in samples:
        os.makedirs(os.path.join(opts.output_directory, id_sample), exist_ok=True)

    # scaffolds_to_bins.tsv
    query_filepath = os.path.join(opts.binning_directory, "scaffolds_to_bins.tsv")
    assert os.path.exists(query_filepath), "--binning_directory must contain scaffolds_to_bins.tsv"
    sample_to_bins = dict()
    for id_sample in samples:
        scaffold_to_bin = pd.read_csv(query_filepath, sep="\t", index_col=0, header=None).iloc[:,0]
        scaffold_to_bin = scaffold_to_bin.loc[sorted(set(scaffold_to_bin.index) & set(sample_to_scaffolds[id_sample]))]
        scaffold_to_bin = scaffold_to_bin.map(lambda x: "{}{}{}".format(id_sample, opts.delimiter, x))
        scaffold_to_bin.to_frame().to_csv(os.path.join(opts.output_directory, id_sample, "scaffolds_to_bins.tsv"), sep="\t", header=None)
        sample_to_bins[id_sample] = set(scaffold_to_bin.unique())

        with open(os.path.join(opts.output_directory, id_sample, "bins.list"), "w") as f_out:
            print(*sorted(scaffold_to_bin.unique()), sep="\n", file=f_out)

        with open(os.path.join(opts.output_directory, id_sample, "unbinned.list"), "w") as f_out:
            unbinned = set(sample_to_scaffolds[id_sample]) - set(scaffold_to_bin.index)
            print(*sorted(unbinned), sep="\n", file=f_out)

    # if opts.remove_composite:
    #     os.remove(query_filepath)

    # binned.list
    query_filepath = os.path.join(opts.binning_directory, "binned.list")
    if os.path.exists(query_filepath):
        # Create sample-specific files
        sample_to_file = dict() 
        for id_sample in samples:
            sample_to_file[id_sample] = open(os.path.join(opts.output_directory, id_sample, "binned.list"), "w")
        # Partition content
        with open(query_filepath, "r") as f_in:
            for line in f_in.readlines():
                id_scaffold = line.strip()
                id_sample = scaffold_to_sample[id_scaffold]
                print(id_scaffold, file=sample_to_file[id_sample])
        # Close sample-specific files
        for id_sample, f_out in sample_to_file.items():
            f_out.close()
    # if opts.remove_composite:
    #     os.remove(query_filepath)

    # unbinned.fasta
    query_filepath = os.path.join(opts.binning_directory, "unbinned.fasta")
    if os.path.exists(query_filepath):
        # Create sample-specific files
        sample_to_file = dict() 
        for id_sample in samples:
            sample_to_file[id_sample] = open(os.path.join(opts.output_directory, id_sample, "unbinned.fasta"), "w")
        # Partition content
        with open(query_filepath, "r") as f_in:
            for id_scaffold, seq in tqdm(SimpleFastaParser(f_in), desc="Identifying unbinned contigs in sample {}".format(id_sample)):
                id_sample = scaffold_to_sample[id_scaffold]
                print(">{}\n{}".format(id_scaffold, seq), file=sample_to_file[id_sample])
        # Close sample-specific files
        for id_sample, f_out in sample_to_file.items():
            f_out.close()

        os.remove(query_filepath)

    # # bins.list
    # query_filepath = os.path.join(opts.binning_directory, "bins.list")
    # if opts.remove_composite:
    #     os.remove(query_filepath)

    # # # genome_statistics.tsv
    # query_filepath = os.path.join(opts.binning_directory, "genome_statistics.tsv")
    # if os.path.exists(query_filepath):
    #     if opts.remove_composite:
    #         print("Please rerun `seqkit stats -a -b -T {}`".format(os.path.join(opts.binning_directory, "*", "*.{}".format(opts.extension))), file=sys.stderr)
    #         os.remove(query_filepath)

    # unbinned.list
    query_filepath = os.path.join(opts.binning_directory, "unbinned.list")
    if os.path.exists(query_filepath):
        os.remove(query_filepath)

    # bins 
    query_filepath = os.path.join(opts.binning_directory, "bins")
    if os.path.exists(query_filepath):
        for id_sample in samples:
            os.makedirs(os.path.join(opts.output_directory, id_sample, "bins"), exist_ok=True)
        # Files
        bin_to_file = dict()
        for fp in glob.glob(os.path.join(opts.binning_directory, "bins", "*.{}".format(opts.extension))):
            id_bin_composite = os.path.split(fp)[-1][:-(len(opts.extension)+1)]
            
            # Get all the samples that have said bin
            samples_in_bin = set()
            with open(fp, "r") as f_in:
                for id_scaffold, seq in tqdm(SimpleFastaParser(f_in), desc="Identifying samples that contributed to bin {}".format(id_bin_composite)):
                    id_scaffold = id_scaffold.split(" ")[0]
                    id_sample = scaffold_to_sample[id_scaffold]
                    samples_in_bin.add(id_sample)
            for id_sample in samples_in_bin:
                id_bin = "{}{}{}".format(id_sample, opts.delimiter, id_bin_composite)
                bin_to_file[id_bin] = open(os.path.join(opts.output_directory, id_sample, "bins", "{}.{}".format(id_bin, opts.extension)), "w")

            with open(fp, "r") as f_in:
                for id_scaffold, seq in tqdm(SimpleFastaParser(f_in), "Partitioning {} sample-specific bins from {}".format(len(samples_in_bin), id_bin_composite)):
                    id_scaffold = id_scaffold.split(" ")[0]
                    id_sample = scaffold_to_sample[id_scaffold]
                    id_bin = "{}{}{}".format(id_sample, opts.delimiter, id_bin_composite)
                    print(">{}\n{}".format(id_scaffold, seq), file=bin_to_file[id_bin])

            for id_bin, f_out in bin_to_file.items():
                f_out.close()

        # if opts.remove_composite:
        #     shutil.rmtree(os.path.join(opts.output_directory, "bins"), ignore_errors=True)


        
if __name__ == "__main__":
    main()
    
                
