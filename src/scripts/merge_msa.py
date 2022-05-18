#!/usr/bin/env python
import sys, os, glob, argparse 
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.05.03"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <msa_directory> -x <ext> -o <output>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--msa_directory", type=str, help = "path/to/msa_directory")
    parser.add_argument("-x","--ext", type=str, default="msa.clipkit", help = "Complete file extension for MSA (anything before will be considered marker name) [Default: msa.clipkit")
    parser.add_argument("-o","--output", type=str,  default="stdout", help = "Output merged multiple sequence alignment [Default: stdout]")
    parser.add_argument("-m","--minimum_genomes_aligned_ratio", type=float,  default=0.5, help = "Minimum ratio of genomes include in alignment [Default: 0.5]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert 0 < opts.minimum_genomes_aligned_ratio <= 1.0, "--minimum_genomes_aligned_ratio must be in (0,1]"

    # Open output file
    if opts.output == "stdout":
        f_out = sys.stdout 
    else:
        f_out = open(opts.output, "w")

    # Get MSA for markers
    marker_to_msa = dict()
    samples = set()
    for fp in tqdm(glob.glob(os.path.join(opts.msa_directory, "*.{}".format(opts.ext))), "Reading MSA files from {}".format(opts.msa_directory), unit=" files"):
        id_marker = fp.split("/")[-1][:-(len(opts.ext)+1)]

        msa = dict() 
        with open(fp, "r") as f:
            for header, seq in SimpleFastaParser(f):
                id = header.split(" ")[0]
                msa[id] = seq 
        # msa = syu.read_fasta(fp, description=False)
        msa = pd.Series(msa)
        samples |= set(msa.index)
        marker_to_msa[id_marker] = msa
    samples = pd.Index(sorted(samples))
    number_of_samples = len(samples)
    minimum_genomes_aligned = int(opts.minimum_genomes_aligned_ratio * number_of_samples)

    # Remove low prevalence markers
    for id_marker, msa in tqdm(list(marker_to_msa.items()), "Removing low prevalence markers", unit = " markers"):
        if len(msa) < minimum_genomes_aligned:
            print("Removing marker {} because it includes only {} genomes out of the total {} genomes included in the union of all alignments".format(id_marker, minimum_genomes_aligned, number_of_samples), file=sys.stderr)
            del marker_to_msa[id_marker]

    # Fill in missing samples for markers
    for id_marker, msa in tqdm(marker_to_msa.items(), "Filling in missing samples from alignment", unit = " markers"):
        lengths = msa.map(len)
        assert lengths.nunique() == 1, "MSA fasta alignments should all be the same length: {}".format(id_marker)
        n =  lengths.values[0]
        marker_to_msa[id_marker] = msa.reindex(samples).fillna("-"*n)

    # Organize MSA by markers
    print("Merging MSA from {} markers".format(len(marker_to_msa)), file=sys.stderr)
    df_markers = pd.DataFrame(marker_to_msa)

    # Merge MSA
    msa_concatenated = df_markers.apply(lambda x: "".join(x), axis=1)

    # Output
    print("Writing merged MSA", file=sys.stderr)
    f_out.writelines(">{}\n{}\n".format(id, seq) for id, seq in msa_concatenated.items())

    if f_out is not sys.stdout: 
        f_out.close()
    

if __name__ == "__main__":
    main()
    
                

