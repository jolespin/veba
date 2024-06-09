#!/usr/bin/env python
import sys, os, glob, argparse 
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.6.21"

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
    parser.add_argument("-g", "--minimum_genomes_aligned_ratio", type=float,  default=0.95, help = "Minimum ratio of genomes include in alignment. This removes markers that are under represented. [Default: 0.95]")
    parser.add_argument("-m", "--minimum_markers_aligned_ratio", type=float,  default=0.2, help = "Minimum ratio of markers aligned. This removes genomes with few markers. Note, this is based on detected markers and NOT total markers in original HMM. [Default: 0.2]")
    parser.add_argument("--prefiltered_alignment_table", type=str,  help = "Alignment table output of (n = genomes, m = markers) where each i,j is a MSA as a (gzipped) tab-separated file")
    parser.add_argument("--boolean_prefiltered_alignment_table", type=str,  help = "Alignment table output of (n = genomes, m = markers) where each i,j is in {0,1} as a (gzipped) tab-separated file")
    parser.add_argument("--alignment_table", type=str,  help = "Alignment table output of (n = genomes, m = markers) where each i,j is a MSA as a (gzipped) tab-separated file")
    parser.add_argument("--boolean_alignment_table", type=str,  help = "Alignment table output of (n = genomes, m = markers) where each i,j is in {0,1} as a (gzipped) tab-separated file")

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
    prefiltered_genomes = set()
    for fp in tqdm(glob.glob(os.path.join(opts.msa_directory, "*.{}".format(opts.ext))), "Reading MSA files from {}".format(opts.msa_directory), unit=" files"):
        id_marker = fp.split("/")[-1][:-(len(opts.ext)+1)]

        msa = dict() 
        with open(fp, "r") as f:
            for header, seq in SimpleFastaParser(f):
                id = header.split(" ")[0]
                msa[id] = seq 
        # msa = syu.read_fasta(fp, description=False)
        msa = pd.Series(msa)
        prefiltered_genomes |= set(msa.index)
        marker_to_msa[id_marker] = msa


    
    # number_of_genomes = len(genomes)
    # number_of_markers = len(marker_to_msa)
    # minimum_genomes_aligned = int(opts.minimum_genomes_aligned_ratio * number_of_genomes)
    # minimum_markers_aligned = int(opts.minimum_markers_aligned_ratio * number_of_markers)

    df_prefiltered_markers = pd.DataFrame(marker_to_msa)
    if opts.prefiltered_alignment_table:
        df_prefiltered_markers.to_csv(opts.prefiltered_alignment_table, sep="\t")
    if opts.boolean_prefiltered_alignment_table:
        df_prefiltered_markers.notnull().astype(int).to_csv(opts.boolean_prefiltered_alignment_table, sep="\t")

    df_bool = df_prefiltered_markers.notnull()
    index_genomes_filtered = df_bool.index[df_bool.mean(axis=1) >= opts.minimum_markers_aligned_ratio]
    print(
        "Removing the following genomes based on --minimum_markers_aligned_ratio {}:\n".format(opts.minimum_markers_aligned_ratio), 
        *sorted(prefiltered_genomes - set(index_genomes_filtered)),  
    sep="\n",
    file=sys.stderr,
    )

    df_bool = df_bool.loc[index_genomes_filtered]
    index_markers_filtered = df_bool.columns[df_bool.mean(axis=0) >= opts.minimum_genomes_aligned_ratio]
    print(
        "Removing the following markers based on --minimum_genomes_aligned_ratio {}:\n".format(opts.minimum_genomes_aligned_ratio), 
        *sorted(set(marker_to_msa.keys()) - set(index_markers_filtered)),  
        sep="\n",
        file=sys.stderr,
    )

    genomes = pd.Index(sorted(df_bool.index))
    markers = pd.Index(sorted(df_bool.columns))

    number_of_genomes = len(genomes)
    number_of_markers = len(markers)

    marker_to_msa = df_prefiltered_markers.loc[index_genomes_filtered, index_markers_filtered].to_dict()
        
    # # Remove low prevalence markers
    # for id_marker, msa in tqdm(list(marker_to_msa.items()), "Removing low prevalence markers", unit = " markers"):
    #     if len(msa) < minimum_genomes_aligned:
    #         print("Removing marker {} because it includes only {} genomes out of the total {} genomes included in the union of all alignments. Missing the following: {}".format(id_marker, len(msa), number_of_genomes, sorted(set(genomes) - set(msa.index))), file=sys.stderr)
    #         del marker_to_msa[id_marker]

    # Fill in missing genomes for markers
    for id_marker, msa in tqdm(marker_to_msa.items(), "Filling in missing genomes from alignment", unit = " markers"):
        msa = pd.Series(msa).dropna()
        lengths = msa.map(len)
        assert lengths.nunique() == 1, "MSA fasta alignments should all be the same length: {}".format(id_marker)
        n =  lengths.values[0]
        marker_to_msa[id_marker] = msa.reindex(genomes).fillna("-"*n)

    # Organize MSA by markers
    print("Merging MSA from {} markers".format(len(marker_to_msa)), file=sys.stderr)
    df_markers = pd.DataFrame(marker_to_msa)
    if opts.alignment_table:
        df_markers.to_csv(opts.alignment_table, sep="\t")
    if opts.boolean_alignment_table:
        df_markers.notnull().astype(int).to_csv(opts.boolean_alignment_table, sep="\t")
    # Merge MSA
    msa_concatenated = df_markers.apply(lambda x: "".join(x), axis=1)

    # Output
    print("Writing merged MSA", file=sys.stderr)
    f_out.writelines(">{}\n{}\n".format(id, seq) for id, seq in msa_concatenated.items())

    if f_out is not sys.stdout: 
        f_out.close()
    

if __name__ == "__main__":
    main()
    
                

