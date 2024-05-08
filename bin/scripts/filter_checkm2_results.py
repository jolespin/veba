#!/usr/bin/env python
import sys, os, glob, argparse
from shutil import copyfile
from typing import OrderedDict
# from collections import OrderedDict
import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.1.25"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <quality_report.tsv> -b <bin_directory> -o <output_directory> -f <scaffolds.fasta> -m 1500".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-i","--checkm2_results", type=str, required=True, help = "path/to/checkm_output/quality_report.tsv")
    parser.add_argument("-b","--bin_directory", type=str, help = "path/to/mag/directory")
    parser.add_argument("-o","--output_directory", type=str, default="checkm_filtered_output", help = "path/to/output_directory [Default: checkm_filtered_output]")
    parser.add_argument("-x","--extension", type=str, default="fa", help = "Fasta file extension for bins [Default: fa]")
    parser.add_argument("-f", "--fasta", type=str, help = "path/to/scaffolds.fasta. Include this only if you want to list of unbinned contigs")
    parser.add_argument("-m", "--minimum_contig_length", type=int, default=1, help="Minimum contig length. [Default: 1]")
    parser.add_argument("-u", "--unbinned", action="store_true", help="Write unbinned fasta sequences to file")
    parser.add_argument("--completeness", type=float, default=50.0, help = "CheckM completeness [Default: 50.0]")
    parser.add_argument("--contamination", type=float, default=10.0, help = "CheckM contamination [Default: 10.0]")
    parser.add_argument("--symlink", action="store_true", help = "Symlink MAGs instead of copying")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Load CheckM
    df_quality_report = pd.read_csv(opts.checkm2_results, sep="\t", index_col=0)

    # Quality assessment on MAGs
    mags = list()

    for id_mag, row in tqdm(df_quality_report.iterrows(), "Filtering MAGs", unit=" MAG"):
        completeness = row["Completeness"]
        contamination = row["Contamination"]

        # Conditions
        conditions = [
            completeness >= opts.completeness,
            contamination < opts.contamination,
        ]
        if all(conditions):
            mags.append(id_mag)

    # Output filtered 
    os.makedirs(opts.output_directory, exist_ok=True)

    df_quality_report.loc[mags,:].to_csv(os.path.join(opts.output_directory,"checkm2_results.filtered.tsv"), sep="\t") # Change to checkm_results to stay consistent with busco_results or other way aorund

    if opts.bin_directory:
        # Output MAGs
        os.makedirs(os.path.join(opts.output_directory,"genomes"), exist_ok=True)

        #bins.list
        with open(os.path.join(opts.output_directory, "bins.list"), "w") as f_bins:
            for id_mag in sorted(mags):
                print(id_mag, file=f_bins)

        # binned.list
        f_binned_list = open(os.path.join(opts.output_directory, "binned.list"), "w")

        binned_contigs = list() 
        scaffold_to_mag = OrderedDict()

        for id_mag in tqdm(mags, "Copying fasta files and writing binned contigs", unit=" MAG"):

            for src in glob.glob(os.path.join(opts.bin_directory, "{}.*".format(id_mag))):
                fn = os.path.split(src)[1]
                dst = os.path.join(opts.output_directory,"genomes", fn)
                # Copy or symlink MAGs                
                if opts.symlink:
                    if os.path.exists(dst):
                        os.remove(dst)
                    os.symlink(os.path.realpath(src), dst)
                else:
                    copyfile(src,dst)

            # Write contigs to list
            path_genome_assembly = os.path.join(opts.bin_directory, "{}.{}".format(id_mag, opts.extension))
            with open(path_genome_assembly, "r") as f_fasta:
                for line in f_fasta.readlines():
                    line = line.strip()
                    if line.startswith(">"):
                        header = line[1:]
                        id_contig = header.split(" ")[0] # Remove fasta description
                        print(id_contig, file=f_binned_list)
                        binned_contigs.append(id_contig)
                        scaffold_to_mag[id_contig] = id_mag
        f_binned_list.close()
        scaffold_to_mag = pd.Series(scaffold_to_mag)
        scaffold_to_mag.to_frame().to_csv(os.path.join(opts.output_directory,"scaffolds_to_bins.tsv"), sep="\t", header=None)

        # Get unbinned contigs
        if opts.fasta:

            from Bio.SeqIO.FastaIO import SimpleFastaParser

            f_unbinned_list = open(os.path.join(opts.output_directory, "unbinned.list"), "w")

            with open(opts.fasta, "r") as f_fasta: # Use stdin?
                for header, seq in tqdm(SimpleFastaParser(f_fasta), "Extracting unbinned contigs", unit=" contig"):
                    id_contig = header.split(" ")[0]
                    conditions = [
                        id_contig not in binned_contigs,
                        len(seq) >= opts.minimum_contig_length,
                    ]
                    if all(conditions):
                        print(id_contig, file=f_unbinned_list)
    
            f_unbinned_list.close()

            if opts.unbinned:
                f_unbinned_fasta = open(os.path.join(opts.output_directory, "unbinned.fasta"), "w")
                with open(opts.fasta, "r") as f_fasta: # Use stdin?
                    for header, seq in tqdm(SimpleFastaParser(f_fasta), "Writing unbinned contigs", unit=" contig"):
                        id_contig = header.split(" ")[0]
                        conditions = [
                            id_contig not in binned_contigs,
                            len(seq) >= opts.minimum_contig_length,
                        ]
                        if all(conditions):
                            print(">{}\n{}".format(id_contig, seq), file=f_unbinned_fasta)
                f_unbinned_fasta.close()

    

if __name__ == "__main__":
    main()
    
                