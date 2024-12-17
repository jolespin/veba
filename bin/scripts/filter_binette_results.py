#!/usr/bin/env python
import sys, os, glob, argparse
from shutil import copyfile
from typing import OrderedDict
# from collections import OrderedDict
import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.12.17"

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
    parser.add_argument("-i","--binette_directory", type=str, required=True, help = "path/to/binette_directory/")
    parser.add_argument("-o","--output_directory", type=str, help = "path/to/output_directory/ [Default: path/to/binette_directory/filtered/]")
    parser.add_argument("-p","--bin_prefix", type=str, default="BINETTE__", help = "Bin prefix [Default: 'BINETTE__']")
    parser.add_argument("-x","--extension", type=str, default="fa", help = "Fasta file extension for bins [Default: fa]")
    parser.add_argument("--completeness", type=float, default=50.0, help = "CheckM2 completeness [Default: 50.0]")
    parser.add_argument("--contamination", type=float, default=10.0, help = "CheckM2 contamination [Default: 10.0]")
    parser.add_argument("-f", "--fasta", type=str, help = "path/to/scaffolds.fasta. Include this only if you want to list of unbinned contigs")
    parser.add_argument("-m", "--minimum_contig_length", type=int, default=1, help="Minimum contig length. [Default: 1]")
    parser.add_argument("-u", "--unbinned", action="store_true", help="Write unbinned fasta sequences to file")
    # parser.add_argument("--symlink", action="store_true", help = "Symlink MAGs instead of copying")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Load CheckM2
    # bin_id  origin  name    completeness    contamination   score   size    N50     contig_count
    df_quality_report = pd.read_csv(os.path.join(opts.binette_directory, "final_bins_quality_reports.tsv"),sep="\t", index_col=0)

    # Quality assessment on MAGs
    magold_to_magnew = dict()
    mags = list()

    for id_mag, row in tqdm(df_quality_report.iterrows(), "Filtering MAGs", unit=" MAG"):
        completeness = row["completeness"]
        contamination = row["contamination"]

        # Conditions
        conditions = [
            completeness >= opts.completeness,
            contamination < opts.contamination,
        ]
        if all(conditions):
            magold_to_magnew[id_mag] = f"{opts.bin_prefix}{id_mag}"
            mags.append(id_mag)

    # Output filtered 
    if not opts.output_directory:
        opts.output_directory = os.path.join(opts.binette_directory, "filtered")
    os.makedirs(opts.output_directory, exist_ok=True)
    os.makedirs(os.path.join(opts.output_directory,"genomes"), exist_ok=True)

    df_quality_report_filtered = df_quality_report.loc[mags,:]
    df_quality_report_filtered.index = df_quality_report_filtered.index.map(lambda x: magold_to_magnew[x])
    df_quality_report_filtered.columns = df_quality_report_filtered.columns.map(str.capitalize)
    df_quality_report_filtered.indexname = "id_genome"
    df_quality_report_filtered.to_csv(os.path.join(opts.output_directory,"filtered","checkm2_results.filtered.tsv"), sep="\t") 


    #bins.list
    with open(os.path.join(opts.output_directory, "bins.list"), "w") as f_bins:
        for id_mag in sorted(mags):
            print(magold_to_magnew[id_mag], file=f_bins)

    # binned.list
    f_binned_list = open(os.path.join(opts.output_directory, "binned.list"), "w")

    binned_contigs = list() 
    scaffold_to_mag = OrderedDict()

    for id_mag in tqdm(mags, "Copying fasta files and writing binned contigs", unit=" MAG"):

        for src in glob.glob(os.path.join(opts.bin_directory, "{}.*".format(magold_to_magnew[id_mag]))):
            fn = os.path.split(src)[1]
            dst = os.path.join(opts.output_directory,"genomes", fn)
            copyfile(src,dst)

        # Write contigs to list
        path_genome_assembly = os.path.join(opts.bin_directory, "{}.{}".format(magold_to_magnew[id_mag], opts.extension))
        with open(path_genome_assembly, "r") as f_fasta:
            for line in f_fasta:
                line = line.strip()
                if line.startswith(">"):
                    header = line[1:]
                    id_contig = header.split(" ")[0] # Remove fasta description
                    print(id_contig, file=f_binned_list)
                    binned_contigs.append(id_contig)
                    scaffold_to_mag[id_contig] = magold_to_magnew[id_mag]
    f_binned_list.close()
    scaffold_to_mag = pd.Series(scaffold_to_mag)
    scaffold_to_mag.to_frame().to_csv(os.path.join(opts.output_directory,"scaffolds_to_bins.tsv"), sep="\t", header=None)

    # Get unbinned contigs
    if opts.fasta:
        import pyfastx

        f_unbinned_list = open(os.path.join(opts.output_directory, "unbinned.list"), "w")

        if opts.unbinned:
            f_unbinned_fasta = open(os.path.join(opts.output_directory, "unbinned.fasta"), "w")
            for id, seq in tqdm(pyfastx.Fasta(f_fasta, build_index=False), "Writing unbinned contigs", unit=" contig"):
                id_contig = header.split(" ")[0]
                conditions = [
                    id_contig not in binned_contigs,
                    len(seq) >= opts.minimum_contig_length,
                ]
                if all(conditions):
                    print(id_contig, file=f_unbinned_list)
                    if all(conditions):
                        print(">{}\n{}".format(id_contig, seq), file=f_unbinned_fasta)
                        
        else:
            for id, seq in tqdm(pyfastx.Fasta(f_fasta, build_index=False), "Writing unbinned contigs", unit=" contig"):
                id_contig = header.split(" ")[0]
                conditions = [
                    id_contig not in binned_contigs,
                    len(seq) >= opts.minimum_contig_length,
                ]
                if all(conditions):
                    print(id_contig, file=f_unbinned_list)
        f_unbinned_list.close()

if __name__ == "__main__":
    main()
    
                