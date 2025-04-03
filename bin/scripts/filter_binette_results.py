#!/usr/bin/env python
import sys, os, glob, argparse, re
from shutil import copyfile
from typing import OrderedDict
# from collections import OrderedDict
import pandas as pd
from tqdm import tqdm
import xxhash

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2025.4.2"

def parse_binette_initial_bin_combinations(combination:str):
    """
    Parses a list of bin combinations to extract unique bins.

    Parameters:
        bin_combinations (list of str): A list of strings containing bin combinations.

    Returns:
        set: A set of unique bins as strings.
    """
    unique_bins = set()

    # Split by delimiters (-, &, |)
    parts = re.split(r"[-&|]", combination)

    for part in parts:
        # Strip whitespace and check if the part is numeric or a bin identifier
        part = part.strip()
        if part:
            unique_bins.add(part)
    return unique_bins

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <binette_directory> -o <output_directory> -f <scaffolds.fasta> -m 1500".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-i","--binette_directory", type=str, required=True, help = "path/to/binette_directory/")
    parser.add_argument("-o","--output_directory", type=str, help = "path/to/output_directory/ [Default: path/to/binette_directory/filtered/]")
    parser.add_argument("-p","--bin_prefix", type=str, default="BINETTE__", help = "Bin prefix [Default: 'BINETTE__']")
    parser.add_argument("-x","--extension", type=str, default="fa", help = "Fasta file extension for bins [Default: fa]")
    parser.add_argument("-f", "--fasta", type=str, help = "path/to/scaffolds.fasta. Include this only if you want to list of unbinned contigs")
    parser.add_argument("-m", "--minimum_contig_length", type=int, default=1, help="Minimum contig length. [Default: 1]")
    parser.add_argument("--completeness", type=float, default=50.0, help = "CheckM2 completeness [Default: 50.0]")
    parser.add_argument("--contamination", type=float, default=10.0, help = "CheckM2 contamination [Default: 10.0]")
    parser.add_argument("-u", "--unbinned", action="store_true", help="Write unbinned fasta sequences to file")
    parser.add_argument("-e", "--exclude",  help="List of genomes to exclude (e.g., eukaryotic genomes)")
    parser.add_argument("-d","--domain_predictions", type=str,  help = "Tab-seperated table of domain predictions [id_genome]<tab>[id_domain], No header.")

    # parser.add_argument("--symlink", action="store_true", help = "Symlink MAGs instead of copying")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    # Exclusion
    exclude_mags = set()
    if opts.exclude:
        with open(opts.exclude, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    exclude_mags.add(line)
    if exclude_mags:
        print("Excluding the following MAGs:", exclude_mags, file=sys.stderr)
        
    # Use hash
    magold_to_hash = dict()
    for filepath in tqdm(glob.glob(os.path.join(opts.binette_directory, "final_bins", "*.fa")), "Creating hashes for each MAG", unit=" MAGs"):
        id_mag = filepath.split("/")[-1][:-3]
        contigs = set()
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    header = line[1:]
                    id_contig = header.split(" ", maxsplit=1)[0]
                    contigs.add(id_contig)
        contigs_repr = repr(sorted(contigs))
        magold_to_hash[id_mag] = xxhash.xxh32(contigs_repr, seed=0).hexdigest()
                
    # Load CheckM2
    # bin_id  origin  name    completeness    contamination   score   size    N50     contig_count
    df_quality_report = pd.read_csv(os.path.join(opts.binette_directory, "final_bins_quality_reports.tsv"),sep="\t", index_col=0)
    df_quality_report.index = df_quality_report.index.map(lambda x: f"bin_{x}")
    
    if opts.domain_predictions:
        print("Adding domain predictions", file=sys.stderr)
        mag_to_domain = pd.read_csv(opts.domain_predictions, sep="\t", index_col=0, header=None).iloc[:,0]
        df_quality_report["domain"] = mag_to_domain
    
    # Quality assessment on MAGs
    magold_to_magnew = dict()
    mags = list()

    for id_mag, row in tqdm(df_quality_report.iterrows(), "Filtering MAGs", unit=" MAGs"):
        origin = row["origin"]
        name = row["name"]
        completeness = row["completeness"]
        contamination = row["contamination"]

        # Conditions
        conditions = [
            completeness >= opts.completeness,
            contamination < opts.contamination,
            id_mag not in exclude_mags,
        ]
        if all(conditions):
            if (origin not in {"diff", "union", "intersec"}) and  (";" not in origin):
                new_mag = name
            else:
                new_mag = f"{opts.bin_prefix}{magold_to_hash[id_mag]}"
                # new_mag = f"{opts.bin_prefix}{id_mag}"
            magold_to_magnew[id_mag] = new_mag
            mags.append(id_mag)

    # Output filtered 
    if not opts.output_directory:
        opts.output_directory = os.path.join(opts.binette_directory, "filtered")
    os.makedirs(opts.output_directory, exist_ok=True)
    os.makedirs(os.path.join(opts.output_directory,"genomes"), exist_ok=True)
    
    # Initial bins
    dataframes = list()
    for filepath_initial in glob.glob(os.path.join(opts.binette_directory, "input_bins_quality_reports", "*.tsv")):
        df = pd.read_csv(filepath_initial, sep="\t", index_col=0)
        dataframes.append(df)
    df_initial = pd.concat(dataframes, axis=0)
    df_initial.index = df_initial.index.map(str)
    initial_to_genome = df_initial["name"].to_dict()

    # Build quality report
    df_quality_report_filtered = df_quality_report.loc[mags,:]
    df_quality_report_filtered.index = df_quality_report_filtered.index.map(lambda x: magold_to_magnew[x])
    
    # Add initial bins to quality report
    genome_to_initial = dict()
    for id_genome, combination in df_quality_report_filtered["name"].items():
        initial_bins = parse_binette_initial_bin_combinations(combination)
        initial_bins_labeled = set()
        for id_bin_number in initial_bins:
            try:
                id_initial_bin = initial_to_genome[id_bin_number]
                initial_bins_labeled.add((id_bin_number, id_initial_bin))
            except KeyError:
                initial_bins_labeled.add(id_bin_number)
        genome_to_initial[id_genome] = initial_bins_labeled
    df_quality_report_filtered["initial_bins"] = pd.Series(genome_to_initial, dtype=object)
    
    # Write quality report
    # df_quality_report_filtered.columns = df_quality_report_filtered.columns.map(str.capitalize)
    df_quality_report_filtered = df_quality_report_filtered.sort_index()
    df_quality_report_filtered.indexname = "id_genome"
    df_quality_report_filtered.to_csv(os.path.join(opts.output_directory,"checkm2_results.filtered.tsv"), sep="\t") 

    #bins.list
    with open(os.path.join(opts.output_directory, "bins.list"), "w") as f_bins:
        for id_mag in sorted(mags):
            print(magold_to_magnew[id_mag], file=f_bins)

    # binned.list
    f_binned_list = open(os.path.join(opts.output_directory, "binned.list"), "w")

    binned_contigs = set() 
    scaffold_to_mag = OrderedDict()
    for id_mag in tqdm(mags, "Copying fasta files and writing binned contigs", unit=" MAGs"):

        for src in glob.glob(os.path.join(opts.binette_directory, "final_bins", "{}.*".format(id_mag))):
            dst = os.path.join(opts.output_directory,"genomes", "{}.fa".format(magold_to_magnew[id_mag]))
            copyfile(src,dst)

        # Write contigs to list
        path_genome_assembly = os.path.join(opts.binette_directory, "final_bins", "{}.{}".format(id_mag, opts.extension))
        with open(path_genome_assembly, "r") as f_fasta:
            for line in f_fasta:
                line = line.strip()
                if line.startswith(">"):
                    header = line[1:]
                    id_contig = header.split(" ")[0].strip() # Remove fasta description
                    print(id_contig, file=f_binned_list)
                    binned_contigs.add(id_contig)
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
            for id_contig, seq in tqdm(pyfastx.Fasta(opts.fasta, build_index=False), "Writing unbinned contigs", unit=" contig"):
                
                conditions = [
                    id_contig not in binned_contigs,
                    len(seq) >= opts.minimum_contig_length,
                ]
               
                if all(conditions):
                    print(id_contig, file=f_unbinned_list)
                    print(">{}\n{}".format(id_contig, seq), file=f_unbinned_fasta)
            f_unbinned_fasta.close()
        else:
            for id_contig, seq in tqdm(pyfastx.Fasta(opts.fasta, build_index=False), "Writing unbinned contigs", unit=" contig"):
                conditions = [
                    id_contig not in binned_contigs,
                    len(seq) >= opts.minimum_contig_length,
                ]
                if all(conditions):
                    print(id_contig, file=f_unbinned_list)
        f_unbinned_list.close()

if __name__ == "__main__":
    main()
