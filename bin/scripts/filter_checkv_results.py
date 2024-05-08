#!/usr/bin/env python
import sys, os, glob, argparse
from shutil import copyfile
# from collections import OrderedDict
import pandas as pd
from tqdm import tqdm
from Bio.SeqIO.FastaIO import SimpleFastaParser

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.2.14"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <checkv_output.tsv> -f <scaffolds.fasta> -o <output_directory>  -m 1000".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-i","--checkv_results", type=str, required=True, help = "path/to/checkv/quality_summary.tsv")
    parser.add_argument("-f", "--fasta", required=True, type=str, help = "path/to/scaffolds.fasta")
    parser.add_argument("-c", "--completeness_tsv", required=True, type=str, help = "path/to/completeness.tsv")
    parser.add_argument("-o","--output_directory", type=str, default="checkm_filtered_output", help = "path/to/output_directory [Default: checkm_filtered_output]")
    parser.add_argument("-x","--extension", type=str, default="fa", help = "Fasta file extension for bins [Default: fa]")
    parser.add_argument("-m", "--minimum_contig_length", type=int, default=1, help="Minimum contig length. [Default: 1]")
    parser.add_argument("-u", "--unbinned", action="store_true", help="Write unbinned fasta sequences to file")
    parser.add_argument("-p", "--viral_prefix", type=str, default="Virus.", help = "Prefix for viral name (e.g. Sample_1__Virus. would make Sample_1__Virus.1, Sample_2__Virus.2, etc.) [Default: Virus.]")
    parser.add_argument("--include_provirus_detection", action="store_true", help="Include provirus viral detection")
    parser.add_argument("--multiplier_viral_to_host_genes", type=int, default=5, help = "Minimum number of viral genes [Default: 5]")
    parser.add_argument("--completeness", type=float, default=50.0, help = "Minimum completeness [Default: 50.0]")
    parser.add_argument("--checkv_quality", type=str, default="High-quality,Medium-quality,Complete", help = "Comma-separated string of acceptable arguments between {High-quality,Medium-quality,Complete} [Default: High-quality,Medium-quality,Complete]")
    parser.add_argument("--miuvig_quality", type=str, default="High-quality,Medium-quality,Complete", help = "Comma-separated string of acceptable arguments between {High-quality,Medium-quality,Complete} [Default: High-quality,Medium-quality,Complete]")
    parser.add_argument("--genomad_virus_taxonomy", type=str, help = "geNomad virus_taxonomy.tsv")
    parser.add_argument("--genomad_virus_summary", type=str, help = "geNomad virus_summary.tsv")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Output directory
    os.makedirs(opts.output_directory, exist_ok=True)
    os.makedirs(os.path.join(opts.output_directory,"genomes"), exist_ok=True)

    # Load CheckV
    df_checkv = pd.read_csv(opts.checkv_results, sep="\t", index_col=0)

    opts.checkv_quality = set(opts.checkv_quality.strip().split(","))
    assert opts.checkv_quality <= {"High-quality", "Medium-quality", "Complete"}, "Please choose some combination of these (comma separated list) {High-quality, Medium-quality, Complete}"
    opts.miuvig_quality = set(opts.miuvig_quality.strip().split(","))
    assert opts.miuvig_quality <= {"High-quality", "Medium-quality", "Complete"}, "Please choose some combination of these (comma separated list) {High-quality, Medium-quality, Complete}"

    df_completeness = pd.read_csv(opts.completeness_tsv, sep="\t", index_col=0)
    for field in ["contig_length", "proviral_length", "kmer_freq"]:
        if field in df_completeness.columns:
            df_completeness = df_completeness.drop([field], axis=1) # Duplicates
    df_completeness = pd.DataFrame(df_completeness.to_dict())

    if opts.include_provirus_detection:
        def filter_checkv(series):
            num_viral, num_host, completeness, checkv_quality, miuvig_quality, contig_length = series[["viral_genes","host_genes", "completeness", "checkv_quality", "miuvig_quality", "contig_length"]]
            conditions = [
                num_viral > 0,
                num_viral >= opts.multiplier_viral_to_host_genes*num_host,
                completeness >= opts.completeness,
                checkv_quality in opts.checkv_quality,
                miuvig_quality in opts.miuvig_quality,
                contig_length >= opts.minimum_contig_length,
            ]
            # print(series.name, num_viral, num_host, completeness, checkv_quality, miuvig_quality, "\t", all(conditions))

            return all(conditions)
    else:
        def filter_checkv(series):
            num_viral, num_host, completeness, checkv_quality, miuvig_quality, contig_length, provirus = series[["viral_genes","host_genes", "completeness", "checkv_quality", "miuvig_quality", "contig_length", "provirus"]]
            conditions = [
                num_viral > 0,
                num_viral >= opts.multiplier_viral_to_host_genes*num_host,
                completeness >= opts.completeness,
                checkv_quality in opts.checkv_quality,
                miuvig_quality in opts.miuvig_quality,
                contig_length >= opts.minimum_contig_length,
                str(provirus).strip() == "No",
            ]
            # print(series.name, num_viral, num_host, completeness, checkv_quality, miuvig_quality, "\t", all(conditions))

            return all(conditions)

    # CheckV Filtered Results
    mask = df_checkv.apply(filter_checkv, axis=1)
    df_checkv = df_checkv.loc[mask]
    df_checkv = pd.concat([ 
        df_checkv,
        df_completeness.loc[df_checkv.index],
    ], axis=1)


    # Contig lists
    f_binned_list = open(os.path.join(opts.output_directory, "binned.list"), "w")
    f_unbinned_list = open(os.path.join(opts.output_directory, "unbinned.list"), "w")
    if opts.unbinned:
        f_unbinned_fasta = open(os.path.join(opts.output_directory, "unbinned.fasta"), "w")
    else:
        f_unbinned_fasta = open(os.devnull, "w")
    if df_checkv.empty:
        with open(opts.fasta, "r") as f_fasta: # Use stdin?
            for header, seq in tqdm(SimpleFastaParser(f_fasta), "Extracting viral and unbinned contigs", unit=" contig"):
                id_scaffold = header.split(" ")[0]
                if len(seq) >= opts.minimum_contig_length:
                    print(id_scaffold, file=f_unbinned_list)
                    print(">{}\n{}".format(header, seq), file=f_unbinned_fasta)
        df_checkv = pd.DataFrame(columns=["id_contig"] + df_checkv.columns.tolist())

    else:
        # Quality assessment on MAGs

        scaffold_to_bin = dict()
        i = 1
        with open(opts.fasta, "r") as f_fasta: # Use stdin?
            for header, seq in tqdm(SimpleFastaParser(f_fasta), "Extracting viral and unbinned contigs", unit=" contig"):
                id_scaffold = header.split(" ")[0]

                if len(seq) >= opts.minimum_contig_length:
                    if id_scaffold in df_checkv.index:
                        id_virus = "{}{}".format(opts.viral_prefix, i)
                        with open(os.path.join(opts.output_directory, "genomes", "{}.fa".format(id_virus)), "w") as f_out:

                            # Contig to virus mapping
                            scaffold_to_bin[id_scaffold] = id_virus   
                            # Add to list
                            print(id_scaffold, file=f_binned_list)
                            print(">{} {}\n{}".format(id_scaffold, id_virus, seq), file=f_out)
                        i += 1
                    else:
                        print(id_scaffold, file=f_unbinned_list)
                        print(">{}\n{}".format(header, seq), file=f_unbinned_fasta)
        scaffold_to_bin = pd.Series(scaffold_to_bin)
        scaffold_to_bin.to_frame().to_csv(os.path.join(opts.output_directory, "scaffolds_to_bins.tsv"), sep="\t", header=None)

        with open(os.path.join(opts.output_directory, "bins.list"), "w") as f_bins:
            for id_virus in scaffold_to_bin.unique():
                print(id_virus, file=f_bins)

        # Add contig
        df_checkv.insert(loc=0, column="id_contig", value=df_checkv.index)

    # Build output table
    df_output = df_checkv 

    if opts.genomad_virus_taxonomy:
        df_genomad_taxonomy = pd.read_csv(opts.genomad_virus_taxonomy, sep="\t", index_col=0)
        df_genomad_taxonomy = df_genomad_taxonomy.reindex(df_checkv.index)
        mask = df_genomad_taxonomy["lineage"].isnull()
        df_genomad_taxonomy.loc[mask,"lineage"] = "Unclassified"
        df_genomad_taxonomy.loc[mask,"taxid"] = -1
        df_genomad_taxonomy["taxid"] = df_genomad_taxonomy["taxid"].astype(int)
        df_genomad_taxonomy.loc[mask,"agreement"] = 1.0
        df_genomad_taxonomy.loc[mask,"n_genes_with_taxonomy"] = 0
        df_genomad_taxonomy["n_genes_with_taxonomy"] = df_genomad_taxonomy["n_genes_with_taxonomy"].astype(int)

        df_output = pd.concat([df_output, df_genomad_taxonomy], axis=1)

    if opts.genomad_virus_summary:
        df_genomad_summary = pd.read_csv(opts.genomad_virus_summary, sep="\t", index_col=0)
        df_genomad_summary = df_genomad_summary.reindex(df_checkv.index)
        df_output = pd.concat([df_output, df_genomad_summary], axis=1)

    df_output.index = df_output.index.map(lambda x: scaffold_to_bin[x])
    df_output.index.name = "id_genome"
    df_output.to_csv(os.path.join(opts.output_directory, "checkv_results.filtered.tsv" ), sep="\t")

    f_binned_list.close()
    f_unbinned_list.close()
    f_unbinned_fasta.close()

if __name__ == "__main__":
    main()
    
                