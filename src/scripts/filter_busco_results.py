#!/usr/bin/env python
import sys, os, glob, argparse
from shutil import copyfile
from collections import OrderedDict
import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.04.04"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <busco_output.tsv> -g <metaeuk_directory> -o <output_directory> -f <scaffolds.fasta> -m 1500 -u".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-i","--busco_results", type=str, required=True, help = "path/to/busco_output/[id_mag]/busco_results.tsv")
    parser.add_argument("-g","--metaeuk_directory", type=str, required=True, help = "path/to/metaeuk/directory")
    parser.add_argument("-o","--output_directory", type=str, default="busco_filtered_output", help = "path/to/output_directory [Default: busco_filtered_output]")
    parser.add_argument("-f", "--fasta", type=str, help = "path/to/scaffolds.fasta. Include this only if you want to list of unbinned contigs")
    parser.add_argument("-m", "--minimum_contig_length", type=int, default=1, help="Minimum contig length. [Default: 1]")
    parser.add_argument("-u", "--unbinned", action="store_true", help="Write unbinned fasta sequences to file")
    parser.add_argument("--completeness", type=float, default=50.0, help = "BUSCO completeness [Default: 50.0]")
    parser.add_argument("--contamination", type=float, default=10.0, help = "BUSCO contamination [Default: 10.0]")

    # parser.add_argument("--genome_extension", type=str, default="fa", help = "Fasta file extension for bins [Default: fa]")
    # parser.add_argument("--protein_extension", type=str, default="faa", help = "Fasta file extension for protein sequences [Default: faa]")
    # parser.add_argument("--cds_extension", type=str, default="ffn", help = "Fasta file extension for cds sequences [Default: ffn]")
    # parser.add_argument("--genemodel_extension", type=str, default="ffn", help = "Fasta file extension for gene models [Default: gff]")

    parser.add_argument("--symlink", action="store_true", help = "Symlink MAGs instead of copying")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Output filtered 
    os.makedirs(opts.output_directory, exist_ok=True)

    # Load BUSCO
    df_busco = pd.read_csv(opts.busco_results, sep="\t", index_col=0, header=[0,1])

    # Genome to completeness
    mag_to_completeness = dict()
    mag_to_contamination = dict()
    # Iterate through all of the genomes
    for id_mag, busco_results in df_busco.iterrows():
        # Is specific or generic available
        busco_specificities = set(busco_results.dropna().index.get_level_values(0))
        # Prefer 'specific' over 'generic'
        if "specific" in busco_specificities:
            busco_notation = busco_results["specific"]["one_line_summary"]

        else:
            busco_notation = busco_results["generic"]["one_line_summary"]
        mag_to_completeness[id_mag] = float(busco_notation.split("%[")[0].split(":")[1])
        mag_to_contamination[id_mag] = float(busco_notation.split("%]")[0].split(":")[-1])
    mag_to_completeness = pd.Series(mag_to_completeness)
    mag_to_contamination = pd.Series(mag_to_contamination)

    df_quality = pd.concat([mag_to_completeness.to_frame("completeness"), mag_to_contamination.to_frame("contamination")], axis=1)
    df_quality.to_csv(os.path.join(opts.output_directory,"busco_results.completeness_contamination.tsv"), sep="\t")


    # Quality Control 
    tmp = mag_to_completeness.copy()
    mag_to_completeness = mag_to_completeness[lambda c: c >= opts.completeness]

    if mag_to_completeness.empty:
        print("No bins had a completeness ≥ {}".format(opts.completeness), file=sys.stderr)
        tmp.to_frame("completeness").to_csv(sys.stderr, sep="\t")
        sys.exit(1)
    else: 
        del tmp

    tmp = mag_to_contamination.copy()
    mag_to_contamination = mag_to_contamination[lambda c: c < opts.contamination]

    if mag_to_contamination.empty:
        print("No bins had a contamination < {}".format(opts.contamination), file=sys.stderr)
        tmp.to_frame("contamination").to_csv(sys.stderr, sep="\t")
        sys.exit(1)
    else: 
        del tmp

    genomes_passed_qc = mag_to_completeness.index.intersection(mag_to_contamination.index)


    df_busco.loc[genomes_passed_qc,:].to_csv(os.path.join(opts.output_directory,"busco_results.filtered.tsv"), sep="\t")
    df_quality.loc[genomes_passed_qc,:].to_csv(os.path.join(opts.output_directory,"busco_results.completeness_contamination.filtered.tsv"), sep="\t")

    # Output MAGs
    os.makedirs(os.path.join(opts.output_directory,"genomes"), exist_ok=True)

    # bins.list
    with open(os.path.join(opts.output_directory, "bins.list"), "w") as f_bins:
        for id_mag in sorted(genomes_passed_qc):
            print(id_mag, file=f_bins)


    # Binned
    f_binned_list = open(os.path.join(opts.output_directory, "binned.list"), "w")

    binned_contigs = list() 
    scaffolds_to_bins = OrderedDict()
    for id_mag in tqdm(genomes_passed_qc, "Copying fasta files and writing binned contigs", unit=" MAG"):

        for src in glob.glob(os.path.join(opts.metaeuk_directory, "genomes", "{}.*".format(id_mag))):
            # Copy or symlink MAG
            # src = os.path.join(opts.mag_directory, "{}.{}".format(id_mag, opts.extension))
            fn = os.path.split(src)[-1]
            dst = os.path.join(opts.output_directory,"genomes", fn)

            if opts.symlink:
                if os.path.exists(dst):
                    os.remove(dst)
                os.symlink(os.path.realpath(src), dst)
            else:
                copyfile(src,dst)

        # Write contigs to list
        with open(os.path.join(opts.metaeuk_directory, "genomes", "{}.fa".format(id_mag)), "r") as f_fasta:
            for line in f_fasta.readlines():
                line = line.strip()
                if line.startswith(">"):
                    id_contig = line[1:]
                    scaffolds_to_bins[id_contig] = id_mag
                    print(id_contig, file=f_binned_list)
                    binned_contigs.append(id_contig)

    f_binned_list.close()
    
    # scaffolds_to_bins.tsv
    scaffolds_to_bins = pd.Series(scaffolds_to_bins)
    scaffolds_to_bins.to_frame().to_csv(os.path.join(opts.output_directory, "scaffolds_to_bins.tsv"), sep="\t", header=None)



    # Identifier Mapping
    f_identifiers = open(os.path.join(opts.output_directory,"genomes", "identifier_mapping.tsv"),"w")

    with open(os.path.join(opts.metaeuk_directory,"genomes", "identifier_mapping.tsv"), "r") as f_in:
        for line in f_in:
            line = line.strip()
            if line:
                id_orf, id_scaffold, id_mag = line.split("\t")
                if id_mag in genomes_passed_qc:
                    print(line, file=f_identifiers)
    f_identifiers.close()


    # identifier_mapping.metaeuk.tsv
    df_metaeuk_identifiers = pd.read_csv(os.path.join(opts.metaeuk_directory, "identifier_mapping.metaeuk.tsv"), sep="\t", index_col=0)
    mask = df_metaeuk_identifiers["C_acc"].map(lambda x: x in binned_contigs).values
    df_metaeuk_identifiers.loc[mask].to_csv(os.path.join(opts.output_directory, "identifier_mapping.metaeuk.tsv"), sep="\t")

    # metaeuk_to_simple.tsv
    metaeuk_to_simple = pd.read_csv(os.path.join(opts.metaeuk_directory,"metaeuk_to_simple.tsv"), sep="\t", index_col=0, header=None).iloc[:,0]
    mask = metaeuk_to_simple.map(lambda x: "_".join(x.split("_")[:-1]) in binned_contigs).values
    metaeuk_to_simple[mask].to_frame().to_csv(os.path.join(opts.output_directory, "metaeuk_to_simple.tsv"), sep="\t", header=None)

    # Get unbinned contigs
    if opts.fasta:

        from Bio.SeqIO.FastaIO import SimpleFastaParser

        f_unbinned_list = open(os.path.join(opts.output_directory, "unbinned.list"), "w")

        with open(opts.fasta, "r") as f_fasta: # Use stdin?
            for header, seq in tqdm(SimpleFastaParser(f_fasta), "Extracting unbinned contigs", unit=" contigs"):
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
                for header, seq in tqdm(SimpleFastaParser(f_fasta), "Writing unbinned contigs", unit=" contigs"):
                    id_contig = header.split(" ")[0]
                    conditions = [
                        id_contig not in binned_contigs,
                        len(seq) >= opts.minimum_contig_length,
                    ]
                    if all(conditions):
                        print(">{}\n{}".format(header, seq), file=f_unbinned_fasta)
            f_unbinned_fasta.close()

    

if __name__ == "__main__":
    main()
    
                