#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, gzip
from collections import OrderedDict
import pandas as pd
from tqdm import tqdm 
from Bio.SeqIO.FastaIO import SimpleFastaParser


pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.6.28"


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -f <genome.fasta> -t <tiara_results.tsv> -n <genome_name> -o <output_directory> -x fa".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline

    parser.add_argument("-f","--fasta",  type=str, default="stdin", help = "path/to/genome_assembly.fasta [Default: stdin]")
    parser.add_argument("-t","--tiara_results",  type=str, required=True, help = "path/to/tiara_results.tsv")
    parser.add_argument("-n","--name",  type=str,  required=True, help = "Name of genome")
    parser.add_argument("-o","--output_directory", type=str, required=True, help = "Output directory")
    parser.add_argument("-u","--unknown_organelle_prediction", type=str, default="nuclear", help = "{unknown, nuclear} [Default: nuclear]")
    parser.add_argument("-x","--extension", type=str, default="fa", help = "Fasta output extension [Default: fa]")
    parser.add_argument("--mitochondrion_suffix", type=str, default="", help = "Mitochondrion suffix [Default: '']")
    parser.add_argument("--plastid_suffix", type=str, default="", help = "Plastid suffix [Default: '']")
    parser.add_argument("--unknown_suffix", type=str, default="", help = "Unknown suffix [Default: '']")
    parser.add_argument("--flat_directory", action="store_true", help = "Use if you don't want a nested directory and for all the files to be flatten in the output directory")
    parser.add_argument("--no_description", action="store_true", help = "Default is to add [GENOME_NAME]:[nuclear/mitochondrion/plastid] as the description")
    parser.add_argument("--no_empty_files", action="store_true", help = "Default is to create file even if it's empty")
    parser.add_argument("--verbose", action="store_true", help = "Print ids and lengths for seqeunces not in Tiara")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert opts.unknown_organelle_prediction in {"unknown", "nuclear"}, "--unknown_organelle_prediction must be one of the following:{unknown, nuclear}"

    # Organelle predicitions
    df_tiara = pd.read_csv(opts.tiara_results, sep="\t", index_col=0)[["class_fst_stage", "class_snd_stage"]]
    df_tiara.index = df_tiara.index.map(lambda x: x.split(" ")[0])
    
    id_to_classification = dict()
    for id, (fst_prediction, snd_prediction) in df_tiara.iterrows():
        if fst_prediction in {"archaea", "bacteria", "eukarya", "unknown"}:
            id_to_classification[id] = "nuclear"
        else:
            if snd_prediction in {"mitochondrion", "plastid"}:
                id_to_classification[id] = snd_prediction
            else:
                id_to_classification[id] = opts.unknown_organelle_prediction
    

    # Output directories
    os.makedirs(opts.output_directory, exist_ok=True)
    output_files = dict()
    if opts.no_empty_files:
        seq_types = set(id_to_classification.values())
    else:
        seq_types = {"nuclear", "mitochondrion", "plastid", opts.unknown_organelle_prediction}

    
    for id_type in seq_types:
        if id_type == "nuclear":
            filepath = os.path.join(opts.output_directory, "{}.{}".format(opts.name, opts.extension))
        else:
            suffix = {"mitochondrion":opts.mitochondrion_suffix, "plastid":opts.plastid_suffix, "unknown":opts.unknown_suffix}[id_type]

            if opts.flat_directory:
                filepath = os.path.join(opts.output_directory, "{}{}.{}".format(opts.name, suffix, opts.extension))
            else:
                os.makedirs(os.path.join(opts.output_directory, id_type), exist_ok=True)
                filepath = os.path.join(opts.output_directory, id_type, "{}{}.{}".format(opts.name, suffix, opts.extension))
        output_files[id_type] = open(filepath, "w")


    # Partition fasta
    fasta_filepath = "stdin"
    if opts.fasta == "stdin":
        opts.fasta = sys.stdin 
    else:
        fasta_filepath = opts.fasta
        if opts.fasta.endswith(".gz"):
            opts.fasta = gzip.open(opts.fasta, "rt")
        else:
            opts.fasta = open(opts.fasta, "r")

    for header, seq in tqdm(SimpleFastaParser(opts.fasta), "Partitioning genome assembly: {}".format(fasta_filepath)):
        id = header.split(" ")[0]
        if id in id_to_classification:
            id_type = id_to_classification[id]
        else:
            id_type = "nuclear"
            if opts.verbose:
                print("Identifier '{}' of length {} was not in Tiara results and will be tagged as 'nuclear'".format(id, len(seq)), file=sys.stdout)
        file = output_files[id_type]
        if opts.no_description:
            print(">{}\n{}".format(id,seq), file=file)
        else:
            description = "{}:{}".format(opts.name, id_type)
            print(">{} {}\n{}".format(id, description, seq), file=file)

    if opts.fasta != sys.stdin:
        opts.fasta.close()

    for id_type, file in output_files.items():
        if opts.verbose:
            print("Closing file type: {}".format(id_type), file=sys.stdout)
        file.close()


if __name__ == "__main__":
    main()
