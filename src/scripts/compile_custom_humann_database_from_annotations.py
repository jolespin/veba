#!/usr/bin/env python
from __future__ import print_function, division
from email import header
import sys, os, argparse, gzip
from collections import defaultdict
import numpy as np
import pandas as pd
from tqdm import tqdm 
from Bio.SeqIO.FastaIO import SimpleFastaParser

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.10.11"


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <identifier_mapping> -a <annotations> -t <taxnomy> -o <output_table>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--identifier_mapping", default="stdin", type=str,  help = "path/to/identifier_mapping.tsv[.gz] [id_protein]<tab>[id_genome] (No header) [Default: stdin]")
    parser.add_argument("-a","--annotations",  type=str,  required=True, help = "path/to/annotations.tsv[.gz] Output from annotations.py. Multi-level header that contains (UniRef, sseqid)")
    parser.add_argument("-t","--taxonomy",  type=str, required=True,  help = "path/to/taxonomy.tsv[.gz] [id_genome]<tab>[classification] (No header).  Use output from `merge_taxonomy_classifications.py` with --no_header and --no_domain")
    parser.add_argument("-s","--sequences",  type=str, required=True,  help = "path/to/proteins.fasta[.gz]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/humann_uniref_annotations.tsv[.gz] [Default: stdout]")
    parser.add_argument("--sep", default=";", help = "Separator for taxonomic levels [Default: ;]")
    # parser.add_argument("--mandatory_taxonomy_prefixes", help = "Comma-separated values for mandatory prefix levels. (e.g., 'c__,f__,g__,s__')")
    # parser.add_argument("--discarded_file",  help = "Proteins that have been discarded due to incomplete lineage")
    parser.add_argument("-g", "--no_append_genome_identifier", action="store_true", help = "Don't add genome to taxonomic lineage")
    parser.add_argument("--genome_prefix", type=str, default="t__", help = "Taxonomic level prefix for genome")
    parser.add_argument("--header", action="store_true", help = "Write header")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    if opts.identifier_mapping == "stdin":
        opts.identifier_mapping = sys.stdin 

    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Proteins to genomes
    protein_to_genome = pd.read_csv(opts.identifier_mapping, sep="\t", index_col=0, header=None).iloc[:,0]
    A1 = set(protein_to_genome.index)
    A2= set(protein_to_genome.values)
    print("--identifier_mapping", opts.identifier_mapping, file=sys.stderr)
    print(" * {} proteins".format(len(A1)), file=sys.stderr)
    print(" * {} genomes".format(len(A2)), file=sys.stderr)

    # Annotations
    df_annotations = pd.read_csv(opts.annotations, sep="\t", index_col=0, header=[0,1])
    assert "UniRef" in df_annotations.columns.get_level_values(0), "--annotations must have a 2 level header (i.e., Pandas MultiIndex with 2 levels) where the first level has 'UniRef' as created by `annotate.py`"
    df_annotations = df_annotations["UniRef"]
    protein_to_uniref = df_annotations["sseqid"].dropna()
    B1 = set(protein_to_uniref.index)
    B2 = set(protein_to_uniref.values)

    print("--annotations", opts.annotations, file=sys.stderr)
    print(" * {} proteins".format(len(B1)), file=sys.stderr)
    print(" * {} UniRef hits".format(len(B2)), file=sys.stderr)
    # Taxonomy
    genome_to_taxonomy = pd.read_csv(opts.taxonomy, sep="\t", index_col=0, header=None).iloc[:,0]
    C1 = set(genome_to_taxonomy.index)
    C2 = set(genome_to_taxonomy.values)

    print("--taxonomy", opts.taxonomy, file=sys.stderr)
    print(" * {} genomes".format(len(C1)), file=sys.stderr)
    print(" * {} taxonomic classifications".format(len(C2)), file=sys.stderr)

    if opts.sequences.endswith(".gz"):
        f = gzip.open(opts.sequences, "rt") 
    else:
        f = open(opts.sequences, "r")
    
    protein_to_length = dict()
    for header, seq in tqdm(SimpleFastaParser(f), "Calculating length of proteins: {}".format(opts.sequences), unit=" Proteins"):
        id = header.split(" ")[0]
        protein_to_length[id] = len(seq)
    protein_to_length = pd.Series(protein_to_length)
    D = set(protein_to_length.index)

    # Checks
    assert A1 >= B1, "Not all proteins in --annotations are in --identifier_mapping."
    assert A2 == C1, "Genomes in --identifier_mapping do not match genomes in --taxonomy.\n\nThe following genomes are specific to --identifier_mapping: {}\n\nThe following genomes are specific to --taxonomy".format("\n".join(A2 - C1), "\n".join(C1 - A2))
    assert B1 <= D, "Not all proteins in --annotations are in --sequences."

    # Append genome to taxonomy
    if not opts.no_append_genome_identifier:
        tmp = dict()
        for id_genome, taxonomy in genome_to_taxonomy.items():
            tmp[id_genome] = "{}{}{}{}".format(taxonomy, opts.sep, opts.genome_prefix, id_genome)
        genome_to_taxonomy = pd.Series(tmp)



    # id_protein, uniref_hit, len, lineage
    df_output = protein_to_uniref.to_frame("UniRef")
    df_output["Length"] = protein_to_length[protein_to_uniref.index]
    df_output["Taxonomy"] = protein_to_uniref.index.map(lambda id_protein: genome_to_taxonomy[protein_to_genome[id_protein]])
    df_output.index.name = "id_protein"


    # if opts.discarded_file:

    df_output.to_csv(opts.output, sep="\t", header=bool(opts.header))

if __name__ == "__main__":
    main()
