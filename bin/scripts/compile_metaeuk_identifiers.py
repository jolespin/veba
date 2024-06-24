#!/usr/bin/env python
import sys, os, glob, argparse 
from collections import OrderedDict
import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.6.20"

def parse_header(header:str, include_strand_in_geneid=True, strand_notation="+/-"):
    """
    This function has been adapted from BUSCO v5. 
    
    Source code: 
    * https://gitlab.com/ezlab/busco/-/blob/master/src/busco/busco_tools/metaeuk.py

    If this is used, please cite BUSCO:
    Mosè Manni, Matthew R Berkeley, Mathieu Seppey, Felipe A Simão, Evgeny M Zdobnov, 
    BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper 
    Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes, 
    Molecular Biology and Evolution, Volume 38, Issue 10, October 2021, 
    Pages 4647–4654, https://doi.org/10.1093/molbev/msab199
    """
    assert strand_notation in {"+/-", "1/-1"}
    header_parts = header.split("|")
    if not header_parts[2] in [
        "+",
        "-",
    ]:  # Deal with sequence IDs that contain the symbol "|"
        try:
            strand_ind = header_parts.index("+")
        except ValueError:
            strand_ind = header_parts.index("-")
        header_parts[0] = "|".join(header_parts[:strand_ind-1])
        for i in range(1, strand_ind-1):
            header_parts.pop(i)

    T_acc = header_parts[0]
    C_acc = header_parts[1]
    strand = header_parts[2]
    bitscore = float(header_parts[3])
    evalue = float(header_parts[4])
    num_exons = int(header_parts[5])
    low_coord = int(header_parts[6])
    high_coord = int(header_parts[7])
    exon_coords = header_parts[8:]

    all_low_exon_coords = []
    all_taken_low_exon_coords = []
    all_high_exon_coords = []
    all_taken_high_exon_coords = []
    all_exon_nucl_len = []
    all_taken_exon_nucl_len = []
    for exon in exon_coords:
        low_exon_coords, high_exon_coords, nucl_lens = exon.split(":")

        low_exon_coord, taken_low_exon_coord = low_exon_coords.split("[")
        all_low_exon_coords.append(int(low_exon_coord))
        taken_low_exon_coord = int(taken_low_exon_coord.strip("]"))

        high_exon_coord, taken_high_exon_coord = high_exon_coords.split("[")
        all_high_exon_coords.append(int(high_exon_coord))
        taken_high_exon_coord = int(taken_high_exon_coord.strip("]"))

        nucl_len, taken_nucl_len = nucl_lens.split("[")
        all_exon_nucl_len.append(int(nucl_len))
        taken_nucl_len = int(taken_nucl_len.strip().rstrip("]"))

        # Need to fix the metaeuk coordinate problem
        if strand == "-":
            if int(taken_high_exon_coord) + int(taken_nucl_len) - 1 != int(
                taken_low_exon_coord
            ):
                taken_low_exon_coord = (
                    int(taken_high_exon_coord) + int(taken_nucl_len) - 1
                )

        all_taken_low_exon_coords.append(taken_low_exon_coord)
        all_taken_high_exon_coords.append(taken_high_exon_coord)
        all_taken_exon_nucl_len.append(taken_nucl_len)
    
    
    if include_strand_in_geneid:
        if strand_notation == "+/-":
            gene_id = "{}_{}:{}({})".format(
                C_acc, 
                low_coord, 
                high_coord,
                strand,
            )
        if strand_notation == "1/-1":
            gene_id = "{}_{}:{}({})".format(
                C_acc, 
                low_coord, 
                high_coord,
                {"+":"1","-":"-1"}[strand], 
            )
    else:
        gene_id = "{}_{}:{}".format(
            C_acc, 
            low_coord, 
            high_coord,
        )
    

    details = {
        "MetaEuk_header":header,
        "T_acc": T_acc,
        "C_acc": C_acc,
        "S": strand,
        "gene_id": gene_id,
        "bitscore": bitscore,
        "e-value": evalue,
        "num_exons": num_exons,
        "low_coord": low_coord,
        "high_coord": high_coord,
        "all_low_exon_coords": all_low_exon_coords,
        "all_taken_low_exon_coords": all_taken_low_exon_coords,
        "all_high_exon_coords": all_high_exon_coords,
        "all_taken_high_exon_coords": all_taken_high_exon_coords,
        "all_exon_nucl_len": all_exon_nucl_len,
        "all_taken_exon_nucl_len": all_taken_exon_nucl_len,
    }
    return details


def main(argv=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -d <output.codon.fas> -a <output.fas> -o <output_directory>".format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    
    # I/O
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-d","--cds", type=str, required=True, help = "path/to/output.codons.fas [Required]")
    parser_io.add_argument("-a","--protein", type=str, required=True, help = "path/to/output.fas [Required]")
    parser_io.add_argument("-o","--output_directory", type=str, default="metaeuk_output", help = "path/to/output_directory [Default: metaeuk_output]")
    parser_io.add_argument("-b", "--basename", type=str, default="gene_models", help = "Base name for gene models [Default: gene_models]")
    parser_io.add_argument("--cds_extension", type=str, default="ffn", help = "CDS fasta extension [Default: ffn]")
    parser_io.add_argument("--protein_extension", type=str, default="faa", help = "Protein fasta extension [Default: faa]")
    parser_io.add_argument("--no_header", action="store_true", help="Specify if MetaEuk header should not be in fasta record as description")

    # GFF
    parser_gff = parser.add_argument_group('GFF arguments')
    parser_gff.add_argument("-e", "--exon_coordinate_type", type=str, default="low/high", help = "Exon coordinate type {'low/high', 'taken_low/high'} [Default: low/high]")
    parser_gff.add_argument("-m", "--include_mrna", action="store_true", help="Include mRNA feature in GFF output")

    # GeneID
    parser_geneid = parser.add_argument_group('Gene identifier arguments')
    # parser_geneid.add_argument("--geneid_delimiter", type=str, default=":", help = "Delimiter for simplified identifiers [id_contig]_[gene_start]<delimiter>[gene_end] [Default: : ]")
    parser_geneid.add_argument("--strand_notation", type=str, default="+/-", help = "Strand notation fo gene_id {'+/-','1/-1'} [Default: '+/-']")
    parser_geneid.add_argument("--no_strand", action="store_true", help="Don't include strand in gene id")

    # Assertions
    parser_assertions = parser.add_argument_group('Assertion arguments')
    parser_assertions.add_argument("--assert_identifiers_are_same_order", action="store_true", help="Assert that --cds and --protein are the same order")
    parser_assertions.add_argument("--assert_no_duplicate_headers", action="store_true", help="If there is an exact overlap of gene positions and strand, the gene with highest bitscore is kept. Choose this if you want to assert there are no duplicates.")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Assertions 
    assert opts.exon_coordinate_type  in {"low/high", "taken_low/high"}
    assert opts.strand_notation  in {"+/-", "1/-1"}

    # Output directory
    os.makedirs(opts.output_directory, exist_ok=True)

    # Load in CDS sequences
    id_to_seq = OrderedDict()
    with open(opts.cds, "r") as f:
        for header, seq in SimpleFastaParser(f):
            id = header.split(" ")[0]
            id_to_seq[id] = seq
    cds_sequences = pd.Series(id_to_seq)

    # Load in protein sequences
    id_to_seq = OrderedDict()
    with open(opts.protein, "r") as f:
        for header, seq in SimpleFastaParser(f):
            id = header.split(" ")[0]
            id_to_seq[id] = seq
    protein_sequences = pd.Series(id_to_seq)

    # Check if CDS sequence identifiers is the same as protein identifiers.
    # If not, then take the union
    cds_and_protein_identifiers_are_the_same_order = True
    if cds_sequences.size == protein_sequences.size:
        if not np.all(cds_sequences.index == protein_sequences.index):
            cds_and_protein_identifiers_are_the_same_order = False 
            if opts.assert_identifiers_are_same_order:
                raise AssertionError("--cds and --protein do not have the same order of identifiers")
                sys.exit(1)
    else:
        cds_and_protein_identifiers_are_the_same_order = False 
        if opts.assert_identifiers_are_same_order:
            raise AssertionError("--cds and --protein are different sizes and, thus, cannot be the same order of identifiers")
            sys.exit(1)
    if cds_and_protein_identifiers_are_the_same_order:
        metaeuk_identifiers = cds_sequences.index
    else:
        metaeuk_identifiers = pd.Index(sorted(set(cds_sequences.index) | set(protein_sequences.index)))

    # Parse MetaEuk headers and simplify gene identifiers
    df_metaeuk_headers = pd.DataFrame(list(map(lambda header: parse_header(
        header=header, 
        strand_notation=opts.strand_notation, 
        include_strand_in_geneid=(not bool(opts.no_strand)),
    ), 
    tqdm(metaeuk_identifiers, desc="Parsing MetaEuk headers", total=len(metaeuk_identifiers), unit=" gene")))).set_index("MetaEuk_header")

    # Handle duplicates
    geneid_value_counts = df_metaeuk_headers["gene_id"].value_counts()
    if geneid_value_counts.max() > 1:
        number_of_genes_with_duplicate_identifiers = np.sum(geneid_value_counts.values > 1)
        genes_with_duplicate_identifiers = geneid_value_counts[geneid_value_counts > 1].index.tolist()
        print("There are {} duplicate gene identifiers:\n{}".format(
            number_of_genes_with_duplicate_identifiers,
            "\n".join(genes_with_duplicate_identifiers),
            ), file=sys.stderr)

        
        if opts.assert_no_duplicate_headers:
            raise AssertionError("There are {} duplicate gene identifiers:\n{}".format( 
                number_of_genes_with_duplicate_identifiers,
                "\n".join(genes_with_duplicate_identifiers),
                ))
        # Get indicies of highest bit scores
        index_highest_bitscores = df_metaeuk_headers.groupby("gene_id")["bitscore"].idxmax()
        # Remove duplicates from MetaEuk headers
        df_metaeuk_headers = df_metaeuk_headers.loc[index_highest_bitscores]

        # Remove duplicates from sequences
        cds_sequences = cds_sequences[cds_sequences.index.intersection(df_metaeuk_headers.index)]
        protein_sequences = protein_sequences[protein_sequences.index.intersection(df_metaeuk_headers.index)]



    # Convert MetaEuk to simple identifiers
    metaeuk_to_simple = df_metaeuk_headers["gene_id"]

    # Relabel identifiers
    if opts.no_header:
        cds_sequences.index = cds_sequences.index.map(lambda id_metaeuk: metaeuk_to_simple[id_metaeuk])
        protein_sequences.index = protein_sequences.index.map(lambda id_metaeuk: metaeuk_to_simple[id_metaeuk])
    else:
        cds_sequences.index = cds_sequences.index.map(lambda id_metaeuk: "{} {}".format(metaeuk_to_simple[id_metaeuk], id_metaeuk))
        protein_sequences.index = protein_sequences.index.map(lambda id_metaeuk: "{} {}".format(metaeuk_to_simple[id_metaeuk], id_metaeuk))

    # Exon coord fields
    if opts.exon_coordinate_type == "low/high":
        low_exon_coord_field = "all_low_exon_coords"
        high_exon_coord_field = "all_high_exon_coords"
    if opts.exon_coordinate_type == "taken_low/high":
        low_exon_coord_field = "all_taken_low_exon_coords"
        high_exon_coord_field = "all_taken_high_exon_coords"
        
    # Create GFF
    gff_output = list()
    for id_metaeuk, fields in tqdm(df_metaeuk_headers.iterrows(), desc="Creating gene, mRNA, and CDS records for GFF", unit=" genes", total=df_metaeuk_headers.shape[0]):
        try:

            # Gene
            id_target = fields["T_acc"]
            id_contig = fields["C_acc"]
            strand = fields["S"]
            id_tcs = "|".join([id_target, id_contig, strand])
            id_gene = fields["gene_id"]

            start_gene = int(fields["low_coord"]) + 1
            end_gene = int(fields["high_coord"]) + 1
            bitscore = float(fields["bitscore"])
            gene_description = "target_id={};tcs_id={};contig_id={};gene_id={};ID={};".format(
                id_target,
                id_tcs,
                id_contig,
                id_gene,
                id_gene,
            )
            gene_fields = [id_contig, "MetaEuk", "gene", start_gene, end_gene, bitscore, strand, ".", gene_description]

            # gff_output.append("\t".join(map(str, gene_fields)))
            gff_output.append(gene_fields)

            # mRNA
            mrna_description = "target_id={};tcs_id={};contig_id={};gene_id={};ID={};Parent={};".format(
                id_target,
                id_tcs,
                id_contig,
                id_gene,
                id_gene,
                id_gene,
            )
            mrna_fields = [id_contig, "MetaEuk", "mRNA", start_gene, end_gene, bitscore, strand, ".", mrna_description]
            if opts.include_mrna:
                gff_output.append(mrna_fields)

            # CDS
            cds_fields = [id_contig, "MetaEuk", "CDS", start_gene, end_gene, bitscore, strand, ".", mrna_description]

            # gff_output.append("\t".join(map(str, mrna_fields)))
            gff_output.append(cds_fields)

            # Exons
            for i, (start_exon, end_exon) in enumerate(zip(fields[low_exon_coord_field], fields[high_exon_coord_field]), start=1):
                start_exon = start_exon + 1
                end_exon = end_exon + 1
                id_exon = "{}.exon_{}".format(id_gene, i)

                exon_description = "target_id={};tcs_id={};contig_id={};gene_id={};ID={};Parent={};exon_id={}".format(
                    id_target,
                    id_tcs,
                    id_contig,
                    id_gene,
                    id_gene,
                    id_gene,
                    id_exon,
                )
                exon_fields = [id_contig, "MetaEuk", "exon", start_exon, end_exon, bitscore, strand, ".", exon_description]
                gff_output.append(exon_fields)
                # output.append("\t".join(map(str, exon_fields)))
        except Exception as e:
            print(e, file=sys.stderr)
            print(id_metaeuk, file=sys.stderr)
            print(fields, file=sys.stderr)
            # sys.exit(1)
    df_gff = pd.DataFrame(gff_output)


    # Write output

    # Identifiers
    df_metaeuk_headers.to_csv(os.path.join(opts.output_directory,"identifier_mapping.metaeuk.tsv"), sep="\t")

    metaeuk_to_simple.to_frame().to_csv(os.path.join(opts.output_directory,"metaeuk_to_simple.tsv"), sep="\t", header=False)

    # CDS
    with open(os.path.join(opts.output_directory,"{}.{}".format(opts.basename, opts.cds_extension)) , "w") as f:
        f.writelines(">{}\n{}\n".format(id, seq) for id, seq in cds_sequences.to_dict(into=OrderedDict).items())

    # Protein
    with open(os.path.join(opts.output_directory,"{}.{}".format(opts.basename, opts.protein_extension)) , "w") as f:
        f.writelines(">{}\n{}\n".format(id, seq) for id, seq in protein_sequences.to_dict(into=OrderedDict).items())

    # GFF
    df_gff.to_csv(os.path.join(opts.output_directory,"{}.gff".format(opts.basename)), sep="\t", index=False, header=False)
    
if __name__ == "__main__":
    main()
    
                

