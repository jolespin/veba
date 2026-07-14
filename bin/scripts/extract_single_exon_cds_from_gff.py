#!/usr/bin/env python
import argparse
import sys
from contextlib import nullcontext

import pyfastx
from loguru import logger
from pyexeggutor import (
    open_file_reader,
    open_file_writer,
)


def parse_attributes(field):
    id_to_value = {}
    for item in field.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            id_to_value[key] = value
    return id_to_value


def reverse_complement(seq):
    complement_table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(complement_table)[::-1]


def main():
    parser = argparse.ArgumentParser(description="Extract nucleotide CDS sequences from a flat GFF3 file (no Parent/exon hierarchy assumed).")
    parser.add_argument("-g", "--gff", required=True, help="Input GFF3 file (.gz supported)")
    parser.add_argument("-f", "--genome", required=True, help="Input genome FASTA file (.gz supported)")
    parser.add_argument("-o", "--output", default=None, help="Output CDS FASTA file (.gz supported) [default: stdout]")
    parser.add_argument("-t", "--feature_type", default="CDS", help="GFF feature type to extract (default: CDS)")
    opts = parser.parse_args()

    logger.info(f"Loading genome: {opts.genome}")
    fasta_genome = pyfastx.Fasta(opts.genome, build_index=False)
    contig_to_seq = {name: seq for name, seq in fasta_genome}

    logger.info(f"Parsing GFF: {opts.gff}")
    n_written = 0
    n_skipped_missing_contig = 0

    if opts.output:
        output_writer = open_file_writer(opts.output)
    else:
        output_writer = nullcontext(sys.stdout)

    with output_writer as f_out:
        with open_file_reader(opts.gff) as f_gff:
            for line in f_gff:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue

                seqid, source, feature_type, start, end, score, strand, frame, attributes = fields
                if feature_type != opts.feature_type:
                    continue

                if seqid not in contig_to_seq:
                    n_skipped_missing_contig += 1
                    continue

                id_to_value = parse_attributes(attributes)
                cds_id = id_to_value.get("ID", f"{seqid}_{start}_{end}")

                seq = contig_to_seq[seqid][int(start) - 1:int(end)]
                if strand == "-":
                    seq = reverse_complement(seq)

                print(f">{cds_id}", file=f_out)
                print(seq, file=f_out)
                n_written += 1

    output_destination = opts.output if opts.output else "stdout"
    logger.info(f"Wrote {n_written} CDS sequences to {output_destination}")
    if n_skipped_missing_contig:
        logger.warning(f"Skipped {n_skipped_missing_contig} CDS features whose seqid was not found in the genome FASTA")


if __name__ == "__main__":
    main()