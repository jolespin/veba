#!/usr/bin/env python3
import argparse
import os
import sys
import pyfastx


def parse_gff3_id_mapping(gff3_path):
    id_to_scaffold = {}
    scaffolds = set()
    with open(gff3_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            scaffold = cols[0]
            scaffolds.add(scaffold)
            for attr in cols[8].rstrip(";").split(";"):
                if attr.startswith("ID="):
                    id_to_scaffold[attr[3:]] = scaffold
                    break
    return id_to_scaffold, scaffolds


def rename_scaffold(name, prefix):
    return prefix + name


def rename_feature_id(feature_id, prefix, id_to_scaffold):
    scaffold = id_to_scaffold.get(feature_id)
    if scaffold is None:
        print(f"WARNING: ID '{feature_id}' not found in GFF3 mapping, "
              f"falling back to prefix-only", file=sys.stderr)
        return prefix + feature_id
    return prefix + scaffold + "_" + feature_id


def process_scaffold_fasta(in_path, out_path, prefix):
    fa = pyfastx.Fasta(in_path, build_index=False)
    with open(out_path, "w") as fout:
        for name, seq in fa:
            new_name = rename_scaffold(name, prefix)
            fout.write(f">{new_name}\n{seq}\n")


def process_seq_fasta(in_path, out_path, prefix, id_to_scaffold):
    headers = []
    fa = pyfastx.Fasta(in_path, build_index=False, full_name=True)
    with open(out_path, "w") as fout:
        for full_name, seq in fa:
            parts = full_name.split()
            renamed = [rename_feature_id(p, prefix, id_to_scaffold)
                       for p in parts]
            header = " ".join(renamed)
            headers.append(header)
            fout.write(f">{header}\n{seq}\n")
    return headers


def process_gff3(in_path, out_path, prefix, id_to_scaffold):
    with open(in_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                fout.write(line)
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                fout.write(line)
                continue
            cols[0] = rename_scaffold(cols[0], prefix)
            new_attrs = []
            for attr in cols[8].rstrip(";").split(";"):
                if attr.startswith("ID="):
                    old_id = attr[3:]
                    new_attrs.append(
                        "ID=" + rename_feature_id(old_id, prefix, id_to_scaffold))
                elif attr.startswith("Parent="):
                    old_parent = attr[7:]
                    new_attrs.append(
                        "Parent=" + rename_feature_id(old_parent, prefix,
                                                      id_to_scaffold))
                else:
                    new_attrs.append(attr)
            cols[8] = ";".join(new_attrs) + ";"
            fout.write("\t".join(cols) + "\n")


def validate_headers_match(protein_headers, cds_headers):
    if protein_headers == cds_headers:
        return True
    errors = []
    if len(protein_headers) != len(cds_headers):
        errors.append(
            f"Count mismatch: {len(protein_headers)} proteins vs "
            f"{len(cds_headers)} CDS entries")
    n = min(len(protein_headers), len(cds_headers))
    mismatches = 0
    for i in range(n):
        if protein_headers[i] != cds_headers[i]:
            if mismatches < 5:
                errors.append(
                    f"  Line {i+1}: protein='{protein_headers[i]}' "
                    f"vs cds='{cds_headers[i]}'")
            mismatches += 1
    if mismatches > 5:
        errors.append(f"  ... and {mismatches - 5} more mismatches")
    print("ERROR: CDS and protein headers do not match:", file=sys.stderr)
    for e in errors:
        print(e, file=sys.stderr)
    return False


def main():
    parser = argparse.ArgumentParser(
        description="Prepend a prefix to scaffold and gene IDs across "
                    "genomic asset files for multi-sample merging.")
    parser.add_argument("-i", "--prefix", required=True,
                        help="Prefix to prepend (include any delimiter, "
                             "e.g. 'sample-name__')")
    parser.add_argument("-a", "--assembly", required=True,
                        help="Input assembly/scaffold FASTA")
    parser.add_argument("-p", "--proteins", required=True,
                        help="Input protein FASTA")
    parser.add_argument("-c", "--cds", required=True,
                        help="Input CDS FASTA")
    parser.add_argument("-g", "--gff3", required=True,
                        help="Input GFF3 annotation file")
    parser.add_argument("-o", "--output_directory", required=True,
                        help="Output directory")
    parser.add_argument("--no-match-cds-protein", action="store_true",
                        default=False, dest="no_match_cds_protein",
                        help="Skip CDS/protein header validation "
                             "(validation is on by default)")
    args = parser.parse_args()

    for path in [args.assembly, args.proteins, args.cds, args.gff3]:
        if not os.path.isfile(path):
            print(f"ERROR: File not found: {path}", file=sys.stderr)
            sys.exit(1)

    os.makedirs(args.output_directory, exist_ok=True)

    print("Parsing GFF3 for ID-to-scaffold mapping...")
    id_to_scaffold, scaffolds = parse_gff3_id_mapping(args.gff3)
    print(f"  Found {len(scaffolds)} scaffolds and {len(id_to_scaffold)} "
          f"feature IDs")

    assembly_out = os.path.join(args.output_directory, os.path.basename(args.assembly))
    print(f"Processing assembly FASTA -> {assembly_out}")
    process_scaffold_fasta(args.assembly, assembly_out, args.prefix)

    proteins_out = os.path.join(args.output_directory, os.path.basename(args.proteins))
    print(f"Processing protein FASTA -> {proteins_out}")
    protein_headers = process_seq_fasta(args.proteins, proteins_out,
                                        args.prefix, id_to_scaffold)

    cds_out = os.path.join(args.output_directory, os.path.basename(args.cds))
    print(f"Processing CDS FASTA -> {cds_out}")
    cds_headers = process_seq_fasta(args.cds, cds_out, args.prefix,
                                    id_to_scaffold)

    gff3_out = os.path.join(args.output_directory, os.path.basename(args.gff3))
    print(f"Processing GFF3 -> {gff3_out}")
    process_gff3(args.gff3, gff3_out, args.prefix, id_to_scaffold)

    if not args.no_match_cds_protein:
        print("Validating CDS/protein header consistency...")
        if not validate_headers_match(protein_headers, cds_headers):
            sys.exit(1)
        print("  CDS and protein headers match.")

    print("Done.")


if __name__ == "__main__":
    main()
