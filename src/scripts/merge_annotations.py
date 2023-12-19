#!/usr/bin/env python
import sys, os, argparse, re, gzip
from collections import defaultdict, OrderedDict
import pandas as pd
import numpy as np
from soothsayer_utils import read_hmmer, pv, get_file_object, assert_acceptable_arguments, format_header, flatten

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.11.24"

# disclaimer = format_header("DISCLAIMER: Lineage predictions are NOT robust and DO NOT USE CORE MARKERS.  Please only use for exploratory suggestions.")

# blacklist_words = {"multispecies:", "multispecies"}
# unknown_labels = {"predicted protein", "uncharacterized", "unknown", "hypothetical protein", "hypothetical protein, partial"}
# predicted_label="Predicted protein"

# def process_nr_name(name, predicted_label):
#     species = None
#     words = name.split(" ")
#     if "." == words[0][-2]:
#         words = words[1:]
#     for i, word in enumerate(words):
#         if word.lower() in blacklist_words:
#             words.pop(i)

#     name = " ".join(words)

#     matches = re.findall("\[.*?\]", name)
#     if matches:
#         for match in matches:
#             if match[1:-1].count(" ") == 1:
#                 species = match[1:-1]
#             name = name.replace(match,"")
#     if species:
#         name = name.replace(" - {}".format(species), "")
#         name = name.replace(species, "")

#     name = name.strip().strip("-")
#     if name in unknown_labels:
#         name = predicted_label
#     return name

def clean_stitle(stitle):
    sseqid, *stitle = stitle.split(" ")
    stitle = " ".join(stitle)
    return stitle 

def parse_uniref(stitle):

    # Define the regular expression pattern
    pattern = r'^(.*?)\s+n=(\d+)\s+Tax=.*?\s+TaxID=(\d+)\s+RepID=(.*)$'

    # Use re.match to search for the pattern in the string
    match = re.match(pattern, stitle)

    # Extract the four fields from the match object
    field1 = match.group(1)
    field1 = np.nan if field1 is None else field1
    field2 = match.group(2)
    field2 = np.nan if field2 is None else field2
    field3 = match.group(3)
    field3 = np.nan if field3 is None else field3
    field4 = match.group(4)
    field4 = np.nan if field4 is None else field4

    # Print the four fields
    return [field1, field2, field3, field4]

def compile_identifiers(df, id_protein_cluster):
    # Organism type
    organism_types = set(df["organism_type"])
    if len(organism_types) == 1:
        organism_types = list(organism_types)[0]

    # # Genomes
    # genomes = set(df["id_genome"])

    # # Samples
    # samples = set(df["sample_of_origin"])

    # Genome clusters
    genome_clusters = set(df["id_genome_cluster"])
    if len(genome_clusters) == 1:
        genome_clusters = list(genome_clusters)[0]

    # data = OrderedDict([
    #     ("id_genome_cluster", genome_clusters),
    #     ("organism_type", organism_types),
    #     ("genomes", genomes),
    #     ("samples_of_origin", samples),
    # ],
    # )
    # data = pd.Series(data, name=id_protein_cluster)
    # data.index = data.index.map(lambda x: ("Identifiers", x))
    data = [genome_clusters, organism_types]#, genomes, samples]
    return data

def compile_uniref(df, id_protein_cluster):
    df = df.dropna(how="all", axis=0)
    unique_identifiers = list(df["sseqid"].unique())
    # data = OrderedDict([
    #     ("number_of_proteins", df.shape[0]),
    #     ("number_of_unique_hits", len(unique_identifiers)),
    #     ("ids", unique_identifiers),
    #     ("names", set(df["product"].unique())),
    # ],
    # )
    # data = pd.Series(data, name=id_protein_cluster)
    # data.index = data.index.map(lambda x: ("UniRef", x))

    data = [df.shape[0], len(unique_identifiers), unique_identifiers, list(df["product"].unique())]
    return data

def compile_nonuniref_diamond(df, id_protein_cluster, label):
    df = df.dropna(how="all", axis=0)
    unique_identifiers = set(df["sseqid"].unique())
    # data = OrderedDict(
    #     [
    #     ("number_of_proteins", df.shape[0]),
    #     ("number_of_unique_hits", len(unique_identifiers)),
    #     ("ids", unique_identifiers),
    #     ("names", np.nan),
    #     ],
    # )
    # data = pd.Series(data, name=id_protein_cluster)
    # data.index = data.index.map(lambda x: (label, x))
    data = [df.shape[0], len(unique_identifiers), list(unique_identifiers)]

    return data

def compile_hmmsearch(df, id_protein_cluster, label):
    df = df.dropna(how="all", axis=0).query("number_of_hits > 0")
    unique_identifiers = flatten(df["ids"], into=list, unique=True)
    unique_names = flatten(df["names"], unique=True)
    
    # data = OrderedDict(
    #     [
    #     ("number_of_proteins", df.shape[0]),
    #     ("number_of_unique_hits", len(unique_identifiers)),
    #     ("ids", unique_identifiers),
    #     ("names", unique_names),
    #     ],
    # )
    # data = pd.Series(data, name=id_protein_cluster)
    # data.index = data.index.map(lambda x: (label, x))
    data = [df.shape[0], len(unique_identifiers), unique_identifiers, unique_names]
    return data


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i  <identifier_mapping> -d <uniref.blast6> -m <mibig.blast6> -b <vfdb.blast6> --hmmsearch_pfam <pfam.hmmsearch.tblout> --hmmsearch_amr <amr.hmmsearch.tblout> --hmmsearch_antifam <antifam.hmmsearch.tblout> -a <gene_models.faa> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_required = parser.add_argument_group('Required I/O arguments')
    parser_required.add_argument("-o","--output_directory", type=str, default="annotation_output",  help = "Output directory for annotations [Default: annotation_output]")
    parser_required.add_argument("-u","--diamond_uniref", type=str, required=True,  help = "path/to/diamond_uniref.blast6 in blast6 format (No header, cannot be empty) with the following fields: \nqseqid sseqid stitle pident length mismatch qlen qstart qend slen sstart send evalue bitscore qcovhsp scovhsp")
    parser_required.add_argument("-m","--diamond_mibig", type=str, required=True,  help = "path/to/diamond_mibig.blast6 in blast6 format (No header) with the following fields: \nqseqid sseqid stitle pident length mismatch qlen qstart qend slen sstart send evalue bitscore qcovhsp scovhsp")
    parser_required.add_argument("-b","--diamond_vfdb", type=str, required=True,  help = "path/to/diamond_vfdb.blast6 in blast6 format (No header) with the following fields: \nqseqid sseqid stitle pident length mismatch qlen qstart qend slen sstart send evalue bitscore qcovhsp scovhsp")
    parser_required.add_argument("-c","--diamond_cazy", type=str, required=True,  help = "path/to/diamond_cazy.blast6 in blast6 format (No header) with the following fields: \nqseqid sseqid stitle pident length mismatch qlen qstart qend slen sstart send evalue bitscore qcovhsp scovhsp")
    parser_required.add_argument("-k","--kofam", type=str, required=True,  help = "path/to/kofamscan.tblout hmmsearch results in tblout format")
    parser_required.add_argument("-p", "--hmmsearch_pfam", type=str, required=True,  help = "path/to/hmmsearch.tblout hmmsearch results in tblout format")
    parser_required.add_argument("-n", "--hmmsearch_amr", type=str, required=True,  help = "path/to/hmmsearch.tblout hmmsearch results in tblout format")
    parser_required.add_argument("-a", "--hmmsearch_antifam", type=str, required=True,  help = "path/to/hmmsearch.tblout hmmsearch results in tblout format")

    parser_optional = parser.add_argument_group('Optional arguments')
    parser_optional.add_argument("-i","--identifier_mapping", type=str, required=False,  help = "Tab-seperated value table (identifier_mapping.proteins.tsv created by cluster.py)  Requirements: 1) contain headers; 2) first column must be protein identifiers; and 3) contains these columns to the right in any order. Format: [id_protein]<tab>[organism_type]<tab>[id_genome]<tab>[sample_of_origin]<tab>[id_genome_cluster]<tab>[id_protein_cluster] with headers [Optional]")
    # parser_optional.add_argument("--genome_cluster_column_label", type=str, default="id_genome_cluster", help = "--genome_cluster_column_label must be in --identifier_mapping header [Default: id_genome_cluster]")
    # parser_optional.add_argument("--protein_cluster_column_label", type=str,  help = "--protein_cluster_column_label must be in --identifier_mapping header")
    parser_optional.add_argument("-j", "--composite_name_joiner", type=str, required=False,  default=";", help = "Composite label separator [Default: ; ]")
    # parser_optional.add_argument("--no_merge_composite_labels", action="store_true", help = "Do not merge composite labels")

    parser_optional.add_argument("-f","--fasta", type=str, required=False,  help = "path/to/gene_models.faa|ffn of ORFs [Optional]")
    # parser_optional.add_argument("-t","--threshold", default=0.5, type=float, help = "taxopy fraction of classifications for consensus for weighted majority vote [Default: 0.5]")
    # parser_optional.add_argument("-V", "--veba_database", type=str, required=False, help = "path/to/VEBA_DATABASE")
    parser_optional.add_argument("--hmmsearch_region", type=str, default="best",  help = "{best,full} Best domain or full sequence [Default: best]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Assertions
    assert_acceptable_arguments(opts.hmmsearch_region, {"best","full"})
    opts.hmmsearch_region = {"best":"best_domain", "full":"full_sequence"}[opts.hmmsearch_region]

    # Make output directories
    os.makedirs(opts.output_directory, exist_ok=True)
    print(" * Output directory:", opts.output_directory, file=sys.stderr)

    # Protein set
    proteins = set()

    # Read Fasta File
    if opts.fasta:
        with get_file_object(opts.fasta, "read", verbose=False) as f:
            for line in pv(f.readlines(), " * Getting protein identifiers from fasta file: {}".format(opts.fasta)):
                line = line.strip()
                if line.startswith(">"):
                    header = line[1:]
                    id = header.split(" ")[0]
                    proteins.add(id)
                    
    # Read identifier mapping table
    if opts.identifier_mapping:
        print(" * Reading identifier mapping table: {}".format(opts.identifier_mapping), file=sys.stderr)
        df_identifier_mapping = pd.read_csv(opts.identifier_mapping, sep="\t", index_col=0)
        fields_required = set(["organism_type", "id_genome", "sample_of_origin", "id_genome_cluster", "id_protein_cluster"])
        columns = set(df_identifier_mapping.columns)
        assert fields_required <= columns, "--identifier_mapping is missing the following columns: {}".format(fields_required - columns)
        df_identifier_mapping.columns = df_identifier_mapping.columns.map(lambda x: ("Identifiers", x))

        # df_identifier_mapping.columns = ["id_protein", "id_contig", "id_genome"]
        # df_identifier_mapping = df_identifier_mapping.set_index("id_protein")
        proteins = proteins | set(df_identifier_mapping.index)

    print(" * Reading Diamond table [UniRef]: {}".format(opts.diamond_uniref), file=sys.stderr)
    columns = ["qseqid","sseqid", "stitle", "pident","length","mismatch","qlen","qstart","qend", "slen", "sstart","send", "evalue","bitscore","qcovhsp","scovhsp"]
    # try:
    df_diamond_uniref = pd.read_csv(opts.diamond_uniref, sep="\t", index_col=None, header=None)
    # except pd.errors.EmptyDataError:
    #     df_diamond_uniref = pd.DataFrame(columns=columns)
    
    df_diamond_uniref.columns = columns
    df_diamond_uniref = df_diamond_uniref.set_index("qseqid")

    # Remove ID from stitle
    df_diamond_uniref["stitle"] = df_diamond_uniref["stitle"].map(clean_stitle)

    stitle_data = list()
    for id_protein, stitle in pv(df_diamond_uniref["stitle"].items(), "Parsing UniRef descriptions"):
        fields = parse_uniref(stitle)
        stitle_data.append(fields)
    df_uniref_stitle = pd.DataFrame(stitle_data, index=df_diamond_uniref.index, columns=["product", "cluster_size", "id_taxon", "id_representative"])

    # Concatenate dataframes and remove redundancy
    df_diamond_uniref = pd.concat([df_diamond_uniref, df_uniref_stitle], axis=1)
    df_diamond_uniref = df_diamond_uniref.drop(["stitle"], axis=1)

    proteins = proteins | set(df_diamond_uniref.index)

    print(" * Reading Diamond table [MIBiG]: {}".format(opts.diamond_mibig), file=sys.stderr)
    columns = ["qseqid","sseqid", "stitle", "pident","length","mismatch","qlen","qstart","qend", "slen", "sstart","send", "evalue","bitscore","qcovhsp","scovhsp"]
    try:
        df_diamond_mibig = pd.read_csv(opts.diamond_mibig, sep="\t", index_col=None, header=None)
    except pd.errors.EmptyDataError:
        df_diamond_mibig = pd.DataFrame(columns=columns)
    df_diamond_mibig.columns = columns
    df_diamond_mibig = df_diamond_mibig.set_index("qseqid")
    df_diamond_mibig = df_diamond_mibig.drop(["stitle"], axis=1)
    proteins = proteins | set(df_diamond_mibig.index)

    print(" * Reading Diamond table [VFDB]: {}".format(opts.diamond_mibig), file=sys.stderr)
    columns = ["qseqid","sseqid", "stitle", "pident","length","mismatch","qlen","qstart","qend", "slen", "sstart","send", "evalue","bitscore","qcovhsp","scovhsp"]
    try:
        df_diamond_vfdb = pd.read_csv(opts.diamond_vfdb, sep="\t", index_col=None, header=None)
    except pd.errors.EmptyDataError:
        df_diamond_vfdb = pd.DataFrame(columns=columns)
    df_diamond_vfdb.columns = columns
    df_diamond_vfdb = df_diamond_vfdb.set_index("qseqid")
    df_diamond_vfdb = df_diamond_vfdb.drop(["stitle"], axis=1)
    proteins = proteins | set(df_diamond_vfdb.index)

    print(" * Reading Diamond table [CAZy]: {}".format(opts.diamond_cazy), file=sys.stderr)
    columns = ["qseqid","sseqid", "stitle", "pident","length","mismatch","qlen","qstart","qend", "slen", "sstart","send", "evalue","bitscore","qcovhsp","scovhsp"]
    try:
        df_diamond_cazy = pd.read_csv(opts.diamond_cazy, sep="\t", index_col=None, header=None)
    except pd.errors.EmptyDataError:
        df_diamond_cazy = pd.DataFrame(columns=columns)
    df_diamond_cazy.columns = columns
    df_diamond_cazy = df_diamond_cazy.set_index("qseqid")
    df_diamond_cazy = df_diamond_cazy.drop(["stitle"], axis=1)
    proteins = proteins | set(df_diamond_cazy.index)

    print(" * Reading HMMSearch table [Pfam]: {}".format(opts.hmmsearch_pfam), file=sys.stderr)
    df_hmmsearch_pfam = read_hmmer(opts.hmmsearch_pfam, program="hmmsearch", format="tblout", verbose=False)
    if not df_hmmsearch_pfam.empty:
        df_hmmsearch_pfam = df_hmmsearch_pfam[["identifier", opts.hmmsearch_region]].droplevel(0, axis=1)
        proteins = proteins | set(df_hmmsearch_pfam["target_name"])
    else:
        print(" *!* HMMSearch table [Pfam] is empty", file=sys.stderr)

    print(" * Reading HMMSearch table [AMR]: {}".format(opts.hmmsearch_amr), file=sys.stderr)
    df_hmmsearch_amr = read_hmmer(opts.hmmsearch_amr, program="hmmsearch", format="tblout", verbose=False)
    if not df_hmmsearch_amr.empty:
        df_hmmsearch_amr = df_hmmsearch_amr[["identifier", opts.hmmsearch_region]].droplevel(0, axis=1)
        proteins = proteins | set(df_hmmsearch_amr["target_name"])
    else:
        print(" *!* HMMSearch table [AMR] is empty", file=sys.stderr)

    print(" * Reading HMMSearch table [AntiFam]: {}".format(opts.hmmsearch_antifam), file=sys.stderr)
    df_hmmsearch_antifam = read_hmmer(opts.hmmsearch_antifam, program="hmmsearch", format="tblout", verbose=False)
    if not df_hmmsearch_antifam.empty:
        df_hmmsearch_antifam = df_hmmsearch_antifam[["identifier", opts.hmmsearch_region]].droplevel(0, axis=1)
        proteins = proteins | set(df_hmmsearch_antifam[ "target_name"])
    else:
        print(" *!* HMMSearch table [AntiFam] is empty", file=sys.stderr)

    print(" * Reading KOFAMSCAN table [KEGG]: {}".format(opts.kofam), file=sys.stderr)
    columns =  ["significance", "id_protein", "id_hmm", "score","score2", "evalue", "name_hmm"]
    try:
        df_kofamscan = pd.read_csv(opts.kofam, sep="\t", index_col=None, header=None)
        df_kofamscan.columns = columns
    except pd.errors.EmptyDataError:
        df_kofamscan = pd.DataFrame(columns=columns)


    if not df_kofamscan.empty:
        proteins = proteins | set(df_kofamscan.iloc[:,1])
    else:
        print(" *!* KOFAMSCAN table [KEGG] is empty", file=sys.stderr)


    # All proteins
    proteins = sorted(proteins)
    print(" * Total proteins: {}".format(len(proteins)), file=sys.stderr)

    # Reindex Diamond
    df_diamond_uniref = df_diamond_uniref.reindex(proteins)

    # Pfam
    if not df_hmmsearch_pfam.empty:
        protein_to_hmmid = defaultdict(list)
        protein_to_hmmname = defaultdict(list)
        protein_to_evalues = defaultdict(list)
        protein_to_scores = defaultdict(list)

        for i, row in pv(df_hmmsearch_pfam.iterrows(), "Formatting HMMSearch [Pfam]", total=df_hmmsearch_pfam.shape[0], unit=" search hits"):
            id_protein, name_hmm, id_hmm, evalue, score = row[["target_name", "query_name", "query_accession", "e-value", "score"]]
            protein_to_hmmid[id_protein].append(id_hmm)
            protein_to_hmmname[id_protein].append(name_hmm)
            protein_to_evalues[id_protein].append(evalue)
            protein_to_scores[id_protein].append(score)

        df_hmms_pfam = pd.DataFrame(
            OrderedDict([
                ("ids", protein_to_hmmid),
                ("names", protein_to_hmmname),
                ("evalues", protein_to_evalues),
                ("scores", protein_to_scores),
            ]
        ))
        df_hmms_pfam = df_hmms_pfam.reindex(proteins)

    else:
        df_hmms_pfam = pd.DataFrame(index=proteins, columns=["hit", "ids", "names", "evalues", "scores"])

    df_hmms_pfam = df_hmms_pfam.applymap(lambda x: x if isinstance(x,list) else [])
    df_hmms_pfam.insert(loc=0, column="number_of_hits", value=df_hmms_pfam["ids"].map(len))
    df_hmms_pfam.index.name = "id_protein"

    # AMR
    if not df_hmmsearch_amr.empty:
        protein_to_hmmid = defaultdict(list)
        protein_to_hmmname = defaultdict(list)
        protein_to_evalues = defaultdict(list)
        protein_to_scores = defaultdict(list)

        for i, row in pv(df_hmmsearch_amr.iterrows(), "Formatting HMMSearch [AMR]", total=df_hmmsearch_amr.shape[0], unit=" search hits"):
            id_protein, name_hmm, id_hmm, evalue, score = row[["target_name", "query_name", "query_accession", "e-value", "score"]]
            protein_to_hmmid[id_protein].append(id_hmm)
            protein_to_hmmname[id_protein].append(name_hmm)
            protein_to_evalues[id_protein].append(evalue)
            protein_to_scores[id_protein].append(score)

        df_hmms_amr = pd.DataFrame(
            OrderedDict([
                ("ids", protein_to_hmmid),
                ("names", protein_to_hmmname),
                ("evalues", protein_to_evalues),
                ("scores", protein_to_scores),
            ]
        ))
        df_hmms_amr = df_hmms_amr.reindex(proteins)
    else:
        df_hmms_amr = pd.DataFrame(index=proteins, columns=["hit", "ids", "names", "evalues", "scores"])
    df_hmms_amr = df_hmms_amr.applymap(lambda x: x if isinstance(x,list) else [])
    df_hmms_amr.index.name = "id_protein"
    df_hmms_amr.insert(loc=0, column="number_of_hits", value=df_hmms_amr["ids"].map(len))


    # AntiFam
    if not df_hmmsearch_antifam.empty:
        protein_to_hmmid = defaultdict(list)
        protein_to_hmmname = defaultdict(list)
        protein_to_evalues = defaultdict(list)
        protein_to_scores = defaultdict(list)

        for i, row in pv(df_hmmsearch_antifam.iterrows(), "Formatting HMMSearch [AntiFam]", total=df_hmmsearch_antifam.shape[0], unit=" search hits"):
            id_protein, name_hmm, id_hmm, evalue, score = row[["target_name", "query_name", "query_accession", "e-value", "score"]]
            protein_to_hmmid[id_protein].append(id_hmm)
            protein_to_hmmname[id_protein].append(name_hmm)
            protein_to_evalues[id_protein].append(evalue)
            protein_to_scores[id_protein].append(score)

        df_hmms_antifam = pd.DataFrame(
            OrderedDict([
                ("ids", protein_to_hmmid),
                ("names", protein_to_hmmname),
                ("evalues", protein_to_evalues),
                ("scores", protein_to_scores),
            ]
        ))
        df_hmms_antifam.insert(loc=0, column="hit", value=True)
        df_hmms_antifam = df_hmms_antifam.reindex(proteins)
    else:
        df_hmms_antifam = pd.DataFrame(index=proteins, columns=["hit", "ids", "names", "evalues", "scores"])
    df_hmms_antifam["hit"] = df_hmms_antifam["hit"].fillna(False)
    df_hmms_antifam = df_hmms_antifam.applymap(lambda x: x if isinstance(x,list) else [])
    df_hmms_antifam.index.name = "id_protein"
    df_hmms_antifam.insert(loc=0, column="number_of_hits", value=df_hmms_antifam["ids"].map(len))

    # KOFAMSCAN
    if not df_kofamscan.empty:
        protein_to_hmmid = defaultdict(list)
        protein_to_hmmname = defaultdict(list)
        protein_to_evalues = defaultdict(list)
        protein_to_scores = defaultdict(list)

        for i, row in pv(df_kofamscan.iterrows(), "Formatting KOFAMSCAN", total=df_kofamscan.shape[0], unit=" search hits"):
            significance, id_protein, id_hmm, score, evalue, name_hmm  = row.iloc[[0, 1,2, 4, 5, 6]]
            if significance == "*":
                protein_to_hmmid[id_protein].append(id_hmm)
                protein_to_hmmname[id_protein].append(name_hmm)
                protein_to_evalues[id_protein].append(evalue)
                protein_to_scores[id_protein].append(score)
                
        df_hmms_kofam = pd.DataFrame(
            OrderedDict([
                ("ids", protein_to_hmmid),
                ("names", protein_to_hmmname),
                ("evalues", protein_to_evalues),
                ("scores", protein_to_scores),
            ]
        ))
        df_hmms_kofam = df_hmms_kofam.reindex(proteins)
    else:
        df_hmms_kofam = pd.DataFrame(index=proteins, columns=["hit", "ids", "names", "evalues", "scores"])
    df_hmms_kofam = df_hmms_kofam.applymap(lambda x: x if isinstance(x,list) else [])
    df_hmms_kofam.index.name = "id_protein"
    df_hmms_kofam.insert(loc=0, column="number_of_hits", value=df_hmms_kofam["ids"].map(len))

    # Output table
    dataframes = list()
    for (name, df) in zip([ "UniRef", "MIBiG", "VFDB", "CAZy", "Pfam", "NCBIfam-AMR", "AntiFam", "KOFAM"],[df_diamond_uniref, df_diamond_mibig, df_diamond_vfdb, df_diamond_cazy, df_hmms_pfam, df_hmms_amr, df_hmms_antifam, df_hmms_kofam]):
        df = df.copy()
        df.columns = df.columns.map(lambda x: (name,x))
        dataframes.append(df)
    df_annotations = pd.concat(dataframes, axis=1)
    df_annotations.index.name = "id_protein"

    # Composite label
    protein_to_labels = dict()
    for id_protein, row in pv(df_annotations.iterrows(), total=df_annotations.shape[0]):
        labels = list()

        name = row["UniRef"]["product"]
        if name not in labels:
            labels.append(name)

        kofam_names = row["KOFAM"]["names"]
        for name in kofam_names:
            if name not in labels:
                labels.append(name)
        
        pfam_names = row["Pfam"]["names"]
        for name in pfam_names:
            if name not in labels:
                labels.append(name)

        composite_name = sorted(filter(lambda x: isinstance(x,str), labels))
        protein_to_labels[id_protein] =  opts.composite_name_joiner.join(composite_name)
    protein_to_labels = pd.Series(protein_to_labels)

    # if not opts.no_merge_composite_labels:
    # protein_to_labels = protein_to_labels.map(lambda x: opts.composite_name_joiner.join(x))
    
    df_annotations.insert(loc=0, column=("Consensus", "composite_name"), value=protein_to_labels)

    if opts.identifier_mapping:
        # df = df_identifier_mapping.copy()
        # df.columns = df.columns.map(lambda x: ("Identifiers", x))
        df_annotations = pd.concat([df_identifier_mapping, df_annotations], axis=1)

    df_annotations.to_csv(os.path.join(opts.output_directory, "annotations.proteins.tsv.gz"), sep="\t")

    if opts.identifier_mapping:
        with gzip.open(os.path.join(opts.output_directory, "annotations.protein_clusters.tsv.gz"), "wt") as f:
            print("\t", 
                  *["Identifiers"]*2,
                  *["Consensus"]*1,

                  *["UniRef"]*4,
                  *["MIBiG"]*3,
                  *["VFDB"]*3,
                  *["CAZy"]*3,
                  *["Pfam"]*4,
                  *["NCBIfam-AMR"]*4,
                  *["KOFAM"]*4,
                  *["AntiFam"]*4,
                  sep="\t", file=f)

            print(
                "id_protein_cluster", 
                *["id_genome_cluster", "organsim_type"], #, "genomes", "samples_of_origin"], # Identifiers
                *["composite_name"], # Consensus
                *["number_of_proteins", "number_of_unique_hits", "ids","names"], # UniRef
                *["number_of_proteins", "number_of_unique_hits", "ids"], # MIBiG
                *["number_of_proteins", "number_of_unique_hits", "ids"], # VFDB
                *["number_of_proteins", "number_of_unique_hits", "ids"], # CAZy
                *["number_of_proteins", "number_of_unique_hits", "ids","names"], # Pfam
                *["number_of_proteins", "number_of_unique_hits", "ids","names"], # NCBIfam-AMR
                *["number_of_proteins", "number_of_unique_hits", "ids","names"], # KOFAM
                *["number_of_proteins", "number_of_unique_hits", "ids","names"], # AntiFam
                sep="\t",
                file=f,
                )
            # Protein clusters
            protein_to_proteincluster = df_annotations[("Identifiers", "id_protein_cluster")]
            protein_cluster_annotations = list()
            for id_protein_cluster, df in pv(df_annotations.groupby(protein_to_proteincluster), description="Compiling consensus annotations for protein clusters", total=protein_to_proteincluster.nunique(), unit=" Protein Clusters"):
                # Identifiers
                data_identifiers = compile_identifiers(df["Identifiers"], id_protein_cluster)
                
                # UniRef
                data_uniref = compile_uniref(df["UniRef"], id_protein_cluster)
                
                # MIBiG
                data_mibig = compile_nonuniref_diamond(df["MIBiG"], id_protein_cluster, "MIBiG")

                # VFDB
                data_vfdb = compile_nonuniref_diamond(df["VFDB"], id_protein_cluster, "VFDB")

                # CAZy
                data_cazy = compile_nonuniref_diamond(df["CAZy"], id_protein_cluster, "CAZy")

                # Pfam
                data_pfam = compile_hmmsearch(df["Pfam"], id_protein_cluster, "Pfam")

                # NCBIfam-AMR
                data_amr = compile_hmmsearch(df["NCBIfam-AMR"], id_protein_cluster, "NCBIfam-AMR")

                # KOFAM
                data_kofam = compile_hmmsearch(df["KOFAM"], id_protein_cluster, "KOFAM")
                
                # AntiFam
                data_antifam = compile_hmmsearch(df["AntiFam"], id_protein_cluster, "AntiFam")

                # Composite name
                composite_name = list()
                composite_name += list(data_uniref[-1]) 
                composite_name += list(data_kofam[-1]) 
                composite_name += list(data_pfam[-1]) 
                composite_name = list(filter(lambda x: isinstance(x, str), composite_name))
                if len(composite_name) > 0:
                    composite_name = opts.composite_name_joiner.join(composite_name)
                else:
                    composite_name = np.nan

                print(
                    id_protein_cluster, 
                    *data_identifiers, 
                    composite_name, 
                    *data_uniref,
                    *data_mibig,
                    *data_vfdb,
                    *data_cazy,
                    *data_pfam,
                    *data_amr,
                    *data_kofam,
                    *data_antifam,
                    sep="\t", 
                    file=f,
                    )
                
                # data_consensus = pd.Series(composite_name, index=[("Consensus", "composite_name")])
                # # Concatenate
                # data_concatenated = pd.concat([
                #     data_identifiers, 
                #     data_consensus,
                #     data_uniref,
                #     data_mibig,
                #     data_vfdb,
                #     data_cazy,
                #     data_pfam,
                #     data_amr,
                #     data_kofam,
                #     data_antifam,
                # ])
                # data_concatenated.name = id_protein_cluster
                # protein_cluster_annotations.append(data_concatenated)
            # df_annotations_proteinclusters = pd.DataFrame(protein_cluster_annotations)
            # df_annotations_proteinclusters.to_csv(os.path.join(opts.output_directory, "annotations.protein_clusters.tsv.gz"), sep="\t")





if __name__ == "__main__":
    main()
    
                
