#!/usr/bin/env python
import sys, os, argparse, re, gzip
from collections import defaultdict, OrderedDict
import pandas as pd
import numpy as np
from soothsayer_utils import read_hmmer, pv, get_file_object, assert_acceptable_arguments, format_header, flatten

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2025.5.26"

DIAMOND_HEADER_FIELDS = "qseqid sseqid stitle pident evalue bitscore qcovhsp scovhsp"
DIAMOND_COLUMNS = list(filter(bool, DIAMOND_HEADER_FIELDS.split(" ")))
    
def clean_stitle(stitle):
    sseqid, *stitle = stitle.split(" ")
    stitle = " ".join(stitle)
    return stitle 

def parse_uniref(stitle):

    # Define the regular expression pattern
    pattern = r'^(.*?)\s+n=(\d+)\s+Tax=.*?\s+TaxID=(\d+)\s+RepID=(.*)$'

    # Use re.match to search for the pattern in the string
    match = re.match(pattern, stitle)
    if match:
        # Extract the four fields from the match object
        field1 = match.group(1)
        field1 = np.nan if field1 is None else field1
        field2 = match.group(2)
        field2 = np.nan if field2 is None else field2
        field3 = match.group(3)
        field3 = np.nan if field3 is None else field3
        field4 = match.group(4)
        field4 = np.nan if field4 is None else field4
    else:
        field1 = np.nan
        field2 = np.nan
        field3 = np.nan
        field4 = np.nan
        
    # Print the four fields
    return [field1, field2, field3, field4]

def compile_identifiers(df):
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
    data = [genome_clusters, organism_types]#, genomes, samples]
    return data

def compile_uniref(df):
    df = df.dropna(how="all", axis=0)
    unique_identifiers = list(df["sseqid"].unique())
    data = [df.shape[0], len(unique_identifiers), unique_identifiers, list(df["product"].unique())]
    return data

def compile_nonuniref_diamond(df):
    df = df.dropna(how="all", axis=0)
    unique_identifiers = set(df["sseqid"].unique())
    data = [df.shape[0], len(unique_identifiers), list(unique_identifiers)]

    return data

def compile_pfam(df):
    df = df.dropna(how="all", axis=0).query("number_of_hits > 0")
    unique_identifiers = flatten(df["ids"], into=list, unique=True)   
    unique_names = flatten(df["names"], unique=True)
    data = [df.shape[0], len(unique_identifiers), unique_identifiers, unique_names]
    return data

def compile_nonpfam_pyhmmsearch(df):
    df = df.dropna(how="all", axis=0).query("number_of_hits > 0")
    unique_identifiers = flatten(df["ids"], into=list, unique=True)    
    data = [df.shape[0], len(unique_identifiers), unique_identifiers]
    return data

def compile_pykofamsearch(df):
    df = df.dropna(how="all", axis=0).query("number_of_hits > 0")
    unique_identifiers = flatten(df["ids"], into=list, unique=True)
    unique_names = flatten(df["names"], unique=True)
    data = [df.shape[0], len(unique_identifiers), unique_identifiers, unique_names]
    return data

def compile_enzymes(df):
    df = df.dropna(how="all", axis=0).query("number_of_hits > 0")
    unique_identifiers = flatten(df["ids"], into=list, unique=True)
    data = [df.shape[0], len(unique_identifiers), unique_identifiers]
    return data


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i  <identifier_mapping> -d <uniref.blast6> -m <mibig.blast6> -b <vfdb.blast6> --pyhmmsearch_pfam <pfam.tsv.gz> --pyhmmsearch_amr <amr.tsv.gz> --pyhmmsearch_antifam <antifam.tsv.gz> --pfam_clans <Pfam-A.clans.tsv.gz> -a <gene_models.faa> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_required = parser.add_argument_group('Required I/O arguments')
    parser_required.add_argument("-o","--output_directory", type=str, default="annotation_output",  help = "Output directory for annotations [Default: annotation_output]")
    parser_required.add_argument("-u","--diamond_uniref", type=str, required=True,  help = f"path/to/diamond_uniref.blast6 in blast6 format (No header, cannot be empty) with the following fields: \n{DIAMOND_COLUMNS}")
    parser_required.add_argument("-m","--diamond_mibig", type=str, required=True,  help = f"path/to/diamond_mibig.blast6 in blast6 format (No header) with the following fields: \n{DIAMOND_COLUMNS}")
    parser_required.add_argument("-b","--diamond_vfdb", type=str, required=True,  help = f"path/to/diamond_vfdb.blast6 in blast6 format (No header) with the following fields: \n{DIAMOND_COLUMNS}")
    parser_required.add_argument("-c","--diamond_cazy", type=str, required=True,  help = f"path/to/diamond_cazy.blast6 in blast6 format (No header) with the following fields: \n{DIAMOND_COLUMNS}")
    parser_required.add_argument("-k","--pykofamsearch", type=str, required=True,  help = "path/to/pykofamsearch_output.tsv with header")
    parser_required.add_argument("-p", "--pyhmmsearch_pfam", type=str, required=True,  help = "path/to/pyhmmsearch_output.tsv with header")
    parser_required.add_argument("-n", "--pyhmmsearch_amr", type=str, required=True,  help = "path/to/pyhmmsearch_output.tsv with header")
    parser_required.add_argument("-a", "--pyhmmsearch_antifam", type=str, required=True,  help = "path/to/pyhmmsearch_output.tsv with header")
    parser_required.add_argument("-q", "--pfam_clans", type=str, required=True,  help = "path/to/Pfam-A.clans.tsv.gz. No header")

    parser_optional = parser.add_argument_group('Optional arguments')
    parser_optional.add_argument("-i","--identifier_mapping", type=str, required=False,  help = "Tab-seperated value table (identifier_mapping.proteins.tsv created by cluster.py)  Requirements: 1) contain headers; 2) first column must be protein identifiers; and 3) contains these columns to the right in any order. Format: [id_protein]<tab>[organism_type]<tab>[id_genome]<tab>[sample_of_origin]<tab>[id_genome_cluster]<tab>[id_protein_cluster] with headers [Optional]")
    # parser_optional.add_argument("--genome_cluster_column_label", type=str, default="id_genome_cluster", help = "--genome_cluster_column_label must be in --identifier_mapping header [Default: id_genome_cluster]")
    # parser_optional.add_argument("--protein_cluster_column_label", type=str,  help = "--protein_cluster_column_label must be in --identifier_mapping header")
    parser_optional.add_argument("-j", "--composite_name_joiner", type=str, required=False,  default=";", help = "Composite label separator [Default: ; ]")

    parser_optional.add_argument("-f","--fasta", type=str, required=False,  help = "path/to/gene_models.faa|ffn of ORFs [Optional]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

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

    # =======
    # Diamond
    # =======
    print(" * Reading Diamond table [UniRef]: {}".format(opts.diamond_uniref), file=sys.stderr)

    
    # try:
    df_diamond_uniref = pd.read_csv(opts.diamond_uniref, sep="\t", index_col=0)
    # except pd.errors.EmptyDataError:
    #     df_diamond_uniref = pd.DataFrame(columns=columns)
    
    # df_diamond_uniref.columns = DIAMOND_COLUMNS
    # df_diamond_uniref = df_diamond_uniref.set_index("qseqid")

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
    try:
        df_diamond_mibig = pd.read_csv(opts.diamond_mibig, sep="\t", index_col=0)
    except pd.errors.EmptyDataError:
        df_diamond_mibig = pd.DataFrame(columns=DIAMOND_COLUMNS)
    # df_diamond_mibig = df_diamond_mibig.set_index("qseqid")
    df_diamond_mibig = df_diamond_mibig.drop(["stitle"], axis=1)
    proteins = proteins | set(df_diamond_mibig.index)

    print(" * Reading Diamond table [VFDB]: {}".format(opts.diamond_mibig), file=sys.stderr)
    try:
        df_diamond_vfdb = pd.read_csv(opts.diamond_vfdb, sep="\t", index_col=0)
    except pd.errors.EmptyDataError:
        df_diamond_vfdb = pd.DataFrame(columns=DIAMOND_COLUMNS)
    # df_diamond_vfdb = df_diamond_vfdb.set_index("qseqid")
    df_diamond_vfdb = df_diamond_vfdb.drop(["stitle"], axis=1)
    proteins = proteins | set(df_diamond_vfdb.index)

    print(" * Reading Diamond table [CAZy]: {}".format(opts.diamond_cazy), file=sys.stderr)
    try:
        df_diamond_cazy = pd.read_csv(opts.diamond_cazy, sep="\t", index_col=0)
    except pd.errors.EmptyDataError:
        df_diamond_cazy = pd.DataFrame(columns=DIAMOND_COLUMNS)
    # df_diamond_cazy = df_diamond_cazy.set_index("qseqid")
    df_diamond_cazy = df_diamond_cazy.drop(["stitle"], axis=1)
    proteins = proteins | set(df_diamond_cazy.index)

    # ==========
    # Pfam Clans
    # ==========
    pfam_to_name = pd.read_csv(opts.pfam_clans, sep="\t", index_col=0, header=None).iloc[:,-1].to_dict()
    
    # ===========
    # PyHMMSearch
    # ===========
    print(" * Reading PyHMMSearch table [Pfam]: {}".format(opts.pyhmmsearch_pfam), file=sys.stderr)
    df_hmms_pfam = pd.read_csv(opts.pyhmmsearch_pfam, sep="\t", index_col=0)
    if not df_hmms_pfam.empty:
        for field in ["ids", "evalues", "scores"]:
            df_hmms_pfam[field] = df_hmms_pfam[field].map(eval)
        protein_to_pfamnames = defaultdict(list)
        for id_protein, hmms in df_hmms_pfam["ids"].items():
            for id_hmm in hmms:
                if "." in id_hmm:
                    id_hmm = id_hmm.split(".")[0]
                if id_hmm in pfam_to_name:
                    name = pfam_to_name[id_hmm]
                    protein_to_pfamnames[id_protein].append(name)
        protein_to_pfamnames = pd.Series(protein_to_pfamnames)
        df_hmms_pfam.insert(2, "names", protein_to_pfamnames)
        proteins = proteins | set(df_hmms_pfam.index)
    else:
        print(" *!* PyHMMSearch table [Pfam] is empty", file=sys.stderr)

    print(" * Reading PyHMMSearch table [AMR]: {}".format(opts.pyhmmsearch_amr), file=sys.stderr)
    df_hmms_amr = pd.read_csv(opts.pyhmmsearch_amr, sep="\t", index_col=0)
    if not df_hmms_amr.empty:
        for field in ["ids", "evalues", "scores"]:
            df_hmms_amr[field] = df_hmms_amr[field].map(eval)
        proteins = proteins | set(df_hmms_amr.index)
    else:
        print(" *!* PyHMMSearch table [AMR] is empty", file=sys.stderr)

    print(" * Reading PyHMMSearch table [AntiFam]: {}".format(opts.pyhmmsearch_antifam), file=sys.stderr)
    df_hmms_antifam = pd.read_csv(opts.pyhmmsearch_antifam, sep="\t", index_col=0)
    if not df_hmms_antifam.empty:
        for field in ["ids", "evalues", "scores"]:
            df_hmms_antifam[field] = df_hmms_antifam[field].map(eval)
        proteins = proteins | set(df_hmms_antifam.index)
    else:
        print(" *!* PyHMMSearch table [AntiFam] is empty", file=sys.stderr)

    # ===========
    # PyKOfamSearch
    # ===========
    print(" * Reading PyKofamSearch table [KEGG]: {}".format(opts.pykofamsearch), file=sys.stderr)
    df_hmms_kofam = pd.read_csv(opts.pykofamsearch, sep="\t", index_col=0)
    # Process enzymes
    df_enzymes = df_hmms_kofam.pop("enzyme_commissions").to_frame("ids")

    if not df_hmms_kofam.empty:
        for field in ["ids", "names", "evalues", "scores"]:
            df_hmms_kofam[field] = df_hmms_kofam[field].map(eval)
        proteins = proteins | set(df_hmms_kofam.index)

    else:
        print(" *!* PyKofamSearch table [KOfam] is empty", file=sys.stderr)
        
    if not df_enzymes.empty:
        for field in ["ids"]:
            df_enzymes[field] = df_enzymes[field].map(eval)
        # proteins = proteins | set(df_enzymes.index)

    else:
        print(" *!* PyKofamSearch table [Enzymes] is empty", file=sys.stderr)

    # All proteins
    proteins = sorted(proteins)
    print(" * Total proteins: {}".format(len(proteins)), file=sys.stderr)

    # Reindex Diamond
    df_diamond_uniref = df_diamond_uniref.reindex(proteins)

    # Reindex HMMSearch Results
    # Pfam
    df_hmms_pfam = df_hmms_pfam.reindex(proteins)
    for field in ["ids", "names", "evalues", "scores"]:
        df_hmms_pfam[field] = df_hmms_pfam[field].map(lambda x: x if isinstance(x,list) else [])
    df_hmms_pfam["number_of_hits"] = df_hmms_pfam["number_of_hits"].fillna(0).astype(int)
    df_hmms_pfam.index.name = "id_protein"

    # AMR
    df_hmms_amr = df_hmms_amr.reindex(proteins)
    for field in ["ids", "evalues", "scores"]:
        df_hmms_amr[field] = df_hmms_amr[field].map(lambda x: x if isinstance(x,list) else [])
    df_hmms_amr["number_of_hits"] = df_hmms_amr["number_of_hits"].fillna(0).astype(int)
    df_hmms_amr.index.name = "id_protein"
        
    # AntiFam
    df_hmms_antifam = df_hmms_antifam.reindex(proteins)
    for field in ["ids", "evalues", "scores"]:
        df_hmms_antifam[field] = df_hmms_antifam[field].map(lambda x: x if isinstance(x,list) else [])
    df_hmms_antifam["number_of_hits"] = df_hmms_antifam["number_of_hits"].fillna(0).astype(int)
    df_hmms_antifam.index.name = "id_protein"

    # KOfam
    df_hmms_kofam = df_hmms_kofam.reindex(proteins)
    for field in ["ids", "names", "evalues", "scores"]:
        df_hmms_kofam[field] = df_hmms_kofam[field].map(lambda x: x if isinstance(x,list) else [])
    df_hmms_kofam["number_of_hits"] = df_hmms_kofam["number_of_hits"].fillna(0).astype(int)
    df_hmms_kofam.index.name = "id_protein"
    
    # Enzymes
    df_enzymes = df_enzymes.reindex(proteins)
    for field in ["ids"]:
        df_enzymes[field] = df_enzymes[field].map(lambda x: x if isinstance(x,list) else [])
    df_enzymes["number_of_hits"] = df_enzymes["ids"].map(len).fillna(0).astype(int)
    df_enzymes.index.name = "id_protein"
       
    # Output table
    dataframes = list()
    for (name, df) in zip([ "UniRef", "MIBiG", "VFDB", "CAZy", "Pfam", "NCBIfam-AMR",  "KOfam", "Enzymes", "AntiFam"],[df_diamond_uniref, df_diamond_mibig, df_diamond_vfdb, df_diamond_cazy, df_hmms_pfam, df_hmms_amr, df_hmms_kofam, df_enzymes, df_hmms_antifam]):
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

        kofam_names = row["KOfam"]["names"]
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
            print(
                  # "\t", 
                  *["Identifiers"]*3,
                  *["Consensus"]*1,
                  *["UniRef"]*4,
                  *["MIBiG"]*3,
                  *["VFDB"]*3,
                  *["CAZy"]*3,
                  *["Pfam"]*4,
                  *["NCBIfam-AMR"]*3,
                  *["KOfam"]*4,
                  *["Enzymes"]*3,
                  *["AntiFam"]*3,
                  sep="\t", file=f)

            print(
                "id_protein_cluster", 
                *["id_genome_cluster", "organsim_type"], #, "genomes", "samples_of_origin"], # Identifiers
                *["composite_name"], # Consensus
                *["number_of_proteins", "number_of_unique_hits", "ids","names"], # UniRef
                *["number_of_proteins", "number_of_unique_hits", "ids"], # MIBiG
                *["number_of_proteins", "number_of_unique_hits", "ids"], # VFDB
                *["number_of_proteins", "number_of_unique_hits", "ids"], # CAZy
                *["number_of_proteins", "number_of_unique_hits", "ids", "names"], # Pfam
                *["number_of_proteins", "number_of_unique_hits", "ids"], # NCBIfam-AMR
                *["number_of_proteins", "number_of_unique_hits", "ids","names"], # KOfam
                *["number_of_proteins", "number_of_unique_hits", "ids"], # Enzymes
                *["number_of_proteins", "number_of_unique_hits", "ids"], # AntiFam
                sep="\t",
                file=f,
                )
            # Protein clusters
            protein_to_proteincluster = df_annotations[("Identifiers", "id_protein_cluster")]
            protein_cluster_annotations = list()
            for id_protein_cluster, df in pv(df_annotations.groupby(protein_to_proteincluster), description="Compiling consensus annotations for protein clusters", total=protein_to_proteincluster.nunique(), unit=" Protein Clusters"):
                # Identifiers
                data_identifiers = compile_identifiers(df["Identifiers"])
                
                # UniRef
                data_uniref = compile_uniref(df["UniRef"])
                
                # MIBiG
                data_mibig = compile_nonuniref_diamond(df["MIBiG"])

                # VFDB
                data_vfdb = compile_nonuniref_diamond(df["VFDB"])

                # CAZy
                data_cazy = compile_nonuniref_diamond(df["CAZy"])

                # Pfam
                data_pfam = compile_pfam(df["Pfam"])

                # NCBIfam-AMR
                data_amr = compile_nonpfam_pyhmmsearch(df["NCBIfam-AMR"])

                # KOfam
                data_kofam = compile_pykofamsearch(df["KOfam"])
                
                # Enzymes
                data_enzymes = compile_enzymes(df["Enzymes"])
                
                # AntiFam
                data_antifam = compile_nonpfam_pyhmmsearch(df["AntiFam"])

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
                    *data_enzymes,
                    *data_antifam,
                    sep="\t", 
                    file=f,
                    )
                
if __name__ == "__main__":
    main()
    
                
