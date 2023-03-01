#!/usr/bin/env python
import sys, os, argparse 
from collections import defaultdict, OrderedDict
import pandas as pd
import numpy as np
from soothsayer_utils import read_hmmer, pv, get_file_object, assert_acceptable_arguments, format_header

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.2.28"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i  <identifier_mapping> -d <nr.blast6> --hmmsearch_pfam <pfam.hmmsearch.tblout> --hmmsearch_amr <amr.hmmsearch.tblout> --hmmsearch_antifam <antifam.hmmsearch.tblout> -a <gene_models.faa> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_required = parser.add_argument_group('Required I/O arguments')
    parser_required.add_argument("-o","--output_directory", type=str, default="annotation_output",  help = "Output directory for annotations [Default: annotation_output]")
    parser_required.add_argument("-d","--diamond", type=str, required=True,  help = "path/to/diamond_nr.blast6 in blast6 format (No header) with the following fields: \nqseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle")
    parser_required.add_argument("-k","--kofam", type=str, required=True,  help = "path/to/kofamscan.tblout hmmsearch results in tblout format")
    parser_required.add_argument("-p", "--hmmsearch_pfam", type=str, required=True,  help = "path/to/hmmsearch.tblout hmmsearch results in tblout format")
    parser_required.add_argument("-n", "--hmmsearch_amr", type=str, required=True,  help = "path/to/hmmsearch.tblout hmmsearch results in tblout format")
    parser_required.add_argument("-a", "--hmmsearch_antifam", type=str, required=True,  help = "path/to/hmmsearch.tblout hmmsearch results in tblout format")

    parser_optional = parser.add_argument_group('Optional arguments')
    parser_optional.add_argument("-i","--identifier_mapping", type=str, required=False,  help = "path/to/identifier_mapping.tsv, Format: [id_orf]<tab>[id_contig]<tab>[id_mag], No header [Optional]")
    parser_optional.add_argument("-f","--fasta", type=str, required=False,  help = "path/to/gene_models.faa|ffn of ORFs [Optional]")
    parser_optional.add_argument("-t","--threshold", default=0.5, type=float, help = "taxopy fraction of classifications for consensus for weighted majority vote [Default: 0.5]")
    parser_optional.add_argument("-V", "--veba_database", type=str, required=False, help = "path/to/VEBA_DATABASE")
    parser_optional.add_argument("--hmmsearch_region", type=str, default="best",  help = "{best,full} Best domain or full sequence [Default: best]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Assertions
    assert_acceptable_arguments(opts.hmmsearch_region, {"best","full"})
    opts.hmmsearch_region = {"best":"best_domain", "full":"full_sequence"}[opts.hmmsearch_region]
    if opts.veba_database:
        import taxopy
        assert os.path.exists(opts.veba_database)
        assert os.path.isdir(opts.veba_database)
        assert os.path.exists(os.path.join(opts.veba_database, "Classify", "NCBITaxonomy", "nodes.dmp"))
        assert os.path.exists(os.path.join(opts.veba_database, "Classify", "NCBITaxonomy", "names.dmp"))
        assert os.path.exists(os.path.join(opts.veba_database, "Classify", "NCBITaxonomy", "merged.dmp"))

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
        df_identifier_mapping = pd.read_csv(opts.identifier_mapping, header=None, sep="\t")
        df_identifier_mapping.columns = ["id_protein", "id_contig", "id_genome"]
        df_identifier_mapping = df_identifier_mapping.set_index("id_protein")
        proteins = proteins | set(df_identifier_mapping.index)
        
    # Read annotations
    print(" * Reading Diamond table [nr]: {}".format(opts.diamond), file=sys.stderr)
    df_dmnd = pd.read_csv(opts.diamond, sep="\t", index_col=None, header=None)
    df_dmnd.columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle".split(" ")
    df_dmnd = df_dmnd.set_index("qseqid")
    df_dmnd = df_dmnd.loc[:,["sseqid", "pident", "length", "evalue", "bitscore", "staxids", "sscinames", "stitle"]]
    proteins = proteins | set(df_dmnd.index)

    print(" * Reading HMMSearch table [Pfam]: {}".format(opts.hmmsearch_pfam), file=sys.stderr)
    df_hmmsearch_pfam = read_hmmer(opts.hmmsearch_pfam, program="hmmsearch", format="tblout", verbose=False)[["identifier", opts.hmmsearch_region]].droplevel(0, axis=1)
    proteins = proteins | set(df_hmmsearch_pfam["target_name"])

    print(" * Reading HMMSearch table [AMR]: {}".format(opts.hmmsearch_amr), file=sys.stderr)
    df_hmmsearch_amr = read_hmmer(opts.hmmsearch_amr, program="hmmsearch", format="tblout", verbose=False)[["identifier", opts.hmmsearch_region]].droplevel(0, axis=1)
    proteins = proteins | set(df_hmmsearch_amr["target_name"])

    print(" * Reading HMMSearch table [AntiFam]: {}".format(opts.hmmsearch_antifam), file=sys.stderr)
    df_hmmsearch_antifam = read_hmmer(opts.hmmsearch_antifam, program="hmmsearch", format="tblout", verbose=False)[["identifier", opts.hmmsearch_region]].droplevel(0, axis=1)
    proteins = proteins | set(df_hmmsearch_antifam[ "target_name"])


    print(" * Reading KOFAMSCAN table [KEGG]: {}".format(opts.kofam), file=sys.stderr)
    df_kofamscan = pd.read_csv(opts.kofam, sep="\t", index_col=None, header=None)
    proteins = proteins | set(df_kofamscan.iloc[:,1])

    # All proteins
    proteins = sorted(proteins)
    print(" * Total proteins: {}".format(len(proteins)), file=sys.stderr)

    # Reindex Diamond
    df_dmnd = df_dmnd.reindex(proteins)

    # Pfam
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
    df_hmms_pfam = df_hmms_pfam.applymap(lambda x: x if isinstance(x,list) else [])
    df_hmms_pfam.index.name = "id_protein"
    df_hmms_pfam.insert(loc=0, column="number_of_hits", value=df_hmms_pfam["ids"].map(len))

    # AMR
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
    df_hmms_amr = df_hmms_amr.applymap(lambda x: x if isinstance(x,list) else [])
    df_hmms_amr.insert(loc=0, column="number_of_hits", value=df_hmms_amr["ids"].map(len))
    df_hmms_amr.index.name = "id_protein"

    # AntiFam
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
    df_hmms_antifam["hit"] = df_hmms_antifam["hit"].fillna(False)
    df_hmms_antifam = df_hmms_antifam.applymap(lambda x: x if isinstance(x,list) else [])
    df_hmms_antifam.insert(loc=0, column="number_of_hits", value=df_hmms_antifam["ids"].map(len))
    df_hmms_antifam.index.name = "id_protein"

    # KOFAMSCAN
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
    df_hmms_kofam = df_hmms_kofam.applymap(lambda x: x if isinstance(x,list) else [])
    df_hmms_kofam.insert(loc=0, column="number_of_hits", value=df_hmms_kofam["ids"].map(len))
    df_hmms_kofam.index.name = "id_protein"

    # Output table
    dataframes = list()
    for (name, df) in zip([ "nr", "Pfam", "NCBIfam-AMR", "AntiFam", "KOFAM"],[df_dmnd, df_hmms_pfam, df_hmms_amr, df_hmms_antifam, df_hmms_kofam]):
        df = df.copy()
        df.columns = df.columns.map(lambda x: (name,x))
        dataframes.append(df)
    df_protein_annotations = pd.concat(dataframes, axis=1)
    df_protein_annotations.index.name = "id_protein"

    if opts.identifier_mapping:
        df = df_identifier_mapping.copy()
        df.columns = df.columns.map(lambda x: ("Identifiers", x))
        df_protein_annotations = pd.concat([df, df_protein_annotations], axis=1)

    if not opts.veba_database:
        df_protein_annotations.to_csv(os.path.join(opts.output_directory, "protein_annotations.tsv.gz"), sep="\t")

    else:
        with open(os.path.join(opts.output_directory, "LINEAGE_DISCLAIMER.txt"), "w") as f:
            disclaimer = format_header("DISCLAIMER: Lineage predictions are NOT robust and DO NOT USE CORE MARKERS.  Please only use for exploratory suggestions.")
            print(disclaimer, file=f)
            print("For genome classification, please use VEBA's classify-prokaryotic.py, classify-eukaryotic.py, or classify-viral.py modules", file=f)
            print("For contig classification, try third party software such as CAT/BAT (https://github.com/dutilh/CAT)", file=f)

        print(" ", file=sys.stderr)
        print(disclaimer, file=sys.stderr)
        print("For genome classification, please use VEBA's classify-prokaryotic.py, classify-eukaryotic.py, or classify-viral.py modules", file=sys.stderr)
        print("For contig classification, try third party software such as CAT/BAT (https://github.com/dutilh/CAT)", file=sys.stderr)
        print(" ", file=sys.stderr)

        # Load taxonomy database
        taxdb = taxopy.TaxDb(
            nodes_dmp=os.path.join(opts.veba_database, "Classify", "NCBITaxonomy", "nodes.dmp"), 
            names_dmp=os.path.join(opts.veba_database, "Classify", "NCBITaxonomy", "names.dmp"),
            merged_dmp=os.path.join(opts.veba_database, "Classify", "NCBITaxonomy", "merged.dmp"),
        )

        # Get NCBI taxon id representatives
        invalid_taxids = set()

        protein_to_taxid = dict()
        for id_protein, id_taxon in pv(df_dmnd["staxids"].items(), description="Getting representative NCBI taxon id for each protein"):
            id_taxon = str(id_taxon).strip()
            if ";" in id_taxon:
                taxa = list()
                for id in id_taxon.split(";"):
                    id = int(id)
                    try:
                        taxon = taxopy.Taxon(id, taxdb)
                        taxa.append(taxon)
                    except taxopy.core.TaxidError:
                        pass
                if len(taxa) < 1:
                    protein_to_taxid[id_protein] = -1
                else:
                    if len(taxa) == 1:
                        protein_to_taxid[id_protein] = taxa[0].taxid
                    else:
                        consensus_taxa = taxopy.find_lca(taxa, taxdb)
                        protein_to_taxid[id_protein] = consensus_taxa.taxid
            else:
                if id_taxon == "nan":
                    id_taxon = -1
                else:
                    id_taxon = int(id_taxon)
                    if id_taxon not in taxdb.taxid2name:
                        invalid_taxids.add(id_taxon)
                        id_taxon = -1
                protein_to_taxid[id_protein] = int(id_taxon)
        protein_to_taxid = pd.Series(protein_to_taxid)
        print(" * The following identifiers are invalid with the NCBI Taxonomy database: {}".format(os.path.join(opts.veba_database, "Classify", "NCBITaxonomy", "*.dmp")), file=sys.stderr)
        x = 25
        if len(invalid_taxids) > x:
            print(str(list(invalid_taxids)[:x])[:-1],"...", file=sys.stderr)
        else:
            print(invalid_taxids, file=sys.stderr)

        # Get lineage for each taxid
        unique_taxids = protein_to_taxid[lambda x: x > -1].unique()
        taxonomic_ranks = ['superkingdom', 'genus', 'phylum', 'class', 'family', 'order','species']

        taxonid_to_lineage = dict()
        for id_taxon in pv(unique_taxids, description="Getting rank and lineage for each NCBI taxon id"):
            taxon = taxopy.Taxon(id_taxon, taxdb)
            lineage = pd.Series(taxon.rank_name_dictionary).reindex(taxonomic_ranks)
            taxonid_to_lineage[id_taxon] = lineage
        taxonid_to_lineage[-1] = pd.Series([np.nan]*len(taxonomic_ranks), index=taxonomic_ranks, dtype=object)
        
        # Get lineage for each protein
        protein_to_lineage = dict()
        for id_protein, id_taxon in pv(protein_to_taxid.items(), description="Getting rank and lineage for each protein id"):
            lineage = taxonid_to_lineage[id_taxon]
            protein_to_lineage[id_protein] = lineage
        df_lineage = pd.DataFrame(protein_to_lineage).T
        df_lineage.insert(loc=0, column="id_taxon", value=protein_to_taxid)
        df_lineage = df_lineage.reindex(proteins)
        df_lineage.loc[df_lineage["id_taxon"].isnull(), "id_taxon"] = -1
        df_lineage["id_taxon"] = df_lineage["id_taxon"].astype(int)
        df_lineage.index.name = "id_protein"

        df = df_lineage.copy()
        df.columns = df.columns.map(lambda x: ("NCBITaxonomy", x))
        df_protein_annotations = pd.concat([df_protein_annotations, df], axis=1)
        df_protein_annotations.to_csv(os.path.join(opts.output_directory, "protein_annotations.tsv.gz"), sep="\t")

        if opts.identifier_mapping:
            # Identifier mappings
            proteins_with_taxids = set(protein_to_taxid[lambda x: x > -1].index)

            # Contigs
            contig_to_proteins = defaultdict(set)
            for id_protein, id_contig in df_identifier_mapping["id_contig"].items():
                contig_to_proteins[id_contig].add(id_protein)

            # Get taxonomy for contigs
            contig_to_taxid = dict()
            contig_to_score = dict()
            contig_to_lineage = dict()
            contig_to_ncomponents = dict()
            for id_contig, components in pv(contig_to_proteins.items(), "Getting lineage for each contig"):

                taxa = list()
                weights = list()
                for id_taxon, score in df_dmnd.loc[list(components & proteins_with_taxids)].groupby(protein_to_taxid)["bitscore"].sum().items():
                    taxon = taxopy.Taxon(id_taxon, taxdb)
                    taxa.append(taxon)
                    weights.append(score)
                    
                if len(taxa) == 0:
                    pass
                else:
                    if len(taxa) == 1:
                        contig_to_taxid[id_contig] = taxa[0].taxid
                        contig_to_score[id_contig] = 1.0
                        lineage = str(taxa[0])
                    else:
                        weighted_majority_vote = taxopy.find_majority_vote(taxa, taxdb, weights=weights, fraction=opts.threshold)
                        contig_to_taxid[id_contig] = weighted_majority_vote.taxid
                        contig_to_score[id_contig] = weighted_majority_vote.agreement
                        lineage = str(weighted_majority_vote)
                    if lineage.startswith("s__"):
                        lineage = "d__" + lineage[1:]
                    contig_to_lineage[id_contig] = lineage
                    contig_to_ncomponents[id_contig] = len(taxa)

            contig_to_taxid = pd.Series(contig_to_taxid)
            contig_to_score = pd.Series(contig_to_score)
            contig_to_lineage = pd.Series(contig_to_lineage)
            contig_to_ncomponents = pd.Series(contig_to_ncomponents)

            df_lineage_contig = pd.DataFrame(
                OrderedDict([
                    ("number_of_proteins_with_lineage", contig_to_ncomponents),
                    ("id_taxon", contig_to_taxid),
                    ("score", contig_to_score),
                    ("lineage", contig_to_lineage),
                ]
            ))
            df_lineage_contig = df_lineage_contig.reindex(sorted(contig_to_proteins.keys()))
            df_lineage_contig.loc[df_lineage_contig["number_of_proteins_with_lineage"].isnull(), "number_of_proteins_with_lineage"] = 0
            df_lineage_contig["number_of_proteins_with_lineage"] = df_lineage_contig["number_of_proteins_with_lineage"].astype(int)
            df_lineage_contig.loc[df_lineage_contig["id_taxon"].isnull(), "id_taxon"] = -1
            df_lineage_contig["id_taxon"] = df_lineage_contig["id_taxon"].astype(int)
            df_lineage_contig.index.name = "id_contig"
            df_lineage_contig.to_csv(os.path.join(opts.output_directory, "lineage.weighted_majority_vote.contigs.tsv.gz"), sep="\t")

            # genomes
            genome_to_proteins = defaultdict(set)
            for id_protein, id_genome in df_identifier_mapping["id_genome"].items():
                genome_to_proteins[id_genome].add(id_protein)
                
                
            # Get taxonomy for genomes
            genome_to_taxid = dict()
            genome_to_score = dict()
            genome_to_lineage = dict()
            genome_to_ncomponents = dict()
            for id_genome, components in pv(genome_to_proteins.items(), "Getting lineage for each genome"):

                taxa = list()
                weights = list()
                for id_taxon, score in df_dmnd.loc[list(components & proteins_with_taxids)].groupby(protein_to_taxid)["bitscore"].sum().items():
                    taxon = taxopy.Taxon(id_taxon, taxdb)
                    taxa.append(taxon)
                    weights.append(score)
                    
                if len(taxa) == 0:
                    pass
                else:
                    if len(taxa) == 1:
                        genome_to_taxid[id_genome] = taxa[0].taxid
                        genome_to_score[id_genome] = 1.0
                        lineage = str(taxa[0])
                    else:
                        weighted_majority_vote = taxopy.find_majority_vote(taxa, taxdb, weights=weights, fraction=opts.threshold)
                        genome_to_taxid[id_genome] = weighted_majority_vote.taxid
                        genome_to_score[id_genome] = weighted_majority_vote.agreement
                        lineage = str(weighted_majority_vote)
                    if lineage.startswith("s__"):
                        lineage = "d__" + lineage[1:]
                    genome_to_lineage[id_genome] = lineage
                    genome_to_ncomponents[id_genome] = len(taxa)

            genome_to_taxid = pd.Series(genome_to_taxid)
            genome_to_score = pd.Series(genome_to_score)
            genome_to_lineage = pd.Series(genome_to_lineage)
            genome_to_ncomponents = pd.Series(genome_to_ncomponents)

            df_lineage_genome = pd.DataFrame(
                OrderedDict([
                    ("number_of_proteins_with_lineage", genome_to_ncomponents),
                    ("id_taxon", genome_to_taxid),
                    ("score", genome_to_score),
                    ("lineage", genome_to_lineage),
                ]
            ))
            df_lineage_genome = df_lineage_genome.reindex(sorted(genome_to_proteins.keys()))
            df_lineage_genome.loc[df_lineage_genome["number_of_proteins_with_lineage"].isnull(), "number_of_proteins_with_lineage"] = 0
            df_lineage_genome["number_of_proteins_with_lineage"] = df_lineage_genome["number_of_proteins_with_lineage"].astype(int)
            df_lineage_genome.loc[df_lineage_genome["id_taxon"].isnull(), "id_taxon"] = -1
            df_lineage_genome["id_taxon"] = df_lineage_genome["id_taxon"].astype(int)
            df_lineage_genome.index.name = "id_genome"
            df_lineage_genome.to_csv(os.path.join(opts.output_directory, "lineage.weighted_majority_vote.genomes.tsv.gz"), sep="\t")




if __name__ == "__main__":
    main()
    
                
