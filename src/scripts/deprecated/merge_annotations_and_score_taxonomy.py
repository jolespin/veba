#!/usr/bin/env python
import sys, os, argparse 
from collections import defaultdict
import pandas as pd
import numpy as np
from soothsayer_utils import read_hmmer
# from ete3 import NCBITaxa
from tqdm import tqdm  

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.2.27"

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
    parser.add_argument("-d","--diamond", type=str, required=True,  help = "path/to/diamond_nr.blast6 in blast6 format (No header) with the following fields: \nqseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle")
    parser.add_argument("-k","--kofam", type=str, required=True,  help = "path/to/kofamscan.tblout hmmsearch results in tblout format")
    parser.add_argument("-p", "--hmmsearch_pfam", type=str, required=True,  help = "path/to/hmmsearch.tblout hmmsearch results in tblout format")
    parser.add_argument("-n", "--hmmsearch_amr", type=str, required=True,  help = "path/to/hmmsearch.tblout hmmsearch results in tblout format")
    parser.add_argument("-a", "--hmmsearch_antifam", type=str, required=True,  help = "path/to/hmmsearch.tblout hmmsearch results in tblout format")
    parser.add_argument("-i","--identifier_mapping", type=str, required=False,  help = "path/to/identifier_mapping.tsv, Format: [id_orf]<tab>[id_contig]<tab>[id_mag], No header [Optional]")
    parser.add_argument("-f","--fasta", type=str, required=False,  help = "path/to/gene_models.faa|ffn of ORFs [Optional]")
    parser.add_argument("-o","--output_directory", type=str, default="veba_output/annotations",  help = "Output directory for annotations [Default: veba_output/annotations]")
    parser.add_argument("--veba_database", type=str, required=False, help = "path/to/VEBA_DATABASE")

    parser.add_argument("--hmmsearch_region", type=str, default="best",  help = "{best,full} Best domain or full sequence [Default: best]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Make directories and file objects
    os.makedirs(opts.output_directory, exist_ok=True)

    # taxa.sqlite
    if opts.database_taxa is None:
        opts.database_taxa = os.path.join(os.environ["VEBA_DATABASE"],  "Classify", "NCBITaxonomy", "taxa.sqlite")

    # Read Fasta File
    orfs = set()
    if opts.fasta:
        from Bio.SeqIO.FastaIO import SimpleFastaParser
        # ===
        # Protein
        # ===
        print(" * Parsing ORF fasta file: {}".format(opts.fasta), file=sys.stderr)
        # Parse
        with open(opts.fasta, "r") as f:
            for header, seq in SimpleFastaParser(f):
                id_orf = header.split(" ")[0]
                orfs.add(id_orf)

    if opts.identifier_mapping:
        print(" * Reading identifier mapping table: {}".format(opts.identifier_mapping), file=sys.stderr)
        df_identifiers = pd.read_csv(opts.identifier_mapping, header=None, sep="\t")
        df_identifiers.columns = ["id_orf", "id_contig", "id_mag"]
        df_identifiers = df_identifiers.set_index("id_orf")

        orfs = orfs | set(df_identifiers.index)

    print(" * Reading Diamond table [nr]: {}".format(opts.diamond), file=sys.stderr)
    df_dmnd = pd.read_csv(opts.diamond, sep="\t", index_col=None, header=None)
    df_dmnd.columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle".split(" ")
    df_dmnd = df_dmnd.set_index("qseqid")
    df_dmnd = df_dmnd.loc[:,["sseqid", "pident", "length", "evalue", "bitscore", "staxids", "sscinames", "stitle"]]
    orfs = orfs | set(df_dmnd.index)

    print(" * Reading HMMSearch table [Pfam]: {}".format(opts.hmmsearch_pfam), file=sys.stderr)
    df_hmmsearch_pfam = read_hmmer(opts.hmmsearch_pfam, program="hmmsearch", format="tblout", verbose=False)
    orfs = orfs | set(df_hmmsearch_pfam[("identifier", "target_name")])

    print(" * Reading HMMSearch table [AMR]: {}".format(opts.hmmsearch_amr), file=sys.stderr)
    df_hmmsearch_amr = read_hmmer(opts.hmmsearch_amr, program="hmmsearch", format="tblout", verbose=False)
    orfs = orfs | set(df_hmmsearch_amr[("identifier", "target_name")])

    print(" * Reading HMMSearch table [AntiFam]: {}".format(opts.hmmsearch_antifam), file=sys.stderr)
    df_hmmsearch_antifam = read_hmmer(opts.hmmsearch_antifam, program="hmmsearch", format="tblout", verbose=False)
    orfs = orfs | set(df_hmmsearch_antifam[("identifier", "target_name")])

    print(" * Reading KOFAMSCAN table [KEGG]: {}".format(opts.kofam), file=sys.stderr)
    df_kofamscan = pd.read_csv(opts.kofam, sep="\t", index_col=None, header=None)
    orfs = orfs | set(df_kofamscan.iloc[:,1])

    # For pandas indexing
    orfs = list(orfs)

    # Process HMMSearch
    # Pfam
    orf_to_hmmid = defaultdict(list)
    orf_to_hmmname = defaultdict(list)
    for i, row in tqdm(df_hmmsearch_pfam["identifier"].iterrows(), "Parsing HMMSearch hits", total=df_hmmsearch_pfam.shape[0], unit=" search hits"):
        id_orf, name_hmm, id_hmm = row[["target_name", "query_name", "query_accession"]]
        orf_to_hmmid[id_orf].append(id_hmm)
        orf_to_hmmname[id_orf].append(name_hmm)
    orf_to_hmmid = pd.Series(orf_to_hmmid).reindex(orfs).map(lambda x: [] if x is np.nan else x)
    orf_to_hmmname = pd.Series(orf_to_hmmname).reindex(orfs).map(lambda x: [] if x is np.nan else x)
    df_hmms_pfam = pd.concat([orf_to_hmmid.to_frame("ids"), orf_to_hmmname.to_frame("descriptions")], axis=1)

    # AMR
    orf_to_hmmid = defaultdict(list)
    orf_to_hmmname = defaultdict(list)
    for i, row in tqdm(df_hmmsearch_amr["identifier"].iterrows(), "Parsing HMMSearch hits", total=df_hmmsearch_amr.shape[0], unit=" search hits"):
        id_orf, name_hmm, id_hmm = row[["target_name", "query_name", "query_accession"]]
        orf_to_hmmid[id_orf].append(id_hmm)
        orf_to_hmmname[id_orf].append(name_hmm)
    orf_to_hmmid = pd.Series(orf_to_hmmid).reindex(orfs).map(lambda x: [] if x is np.nan else x)
    orf_to_hmmname = pd.Series(orf_to_hmmname).reindex(orfs).map(lambda x: [] if x is np.nan else x)
    df_hmms_amr = pd.concat([orf_to_hmmid.to_frame("ids"), orf_to_hmmname.to_frame("descriptions")], axis=1)

    # AntiFam
    orf_to_hmmid = defaultdict(list)
    orf_to_hmmname = defaultdict(list)
    for i, row in tqdm(df_hmmsearch_antifam["identifier"].iterrows(), "Parsing HMMSearch hits", total=df_hmmsearch_antifam.shape[0], unit=" search hits"):
        id_orf, name_hmm, id_hmm = row[["target_name", "query_name", "query_accession"]]
        orf_to_hmmid[id_orf].append(id_hmm)
        orf_to_hmmname[id_orf].append(name_hmm)
    orf_to_hmmid = pd.Series(orf_to_hmmid).reindex(orfs).map(lambda x: [] if x is np.nan else x)
    orf_to_hmmname = pd.Series(orf_to_hmmname).reindex(orfs).map(lambda x: [] if x is np.nan else x)
    df_hmms_antifam = pd.concat([orf_to_hmmid.to_frame("ids"), orf_to_hmmname.to_frame("descriptions")], axis=1)

    # Process KOFAMSCAN
    orf_to_koid = defaultdict(list)
    orf_to_koname = defaultdict(list)
    for i, row in tqdm(df_kofamscan.iterrows(), "Parsing KOFAMSCAN hits", total=df_kofamscan.shape[0], unit=" search hits"):
        significance, id_orf, name_hmm, id_hmm = row.iloc[[0, 1,6,2]]
        if significance == "*":
            orf_to_koid[id_orf].append(id_hmm)
            orf_to_koname[id_orf].append(name_hmm)
    orf_to_koid = pd.Series(orf_to_koid).reindex(orfs).map(lambda x: [] if x is np.nan else x)
    orf_to_koname = pd.Series(orf_to_koname).reindex(orfs).map(lambda x: [] if x is np.nan else x)
    df_kofam = pd.concat([orf_to_koid.to_frame("ids"), orf_to_koname.to_frame("descriptions")], axis=1)

    # Taxonomy
    ncbi = NCBITaxa(dbfile=opts.database_taxa)

    # Get lineage from taxon identifier
    def get_lineage_from_taxonid(id_taxon, ncbi=ncbi):
        # Force integer
        id_taxon = int(id_taxon)

        try:
            # Get lineage in taxon ids
            lineage = ncbi.get_lineage(id_taxon)

            # Get ranks for taxon ids
            ranks = dict(filter(lambda x: x[1] != "no rank", ncbi.get_rank(lineage).items()))

            # Translate ranks
            lineage_translated = pd.Series(ncbi.get_taxid_translator(ranks.keys()))
            lineage_translated.index = lineage_translated.index.map(lambda x:ranks[x])

            # Return lineage
            return lineage_translated.reindex(['superkingdom', 'genus', 'phylum', 'class', 'family', 'order','species'])
        except ValueError:
            return np.nan

    # Taxa identifiers
    taxa_ids = pd.Index(df_dmnd["staxids"].dropna().map(lambda x: x.split(";")[0]).unique())
    taxonid_to_lineage = dict()
    for id_taxon in tqdm(taxa_ids, "Querying lineage ranks", unit=" Taxa"):
        taxonid_to_lineage[id_taxon] = get_lineage_from_taxonid(id_taxon)  
    taxonid_to_lineage = pd.DataFrame(taxonid_to_lineage).T

    # Getting lineage for orf hits"
    df_lineage = dict()
    orf_to_taxonid = df_dmnd["staxids"].dropna().map(lambda x: x.split(";")[0])
    for id_orf, id_taxon in tqdm(orf_to_taxonid.items(), "Getting lineage for orf hits", unit=" orf", total=orf_to_taxonid.size):
        df_lineage[id_orf] = taxonid_to_lineage.loc[id_taxon]
    df_lineage = pd.DataFrame(df_lineage).T

    if not opts.identifier_mapping: 
        # Annotation
        dataframes = list()
        for (name, df) in zip([ "nr", "Pfam", "NCBIfam-AMR", "AntiFam", "KOFAM", "NCBITaxon"],[df_dmnd, df_hmms_pfam, df_hmms_amr, df_hmms_antifam, df_kofam, df_lineage]):
            df = df.copy()
            df.columns = df.columns.map(lambda x: (name,x))
            dataframes.append(df)
        df_annot_orfs = pd.concat(dataframes, axis=1)
        df_annot_orfs.index.name = "id_orf"
        df_annot_orfs.to_csv(os.path.join(opts.output_directory, "annotations.orfs.tsv.gz"), sep="\t")

    else: # Change this to taxopy instead of ETE3
        # Annotation
        dataframes = list()
        for (name, df) in zip(["Identifiers", "nr", "Pfam",  "NCBIfam-AMR", "AntiFam", "KOFAM", "NCBITaxon"],[df_identifiers,  df_hmms_pfam, df_hmms_amr, df_hmms_antifam, df_kofam, df_lineage]):
            df = df.copy()
            df.columns = df.columns.map(lambda x: (name,x))
            dataframes.append(df)
        df_annot_orfs = pd.concat(dataframes, axis=1)
        df_annot_orfs.index.name = "id_orf"
        df_annot_orfs.to_csv(os.path.join(opts.output_directory, "annotations.orfs.tsv.gz"), sep="\t")

        # Predict taxonomy for contigs and MAGs from NR
        contig_dataframes = list()
        mag_dataframes = list()

        for taxon_level in df_lineage.columns[::-1]:
            df = df_annot_orfs[[("Identifiers", "id_contig"), ("Identifiers", "id_mag"), ("NCBITaxon", taxon_level), ("nr", "bitscore")]].dropna(how="any", axis=0)
            
            contig_to_taxa_to_score = defaultdict(lambda: defaultdict(float))
            mag_to_taxa_to_score = defaultdict(lambda: defaultdict(float))
            
            # Scoring each taxa
            for id_orf, (id_contig, id_mag, taxa, bitscore) in tqdm(df.iterrows(), "Scoring taxa for contigs and MAGs [{}]".format(taxon_level), unit=" ORFs", total=df.shape[0]):
                
                if pd.notnull(taxa):
                    contig_to_taxa_to_score[id_contig][taxa] += bitscore
                    mag_to_taxa_to_score[id_mag][taxa] += bitscore
                    
            # Predicting taxonomy for contigs
            contig_scores = dict()
            for id_contig, scores in tqdm(contig_to_taxa_to_score.items(), "Predicting taxonomy for contigs [{}]".format(taxon_level)):
                scores = pd.Series(scores)
                total_scores = scores.sum()
                normalized_scores = scores/total_scores
                max_score = normalized_scores.max()
                idx_max = normalized_scores.index[normalized_scores == max_score]
                if idx_max.size == 1:
                    taxa_max = idx_max[0]
                else:
                    taxa_max = "INCONCLUSIVE[{}]".format("|".join(idx_max))
                contig_scores[id_contig] = {
                    "predicted_taxa":taxa_max,
                    "score":max_score,
                    "sum_bitscores":total_scores,
                    "info":normalized_scores.to_dict()
                }

            contig_scores = pd.DataFrame(contig_scores).T
            contig_scores.index.name = "id_contig"
            contig_scores.columns = contig_scores.columns.map(lambda x: (taxon_level, x))
            contig_dataframes.append(contig_scores)
            
            # Predicting taxonomy for MAGs
            mag_scores = dict()
            for id_mag, scores in tqdm(mag_to_taxa_to_score.items(), "Predicting taxonomy for MAGs [{}]".format(taxon_level)):
                scores = pd.Series(scores)
                total_scores = scores.sum()
                normalized_scores = scores/total_scores
                max_score = normalized_scores.max()
                idx_max = normalized_scores.index[normalized_scores == max_score]
                if idx_max.size == 1:
                    taxa_max = idx_max[0]
                else:
                    taxa_max = "INCONCLUSIVE[{}]".format("|".join(idx_max))
                mag_scores[id_mag] = {
                    "predicted_taxa":taxa_max,
                    "score":max_score,
                    "sum_bitscores":total_scores,
                    "info":normalized_scores.to_dict()
                }

            mag_scores = pd.DataFrame(mag_scores).T
            mag_scores.index.name = "id_mag"
            mag_scores.columns = mag_scores.columns.map(lambda x: (taxon_level, x))
            mag_dataframes.append(mag_scores)
            
            
        df_annot_contigs = pd.concat(contig_dataframes, axis=1)[df_lineage.columns]
        df_annot_contigs.columns.names = ["taxonomic_level", None]
        df_annot_contigs.columns = df_annot_contigs.columns.map(lambda x: x[::-1])
        df_annot_contigs = df_annot_contigs[["predicted_taxa", "score", "sum_bitscores", "info"]]

        df_annot_mags= pd.concat(mag_dataframes, axis=1)[df_lineage.columns]
        df_annot_mags.columns.names = ["taxonomic_level", None]
        df_annot_mags.columns = df_annot_mags.columns.map(lambda x: x[::-1])
        df_annot_mags = df_annot_mags[["predicted_taxa", "score", "sum_bitscores", "info"]]

        df_annot_contigs.to_csv(os.path.join(opts.output_directory, "lineage_predictions.contigs.tsv.gz"), sep="\t")
        df_annot_mags.to_csv(os.path.join(opts.output_directory, "lineage_predictions.mags.tsv.gz"), sep="\t")
            



    # else:
    #     # Annotation
    #     dataframes = list()
    #     for (name, df) in zip(["nr", "Pfam", "KOFAM"],[df_dmnd, df_pfam, df_kofam, df_lineage]):
    #         df = df.copy()
    #         df.columns = df.columns.map(lambda x: (name,x))
    #         dataframes.append(df)
    #     df_annot_orfs = pd.concat(dataframes, axis=1)
    #     df_annot_orfs.index.name = "id_orf"



if __name__ == "__main__":
    main()
    
                
