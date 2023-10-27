#!/usr/bin/env python
from asyncio import protocols
from audioop import reverse
import sys, os, argparse, gzip, glob
from collections import OrderedDict, defaultdict
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.9.17"

def gc_content(seq):
    seq = seq.upper()
    number_of_gc = seq.count("G") + seq.count("C")
    return number_of_gc/len(seq)

# Reverse Complement

def reverse_complement(sequence, conversion = str.maketrans('ACGTacgt','TGCAtgca')):
    """
    Inputs a sequence of DNA and returns the reverse complement se$
    Outputs a reversed string of the input where A>T, T>A, G>C, an$
    """
    complement = sequence.translate(conversion)
    return complement[::-1]

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <antismash_directory> -o <output_directory>".format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--antismash_directory", type=str, help = "path/to/antismash_directory")
    parser.add_argument("-n","--name", type=str, help = "Genome name")
    parser.add_argument("-o","--output_directory", type=str, help = "path/to/output_directory [Default: --antismash_directory]")
    parser.add_argument("-e","--exclude_contig_edges", action="store_true", help = "Exclude BGCs that are on contig edges")
    parser.add_argument("--no_protein_fasta", type=str, help = "No protein fasta file output")
    parser.add_argument("--no_nucleotide_fasta", type=str, help = "No nucleotide fasta file output")
    parser.add_argument("--sample", type=str, help = "Sample of origin")
    parser.add_argument("--use_original_gene_ids", action="store_true", help = "Use original gene ids for the proteins in fasta/bgcs.faa.gz file")
    parser.add_argument("-g", "--gzipped_genbanks", action="store_true", help = "Use if genbanks are gzipped")

    # parser.add_argument("--separator_in_protein_header", type=str, default=" ", help = "Seperator between [id_gene]<sep>[bgc_description].  The space makes it a id and description.  If duplicate identifiers, you can separate with underscores (e.g., '__') [Default: <space>]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Output
    if not opts.output_directory:
        opts.output_directory = opts.antismash_directory

    os.makedirs(opts.output_directory, exist_ok=True)

    # Wildcard genbank files
    bgc_to_dna = dict()
    if opts.gzipped_genbanks:
        genbank_region_filepaths = glob.glob(os.path.join(opts.antismash_directory, "*region*.gbk.gz"))
    else:
        genbank_region_filepaths = glob.glob(os.path.join(opts.antismash_directory, "*region*.gbk"))

    if len(genbank_region_filepaths):
        output = list()
        for fp in genbank_region_filepaths:
            # Get genome and region identifiers
            if opts.name:
                id_genome = opts.name
            else:
                id_genome = fp.split("/")[-2]
            
            if opts.gzipped_genbanks:
                id_region = fp.split(".")[-3]
            else:
                id_region = fp.split(".")[-2]


            # id_contig = fp[:-14] # Won't handle ids that have illegal characters

            # Iterate through sequence records
            if opts.gzipped_genbanks:
                f = gzip.open(fp, "rt")
            else:
                f = open(fp, "r")
            for seq_record in SeqIO.parse(f, "genbank"):
                # Iterate through sequence features
                id_contig = seq_record.name.strip() #seq_record.id.strip() corrupts the ID (https://github.com/antismash/antismash/issues/651)
                cds_features = list()
                product = None
                contig_edge = None

                id_bgc = "|".join([id_genome, id_contig, id_region])
                
                for feature in tqdm(seq_record.features, "Parsing genbank file: {}".format(fp), unit=" features"):
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    strand = {1:"+", -1:"-"}[feature.location.strand]

                    if feature.type == "region":
                        product = feature.qualifiers["product"][0]
                        contig_edge = feature.qualifiers["contig_edge"][0]


                        sequence = str(seq_record.seq)[start:end]
                        if strand == "-":
                            sequence = reverse_complement(sequence)
                        bgc_to_dna[id_bgc] = sequence
                    else:
                        if feature.type == "CDS":
                            data = {"genome_id":id_genome, "contig_id":id_contig, "region_id":id_region, "start":start, "end":end, "strand":strand}
                            for k,v in feature.qualifiers.items():

                                if isinstance(v,list):
                                    if len(v) == 1:
                                        v = v[0]

                                data[k] = v
                            cds_features.append(pd.Series(data))


                df = pd.DataFrame(cds_features)

                # Fix missing values on "allorfs" genes predicted within antiSMASH that are not in GFF
                if "locus_tag" in df.columns:
                    locus_tags_with_allorf = df["locus_tag"].dropna()
                    index_with_allorf = locus_tags_with_allorf.index[locus_tags_with_allorf.map(lambda x: x.startswith("allorf"))]

                    # Excluding the allorfs, which gene name fields are in the genbank
                    if not index_with_allorf.empty:
                        df_noallorfs = df.drop(index_with_allorf, axis=0).dropna(how="all", axis=1)
                        gene_name_columns = df_noallorfs.loc[:,list(set(df_noallorfs.columns) & set([ "Name", "ID", "locus_tag", "protein_id", "gene", "gene_id"]))].columns
                        # missing_gene_name_columns = 
                        for i, row in df.loc[index_with_allorf].iterrows():
                            id_allorf = row["locus_tag"]
                            _, start, end = id_allorf.split("_")
                            # strand = {"1":"+", "-1":"-"}[str(row["strand"])]
                            strand = row["strand"]
                            id_gene = "{}_antiSMASH-{}:{}({})".format(id_contig, start, end, strand)
                            for id_field in gene_name_columns:
                                df.loc[i,id_field] = id_gene

                # NCBI Genbank
                # Add gene id preferentially if one doesn't exist: gene_id, Name, ID, and locus_tag in that order.
                if "gene_id" not in df.columns:
                    genes = None
                    if "Name" in df.columns:
                        genes = df["Name"].strip()
                    elif "ID" in df.columns:
                        genes = df["ID"].strip()
                    elif "locus_tag" in df.columns:
                        genes = df["locus_tag"].strip()
                    elif "protein_id" in df.columns:
                        genes = df["protein_id"].strip()
                    elif "gene" in df.columns:
                        genes = df["gene"].strip()
                    assert genes is not None, "Cannot identify gene identifiers for contig: {}".format(id_contig)
                    df["gene_id"] = genes

                # # If no contig identifiers, add it.  This prioritizes gff fields.
                # if "contig_id" not in df.columns:
                #     df["contig_id"] = id_contig.strip()

                df.insert(0, "protocluster_type", product)
                df.insert(1, "cluster_on_contig_edge", contig_edge)

                output.append(df)

            f.close()

        df_components = pd.concat(output, axis=0).sort_values(["genome_id", "contig_id", "start", "end"])

        for field in ["genome_id", "contig_id", "region_id", "gene_id"]:
            df_components[field] = df_components[field].map(lambda x: x.strip().replace(" ",""))

        for field in ["start","end"]:
            df_components[field] = df_components[field].astype(int)

        df_components = df_components.set_index(["genome_id", "contig_id", "region_id", "gene_id"]).sort_index()
        df_components["translation"] = df_components.pop("translation")

        # Add position on BGC
        i_to_position = OrderedDict()
        for (id_genome, id_contig, id_region), df in df_components.groupby(["genome_id", "contig_id", "region_id"]):
            df = df.sort_values(["start", "end"])
            d = dict(zip(df.index, range(1, df.shape[0] + 1)))
            i_to_position.update(d)
        i_to_position = pd.Series(i_to_position).astype(int)
        j = df_components.columns.get_loc("start")
        df_components.insert(loc=j, column="position_in_bgc", value=i_to_position)

        # BGC ID
        i_to_bgc = OrderedDict()
        for (id_genome, id_contig, id_region), df in df_components.groupby(["genome_id", "contig_id", "region_id"]):
            id_bgc = "|".join([id_genome, id_contig, id_region])
            d = dict(zip(df.index, [id_bgc]*df.shape[0]))
            i_to_bgc.update(d)
        i_to_bgc = pd.Series(i_to_bgc)
        df_components.insert(loc=0, column="bgc_id", value=i_to_bgc)

        # BGC to number of components
        bgc_to_numberofcomponents = df_components["bgc_id"].value_counts()

        # BGC to contig edge
        bgc_to_contigedge = dict(zip(df_components["bgc_id"], df_components["cluster_on_contig_edge"]))

        # Component
        i_to_component = OrderedDict()
        for i, row in df_components.iterrows():
            id_component = "{}_{}|{}:{}({})".format(*row[["bgc_id", "position_in_bgc", "start", "end"]], row["strand"])
            i_to_component[i] = id_component 
        i_to_component = pd.Series(i_to_component)
        df_components.insert(loc=1, column="component_id", value=i_to_component)

        # Add sample of origin if provided
        if opts.sample:
            df_components["sample_of_origin"] = opts.sample

        # Protocluster-type
        # ==================
        # With BGCs on edge
        tmp = df_components.reset_index().set_index(["genome_id", "contig_id", "region_id", "protocluster_type"]).index.unique().value_counts()
        
        value_counts = defaultdict(int)
        for (id_genome, _, _, protocluster_type), count in tmp.items():
            value_counts[(id_genome, protocluster_type)] += 1
        value_counts = pd.Series(value_counts)
            
        df_typecounts_with_edges = pd.DataFrame([[id_genome, protocluster_type, count] for ((id_genome, protocluster_type), count) in value_counts.items()], columns=["id_genome", "protocluster_type", "number_of_bgcs"]).set_index(["id_genome", "protocluster_type"]).sort_values("number_of_bgcs", ascending=False)

        # With only BGCs not on edge
        tmp = df_components.loc[~df_components["cluster_on_contig_edge"].map(eval).values].reset_index().set_index(["genome_id", "contig_id", "region_id", "protocluster_type"]).index.unique().value_counts()

        value_counts = defaultdict(int)
        for (id_genome, _, _, protocluster_type), count in tmp.items():
            value_counts[(id_genome, protocluster_type)] += 1
        value_counts = pd.Series(value_counts)

        df_typecounts_noedges = pd.DataFrame([[id_genome, protocluster_type, count] for ((id_genome, protocluster_type), count) in value_counts.items()], columns=["id_genome", "protocluster_type", "number_of_bgcs(not_on_edge)"]).set_index(["id_genome", "protocluster_type"]).sort_values("number_of_bgcs(not_on_edge)", ascending=False)

        # Merging
        df_type_counts = pd.concat([df_typecounts_with_edges, df_typecounts_noedges], axis=1).fillna(0).astype(int)

        # Add sample of origin if provided
        if opts.sample:
            df_type_counts["sample_of_origin"] = opts.sample
            
        df_type_counts.to_csv(os.path.join(opts.output_directory, "bgc_protocluster-types.tsv.gz"), sep="\t")

        # Output BGC summary
        # ==================
        df = df_components.reset_index(drop=False).set_index("bgc_id")

        df_bgcs = df_components.reset_index(drop=False).set_index("component_id")["bgc_id"].value_counts().to_frame("number_of_genes")
        for field in ["genome_id", "contig_id", "region_id", "protocluster_type", "cluster_on_contig_edge"]:
            df_bgcs[field] = pd.Series(df[field].to_dict())
        df_bgcs.index.name = "bgc_id"

        fields = ['genome_id', 'contig_id', 'region_id', 'protocluster_type','cluster_on_contig_edge', 'number_of_genes']
        df_bgcs = df_bgcs.loc[:,fields]
        
        # Add sample of origin if provided
        if opts.sample:
            df_bgcs["sample_of_origin"] = opts.sample

        # Exclude contig edges
        if opts.exclude_contig_edges:
            df_bgcs = df_bgcs.loc[~df_bgcs["cluster_on_contig_edge"].map(eval).values]

        print("Number of BGCs:", df_bgcs.shape[0], file=sys.stderr)
        df_bgcs.to_csv(os.path.join(opts.output_directory, "identifier_mapping.bgcs.tsv.gz"), sep="\t")


        # Output components summary
        # =========================
        # Exclude contig edges
        if opts.exclude_contig_edges:
            df_components = df_components.loc[~df_components["cluster_on_contig_edge"].map(eval).values]

        print("Number of components:", df_components.shape[0], file=sys.stderr)
        df_components.reset_index(drop=False).set_index("component_id").to_csv(os.path.join(opts.output_directory, "identifier_mapping.components.tsv.gz"), sep="\t")

        # Fasta (Protein)
        # ===============
        if not opts.no_protein_fasta:
            os.makedirs(os.path.join(opts.output_directory, "fasta"), exist_ok=True)

            fp = os.path.join(opts.output_directory, "fasta", "components.faa.gz")

            f = gzip.open(fp, "wt")

            if opts.use_original_gene_ids:
                for (id_genome, id_contig, id_region, id_gene), row in tqdm(df_components.iterrows(), "Writing genes to fasta file: {}".format(fp), total=df_components.shape[0]):
                    id_component = row["component_id"]
                    seq = row["translation"]

                    header = "{}{}{}".format(
                        id_gene,
                        opts.separator_in_protein_header, 
                        id_component,
                    )
                    print(">{}\n{}".format(header, seq), file=f)

            else:
                for (id_genome, id_contig, id_region, id_gene), row in tqdm(df_components.iterrows(), "Writing genes to fasta file: {}".format(fp), total=df_components.shape[0]):
                    id_component = row["component_id"]
                    seq = row["translation"]

                    header = "{} {}".format(
                        id_component,
                        id_gene,
                    )
                    print(">{}\n{}".format(header, seq), file=f)
            f.close()

        # Fasta (Nucleotide)
        # ==========
        if not opts.no_nucleotide_fasta:
            os.makedirs(os.path.join(opts.output_directory, "fasta"), exist_ok=True)
            fp = os.path.join(opts.output_directory, "fasta", "bgcs.fasta.gz")

            f = gzip.open(fp, "wt")

            for (id_bgc, seq) in tqdm(bgc_to_dna.items(), "Writing nucleotides to fasta file: {}".format(fp)):

                description = "len={};gc={:.3};n_genes={};edge={}".format(
                        len(seq),
                        gc_content(seq),
                        bgc_to_numberofcomponents[id_bgc],
                        bgc_to_contigedge[id_bgc],
                )

                header = "{} {}".format(
                    id_bgc,
                    description,
                )
                print(">{}\n{}".format(header, seq), file=f)
            f.close()
    else:
        print("No BGC regions detected: {}".format(opts.antismash_directory), file=sys.stderr)

if __name__ == "__main__":
    main()
    
                
