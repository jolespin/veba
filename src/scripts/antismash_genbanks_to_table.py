#!/usr/bin/env python
import sys, os, argparse, gzip, glob
from collections import OrderedDict
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.3.9"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <antismash_directory> -o <output.tsv[.gz]>".format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--antismash_directory", type=str, help = "path/to/antismash_directory")
    parser.add_argument("-n","--name", type=str, help = "Genome name")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")
    parser.add_argument("-e","--exclude_contig_edges", action="store_true", help = "Exclude clusters that are on contig edges")
    parser.add_argument("-s","--synopsis", type=str, help = "Summary file output")
    parser.add_argument("-t","--type_counts", type=str, help = "Type counts summary file output")
    parser.add_argument("-f","--fasta_output", type=str, help = "Fasta file output")
    parser.add_argument("--sep", type=str, default=" ", help = "Seperator between [id_gene]<sep>[bgc_description].  The space makes it a id and description.  If duplicate identifiers, you can separate with underscores (e.g., '__') [Default: <space>]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Output
    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Wildcard genbank files
    output = list()
    for fp in glob.glob(os.path.join(opts.antismash_directory, "*region*.gbk")):
        # Get genome and region identifiers
        if opts.name:
            id_genome = opts.name
        else:
            id_genome = fp.split("/")[-2]
        
        id_region = fp.split(".")[-2]
        # Iterate through sequence records
        for seq_record in SeqIO.parse(fp, "genbank"):
            # Iterate through sequence features
            id_contig = seq_record.id
            cds_features = list()
            product = None
            contig_edge = None
            for feature in tqdm(seq_record.features, "Parsing genbank file: {}".format(fp), unit=" features"):
                if feature.type == "CDS":
                    data = {"genome_id":id_genome, "region_id":id_region, "start":int(feature.location.start), "end":int(feature.location.end), "strand":feature.location.strand}
                    for k,v in feature.qualifiers.items():
                        if isinstance(v,list):
                            if len(v) == 1:
                                v = v[0]
                        data[k] = v
                    cds_features.append(pd.Series(data))
                elif feature.type == "region":                    
                    product = feature.qualifiers["product"][0]
                    contig_edge = feature.qualifiers["contig_edge"][0]
            df = pd.DataFrame(cds_features)
            
            # NCBI Genbank
            # Add gene id preferentially if one doesn't exist: gene_id, Name, ID, and locus_tag in that order.
            if "gene_id" not in df.columns:
                genes = None
                if "Name" in df.columns:
                    genes = df["Name"]
                elif "ID" in df.columns:
                    genes = df["ID"]
                elif "locus_tag" in df.columns:
                    genes = df["locus_tag"]
                elif "protein_id" in df.columns:
                    genes = df["protein_id"]
                elif "gene" in df.columns:
                    genes = df["gene"]

                df["gene_id"] = genes
            # If no contig identifiers, add it.  This prioritizes gff fields.
            if "contig_id" not in df.columns:
                df["contig_id"] = id_contig

            df.insert(0, "bgc_type", product)
            df.insert(1, "cluster_on_contig_edge", contig_edge)
            output.append(df)
    df_bgcs = pd.concat(output, axis=0).sort_values(["genome_id", "contig_id", "start", "end"])
    # for field in ["start","end"]:
    #     df_bgcs[field] = df_bgcs[field].astype(int)
    df_bgcs = df_bgcs.set_index(["genome_id", "contig_id", "region_id", "gene_id"]).sort_index()
    df_bgcs["translation"] = df_bgcs.pop("translation")

    # Add position on BGC
    i_to_position = OrderedDict()
    for (id_genome, id_contig, id_region), df in tqdm(df_bgcs.groupby(["genome_id", "contig_id", "region_id"])):
        df = df.sort_values(["start", "end"])
        d = dict(zip(df.index, range(1, df.shape[0] + 1)))
        i_to_position.update(d)
    i_to_position = pd.Series(i_to_position).astype(int)
    j = df_bgcs.columns.get_loc("start")
    df_bgcs.insert(loc=j, column="position_in_bgc", value=i_to_position)

    # BGC ID
    i_to_bgc = OrderedDict()
    for (id_genome, id_contig, id_region), df in tqdm(df_bgcs.groupby(["genome_id", "contig_id", "region_id"])):
        id_bgc = "|".join([id_genome, id_contig, id_region])
        d = dict(zip(df.index, [id_bgc]*df.shape[0]))
        i_to_bgc.update(d)
    i_to_bgc = pd.Series(i_to_bgc)
    df_bgcs.insert(loc=0, column="bgc_id", value=i_to_bgc)

    # Component
    i_to_component = OrderedDict()
    for i, row in df_bgcs.iterrows():
        id_component = "{}_{}|{}-{}({})".format(*row[["bgc_id", "position_in_bgc", "start", "end"]], {1:"+",-1:"-"}[row["strand"]])
        i_to_component[i] = id_component 
    i_to_component = pd.Series(i_to_component)
    df_bgcs.insert(loc=1, column="component_id", value=i_to_component)


    # df_bgcs = df_bgcs.set_index(drop=False)



    if opts.type_counts:
        tmp = df_bgcs.reset_index().set_index(["genome_id", "contig_id", "region_id", "bgc_type"]).index.unique().value_counts()
        value_counts = tmp.groupby(lambda x: (x[0], x[3])).sum()
        df_typecounts_with_edges = pd.DataFrame([[*x[0], x[1]] for x in value_counts.items()], columns=["id_genome", "bgc_type", "number_of_bgcs"]).set_index(["id_genome", "bgc_type"]).sort_values("number_of_bgcs", ascending=False)

        tmp = df_bgcs.loc[~df_bgcs["cluster_on_contig_edge"].map(eval).values].reset_index().set_index(["genome_id", "contig_id", "region_id", "bgc_type"]).index.unique().value_counts()
        value_counts = tmp.groupby(lambda x: (x[0], x[3])).sum()
        df_typecounts_noedges = pd.DataFrame([[*x[0], x[1]] for x in value_counts.items()], columns=["id_genome", "bgc_type", "number_of_bgcs(not_on_edge)"]).set_index(["id_genome", "bgc_type"]).sort_values("number_of_bgcs(not_on_edge)", ascending=False)

        df_type_counts = pd.concat([df_typecounts_with_edges, df_typecounts_noedges], axis=1).fillna(0).astype(int)
        df_type_counts.to_csv(opts.type_counts, sep="\t")

    if opts.synopsis:
        df = df_bgcs.reset_index(drop=False).set_index("bgc_id")

        df_synopsis = df_bgcs.reset_index(drop=False).set_index("component_id")["bgc_id"].value_counts().to_frame("number_of_genes")
        for field in ["genome_id", "contig_id", "region_id", "bgc_type", "cluster_on_contig_edge"]:
            df_synopsis[field] = pd.Series(df[field].to_dict())
        df_synopsis.index.name = "bgc_id"

        fields = ['genome_id', 'contig_id', 'region_id', 'bgc_type','cluster_on_contig_edge', 'number_of_genes']
        df_synopsis = df_synopsis.loc[:,fields]
        df_synopsis.to_csv(opts.synopsis, sep="\t")


    # Exclude contig edges
    if opts.exclude_contig_edges:
        df_bgcs = df_bgcs.loc[~df_bgcs["cluster_on_contig_edge"].map(eval).values]

    # Output
    df_bgcs.reset_index(drop=False).set_index(["bgc_id", "component_id"]).to_csv(opts.output, sep="\t")

    # Fasta
    if opts.fasta_output:
        if opts.fasta_output.endswith(".gz"):
            f = gzip.open(opts.fasta_output, "wt")
        else:
            f = open(opts.fasta_output, "w")
        for (id_genome, id_contig, id_region, id_gene), row in tqdm(df_bgcs.iterrows(), "Writing genes to fasta file: {}".format(opts.fasta_output)):
            id_component = row["component_id"]
            seq = row["translation"]

            header = "{}{}{}".format(
                id_gene,
                opts.sep, 
                id_component,
            )
            print(">{}\n{}".format(header, seq), file=f)
        f.close()

if __name__ == "__main__":
    main()
    
                
