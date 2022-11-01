#!/usr/bin/env python
import sys, os, argparse, gzip, glob
from collections import OrderedDict
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.11.01"

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
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")
    parser.add_argument("-e","--exclude_contig_edges", action="store_true", help = "Exclude clusters that are on contig edges")

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
        id_genome = fp.split("/")[-2]
        id_region = fp.split(".")[-2]
        # Iterate through sequence records
        for seq_record in SeqIO.parse(fp, "genbank"):
            # Iterate through sequence features
            id_contig = seq_record.id
            cds_features = list()
            product = None
            contig_edge = None
            for feature in tqdm(seq_record.features, "Parsing genbanke file: {}".format(fp), unit=" features"):
                if feature.type == "CDS":
                    data = {"genome_id":id_genome, "region_id":id_region}
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
            df.insert(0, "bgc_type", product)
            df.insert(1, "cluster_on_contig_edge", contig_edge)
            output.append(df)
    df_bgcs = pd.concat(output, axis=0).set_index(["genome_id", "contig_id", "region_id", "gene_id"]).sort_index()
    df_bgcs["translation"] = df_bgcs.pop("translation")

    # Exclude contig edges
    if opts.exclude_contig_edges:
        df_bgcs = df_bgcs.loc[~df_bgcs["cluster_on_contig_edge"].values.astype(bool)]

    # Output
    df_bgcs.to_csv(opts.output, sep="\t")

if __name__ == "__main__":
    main()
    
                
