#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
import pandas as pd
from ete3 import NCBITaxa
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.08.09"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__

    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <diamond> -t <orf_to_mag> -o <output> --heatmap_output <pdf> ".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i","--diamond", type=str, required=True, help = "Diamond blast6 format. No header")
    parser_io.add_argument("-t","--orf_to_mag", type=str, required=True, help = "Tab-separated values file with [orf]<tab>[mag].  No header.")
    parser_io.add_argument("-d", "--database", type=str, default=os.path.join(script_directory, "../db"), help="Diamond | Marker database [Default: install_directory']")
    parser_io.add_argument("-o", "--output", type=str, default="stdout", help="Prediction output [Default: stdout']")
    parser_io.add_argument("-p", "--prokaryotes", type=str, help="Prokaryotes list output")
    parser_io.add_argument("-e", "--eukaryotes", type=str, help="Eukaryotes list output")
    parser_io.add_argument( "--diamond_header", 
                            type=str, 
                            default="qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,staxids,sscinames,stitle", 
                            help="Diamond | Comma separated values for diamond blast6 format headers [Default: qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,staxids,sscinames,stitle]",
    )
    

    parser_plotting = parser.add_argument_group('Required plotting arguments')
    parser_plotting.add_argument("--figsize", type=str, default="infer", help="Matplotlib | figsize [Default: 'infer']")
    parser_plotting.add_argument("--colormap", type=str, default="Greys", help="Matplotlib | colormap [Default: 'Gray']")
    parser_plotting.add_argument("--heatmap_output", type=str, help="Matplotlib | Heatmap output")
    parser_plotting.add_argument("--title", type=str, help="Matplotlib | Heatmap title")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Get NCBI Taxonomy server
    print("Loading NCBI Taxonomy Database:", os.path.join(opts.database, "NCBITaxonomy", "taxa.sqlite"), file=sys.stderr)
    ncbi = NCBITaxa(
        dbfile=os.path.join(opts.database, "NCBITaxonomy", "taxa.sqlite"), 
    )

    # Superkingdom from lineage
    def get_superkingdom_from_taxonid(id_taxon, ncbi=ncbi):
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

            # Return superkingdom
            return lineage_translated["superkingdom"]
        except ValueError:
            return np.nan

    # ID Mapping
    print("Reading and formatting ORFs to MAGs table:", opts.orf_to_mag, file=sys.stderr)
    orf_to_mag = pd.read_csv(opts.orf_to_mag, index_col=0, sep="\t", header=None).iloc[:,0]

    # Format Diamond results
    print("Reading and formatting Diamond results:", opts.diamond, file=sys.stderr)
    df_dmnd = pd.read_csv(opts.diamond, sep="\t", index_col=None, header=None)
    opts.diamond_header = opts.diamond_header.strip().split(",")
    assert {"qseqid", "bitscore", "staxids"} <= set(opts.diamond_header), "Diamond header must have at least 'qseqid','bitscore', and 'staxids'"
    df_dmnd.columns = opts.diamond_header
    df_dmnd = df_dmnd.set_index("qseqid").loc[:,["bitscore", "staxids"]].dropna(how="any", axis=0)
    df_dmnd["staxids_representative"] = df_dmnd["staxids"].map(lambda staxids: staxids.split(";")[0])


    # Get superkingdom
    unique_taxonids = df_dmnd["staxids_representative"].unique()
    taxonid_to_superkingdom = dict() #pd.Series(data=unique_taxonids,index=unique_taxonids) # Hacky but just want to get a container so I can update the labels w/ superkingom and keep the order
    # taxonid_to_superkingdom = taxonid_to_superkingdom.map(lambda id_taxon: get_superkingdom_from_taxonid(id_taxon))
    for id_taxon in tqdm(unique_taxonids, "Getting superkingdom for each NCBITaxonID", unit=" NCBITaxonID"):
        taxonid_to_superkingdom[id_taxon] = get_superkingdom_from_taxonid(id_taxon)
    taxonid_to_superkingdom = pd.Series(taxonid_to_superkingdom)

    # Update the Diamond results with superkingdom
    df_dmnd["superkingdom"] = df_dmnd["staxids_representative"].map(lambda id_taxon: taxonid_to_superkingdom[id_taxon])

    # Get scores
    mag_to_scores = dict()
    for id_mag, df in tqdm(df_dmnd.groupby(orf_to_mag), "Getting superkingdom scores for each MAG", unit=" MAG"):
        mag_to_scores[id_mag] = df.groupby("superkingdom").sum().reindex(["Archaea", "Bacteria", "Eukaryota"]).fillna(0.0)["bitscore"]
    df_scores = pd.DataFrame(mag_to_scores).T.astype(float)

    # Get predictions
    print("Calculating probabilities for each superkingdom", file=sys.stderr)
    Y_hat = df_scores/df_scores.sum(axis=1).values.reshape((-1,1))
    y_hat = Y_hat.idxmax(axis=1)

    # Add multiindex
    df_scores.columns = df_scores.columns.map(lambda x: ("Scores", x))
    Y_hat.columns = Y_hat.columns.map(lambda x: ("Probabilities", x))

    # Output
    print("Creating prediction output file:", opts.output, file=sys.stderr)
    df_output = pd.concat([
        y_hat.to_frame(("Prediction", "Superkingdom")),
        Y_hat,
        df_scores,
    ], axis=1)
    df_output.index.name = "Genomes"

    df_output.to_csv({True:sys.stdout, False:opts.output}[opts.output == "stdout"], sep="\t")

    # Prokaryotes
    if opts.prokaryotes:
        print("Creating Prokaryote list file:", opts.prokaryotes, file=sys.stderr)
        with open(opts.prokaryotes, "w") as f:
            for id_mag in y_hat[y_hat.map(lambda y: y in {"Archaea", "Bacteria"})].index:
                print(id_mag, file=f)
    # Prokaryotes
    if opts.eukaryotes:
        print("Creating Eukaryote list file:", opts.eukaryotes, file=sys.stderr)
        with open(opts.eukaryotes, "w") as f:
            for id_mag in y_hat[y_hat.map(lambda y: y in { "Eukaryota"})].index:
                print(id_mag, file=f)

    # Heatmap
    if opts.heatmap_output:
        print("Creating heatmap of prediction probabilities:", opts.heatmap_output, file=sys.stderr)

        import matplotlib.pyplot as plt 
        import seaborn as sns 

        df_heatmap = df_output["Probabilities"]
        if opts.figsize == "infer":
            n, m = df_heatmap.shape
            opts.figsize = (m*1.618, n*0.5)
        else:
            opts.figsize = opts.figsize.strip().split(",")
            opts.figsize = (float(opts.figsize[0]), float(opts.figsize([1])))

        with plt.style.context("seaborn-white"):
            fig, ax = plt.subplots(figsize=opts.figsize)
            sns.heatmap(df_heatmap, vmin=0, vmax=1, cmap=opts.colormap, annot=True, ax=ax, edgecolor="white", linewidth=1, cbar_kws={"label":"Probability"})
            ax.set_xticklabels(ax.get_xticklabels(), fontsize=15)
            ax.set_yticklabels(ax.get_yticklabels(), fontsize=12)
            ax.set_xlabel("Superkingdom", fontsize=15, fontweight="bold")
            ax.set_ylabel("Genomes", fontsize=15, fontweight="bold")
            if opts.title:
                ax.set_title(opts.title, fontsize=15, fontweight="bold")
            fig.savefig(opts.heatmap_output, format="pdf", bbox_inches="tight", dpi=300)

if __name__ == "__main__":
    main()
