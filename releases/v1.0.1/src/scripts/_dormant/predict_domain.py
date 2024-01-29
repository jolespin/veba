#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, warnings
import pandas as pd
from ete3 import NCBITaxa

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.10.08"

DATABASE="/usr/local/scratch/CORE/jespinoz/db/veba/v2021.10.08"
TAXASQLITE = os.path.join(DATABASE, "NCBITaxonomy", "taxa.sqlite")


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__

    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <mags> -o <output> --heatmap_output <pdf> ".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i","--comparesketch_directory", type=str, help = "Directory of comparesketch tables")
    parser_io.add_argument("-x","--extension", type=str, default="txt", help = "Extension of comparesketch.sh output [Default: txt]")
    parser_io.add_argument("-o","--output_directory", type=str, default="domain_predictions_output", help = "path/to/output_directory/ [Default: domain_predictions_output]")
    parser_io.add_argument("-d", "--taxa_sqlite", type=str, default=TAXASQLITE, help="path/to/taxa.sqlite [Default: {}']".format(TAXASQLITE))
    parser_io.add_argument("-m", "--metric", type=str, default="WKID", help="Metric to use for weight.  Choose between {WKID,KID,ANI} [Default: WKID']")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Plotting
    parser_plotting = parser.add_argument_group('Required plotting arguments')
    parser_plotting.add_argument("--no_heatmap", action="store_true", help="Matplotlib | Heatmap output")
    parser_plotting.add_argument("--figsize", type=str, default="infer", help="Matplotlib | figsize [Default: 'infer']")
    parser_plotting.add_argument("--colormap", type=str, default="Greys", help="Matplotlib | colormap [Default: 'Gray']")
    parser_plotting.add_argument("--style", type=str, default="seaborn-white", help="Matplotlib | Style [Default: 'seaborn-white']")
    parser_plotting.add_argument("--title", type=str, help="Matplotlib | Heatmap title")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Output directory
    os.makedirs(opts.output_directory, exist_ok=True)


    # Get NCBI Taxonomy server
    print("Loading NCBI Taxonomy Database:", opts.taxa_sqlite, file=sys.stderr)
    ncbi = NCBITaxa(
        dbfile=opts.taxa_sqlite, 
    )

    # Domain from lineage
    def get_domain_from_taxonid(id_taxon, ncbi=ncbi):
        # Force integer
        id_taxon = int(id_taxon)

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
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

    # Input filepaths 
    filepaths = glob.glob(os.path.join(opts.comparesketch_directory, "*.{}".format(opts.extension)))
    assert len(filepaths), "No sketch files using extension='{}' available in {}".format(opts.extension, opts.comparesketch_directory)
    # Get weights
    mag_to_weights = dict()
    for fp in filepaths:
        id_mag = fp.split("/")[-1][:-1*(len(opts.extension) + 1)]
        domain_to_weight = {
            "Archaea":0.0,
            "Bacteria":0.0,
            "Eukaryota":0.0,
            "Viruses":0.0,
        }
        df = pd.read_csv(fp, sep="\t",  skiprows=2)
        df["Domain"] = df["TaxID"].map(lambda id_taxon: get_domain_from_taxonid(id_taxon))
        df[opts.metric] = df[opts.metric].map(lambda x: float(x[:-1]))
        for i, row in df.iterrows():
            id_domain = row["Domain"]
            weight = row[opts.metric]
            domain_to_weight[id_domain] += weight
        mag_to_weights[id_mag] = pd.Series(domain_to_weight)
    df_weights = pd.DataFrame(mag_to_weights).T
    df_weights.index.name = "Genomes"

    # Get predictions
    print("Calculating probabilities for each domain", file=sys.stderr)
    Y_hat = df_weights/df_weights.sum(axis=1).values.reshape((-1,1))
    y_hat = Y_hat.idxmax(axis=1)

    # Add multiindex
    df_weights.columns = df_weights.columns.map(lambda x: ("Weights", x))
    Y_hat.columns = Y_hat.columns.map(lambda x: ("Probabilities", x))

    # Output
    path_predictions = os.path.join(opts.output_directory, "predictions.tsv")
    print("Creating prediction output file:", path_predictions, file=sys.stderr)
    df_output = pd.concat([
        y_hat.to_frame(("Prediction", "Domain")),
        Y_hat,
        df_weights,
    ], axis=1).sort_index()
    df_output.index.name = "Genomes"

    df_output.to_csv(path_predictions, sep="\t")

    # Prokaryotes
    path_prokaryota = os.path.join(opts.output_directory, "prokaryota.list")
    prokaryota = y_hat[y_hat.map(lambda y: y in {"Archaea", "Bacteria"})].index
    print("Creating Prokaryota list file [N = {} genomes]:".format(len(prokaryota)), path_prokaryota, file=sys.stderr)
    with open(path_prokaryota, "w") as f:
        for id_mag in prokaryota:
            print(id_mag, file=f)

    # Eukaryotes
    path_eukaryota = os.path.join(opts.output_directory, "eukaryota.list")
    eukaryota = y_hat[y_hat.map(lambda y: y in {"Eukaryota"})].index
    print("Creating Eukaryota list file [N = {} genomes]:".format(len(eukaryota)), path_eukaryota, file=sys.stderr)
    with open(path_eukaryota, "w") as f:
        for id_mag in eukaryota:
            print(id_mag, file=f)

    # Viruses
    path_viruses = os.path.join(opts.output_directory, "viruses.list")
    viruses = y_hat[y_hat.map(lambda y: y in {"Viruses"})].index
    print("Creating Viruses list file [N = {} genomes]:".format(len(viruses)), path_viruses, file=sys.stderr)
    with open(path_viruses, "w") as f:
        for id_mag in viruses:
            print(id_mag, file=f)

    # Heatmap
    if not opts.no_heatmap:
        path_heatmap = os.path.join(opts.output_directory, "heatmap.pdf")
        print("Creating heatmap of prediction probabilities:", path_heatmap, file=sys.stderr)

        import matplotlib.pyplot as plt 
        import seaborn as sns 

        df_heatmap = df_output["Probabilities"]
        if opts.figsize == "infer":
            n, m = df_heatmap.shape
            opts.figsize = (m*1.618, n*0.5)
        else:
            opts.figsize = opts.figsize.strip().split(",")
            opts.figsize = (float(opts.figsize[0]), float(opts.figsize([1])))

        with plt.style.context(opts.style):
            fig, ax = plt.subplots(figsize=opts.figsize)
            sns.heatmap(df_heatmap, vmin=0, vmax=1, cmap=opts.colormap, annot=True, ax=ax, edgecolor="white", linewidth=1, cbar_kws={"label":"Probability"})
            ax.set_xticklabels(ax.get_xticklabels(), fontsize=15)
            ax.set_yticklabels(ax.get_yticklabels(), fontsize=12)
            ax.set_xlabel("Domain", fontsize=15, fontweight="bold")
            ax.set_ylabel("Genomes", fontsize=15, fontweight="bold")
            if opts.title:
                ax.set_title(opts.title, fontsize=15, fontweight="bold")
            fig.savefig(path_heatmap, format="pdf", bbox_inches="tight", dpi=300)

if __name__ == "__main__":
    main()
