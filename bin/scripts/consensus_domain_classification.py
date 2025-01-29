#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse
import pandas as pd
from scipy.special import softmax

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.12.27"

def main(argv=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <scaffolds_to_bins.tsv> -t <tiara_predictions.tsv> -o <output_directory>".format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-i","--scaffolds_to_bins", type=str, required=True, help = "path/to/scaffolds_to_bins.tsv, [id_scaffold]<tab>[id_bin], No header.")
    parser.add_argument("-t","--tiara_predictions", type=str, required=True, help = "Tab-seperated table as output from Tiara using --probabilities")
    parser.add_argument("-o","--output_directory", type=str, default="consensus_domain_classification_output", help = "Output directory [Default: consensus_domain_classification_output]")
    parser.add_argument("-l", "--logit_transform", type=str, default="softmax", help = "Transformation: {softmax, tss} [Default: softmax]")
    parser.add_argument("-m", "--merge_prokaryota", action="store_true", help = "Merge Archaea + Bacteria = Prokaryota")

    # Options
    opts = parser.parse_args(argv)

    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    os.makedirs(opts.output_directory, exist_ok=True)
    assert opts.logit_transform in {"softmax", "tss"}, "--logit_transform must be either softmax or tss"
    # Get scaffolds to bins
    scaffold_to_bin = pd.read_csv(opts.scaffolds_to_bins, sep="\t", index_col=0, header=None).iloc[:,0]
    scaffold_to_bin.index = scaffold_to_bin.index.map(lambda x: x.split(" ")[0])
    scaffold_to_bin.name = "id_genome"

    # Read Tiara
    df_tiara = pd.read_csv(opts.tiara_predictions, sep="\t", index_col=0)[["arc", "bac","euk", "org", "unk1"]].fillna(0)
    df_tiara.index = df_tiara.index.map(lambda x: x.split(" ")[0])
    df_tiara.columns = ["Archaea", "Bacteria", "Eukaryota", "Organelle", "Unknown"]

    if opts.merge_prokaryota:
        df_tiara = df_tiara.groupby(lambda x: {"Archaea":"Prokaryota", "Bacteria":"Prokaryota"}.get(x,x), axis=1).sum()

    # Logits
    df_logits = df_tiara.groupby(scaffold_to_bin).sum()

    # Predictions
    if opts.logit_transform == "softmax":
        df_predict_proba = pd.DataFrame(softmax(df_logits.values, axis=1), index=df_logits.index, columns=df_logits.columns)
    else:
        df_predict_proba = df_logits/df_logits.sum(axis=1).values.reshape(-1,1)

    predictions = df_predict_proba.idxmax(axis=1)

    # Output
    df_logits.to_csv(os.path.join(opts.output_directory, "prediction_logits.tsv.gz"), sep="\t")
    df_predict_proba.to_csv(os.path.join(opts.output_directory, "prediction_probabilities.tsv.gz"), sep="\t")
    predictions.to_frame().to_csv(os.path.join(opts.output_directory, "predictions.tsv"), sep="\t", header=None)
    
    for id_domain in df_predict_proba.columns:
        subset = sorted(predictions[lambda y_hat: y_hat == id_domain].index)
        with open(os.path.join(opts.output_directory, "{}.list".format(id_domain.lower())), "w") as f:
            for id_genome in subset:
                print(id_genome, file=f)




if __name__ == "__main__":
    main()
