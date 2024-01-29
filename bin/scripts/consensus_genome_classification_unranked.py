#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse
from collections import OrderedDict, defaultdict
import pandas as pd
import numpy as np
import taxopy

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.12.30"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <taxonkit-lineage> -o <formatted-lineage>".format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", default="stdin", type=str, help = "path/to/mag_to_classification.tsv  [id_mag]<tab>[id_slc]<tab>[classification]|OPTIONAL:<tab>[weight]| [Default: stdin]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "Output table with consensus classification [Default: stdout]")
    parser.add_argument("-t","--threshold", default=0.5, type=float, help = "taxopy fraction of classifications for consensus for weighted majority vote [Default: 0.5]")
    parser.add_argument("-d", "--delimiter", type=str, default=";", help = "Taxonomic delimiter [Default: ; ]")
    parser.add_argument("-x", "--blacklist", type=str, default="unknown,uncharacterized,unclassified", help = "Black list labels.  Comma separated. [Default: unknown,uncharacterized,unclassified ]")
    parser.add_argument("--unclassified_label", type=str, default="Unclassified", help = "Label to use for unclassified taxa [Default: Unclassified]")
    parser.add_argument("--unclassified_taxonid", type=int, default=-1, help = "Integer to use for unclassified taxon ids [Default: -1]")
    # parser.add_argument("--remove_missing_classifications", action="store_true", help = "Remove all classifications and weights that are  null.  For viruses this could cause an error if this isn't selected.")
    parser.add_argument("--veba_database", type=str,  help=f"VEBA database location.  [Default: $VEBA_DATABASE environment variable]")
    parser.add_argument("--verbose", action="store_true", help = "Verbose")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # I/O
    if opts.input == "stdin":
        opts.input = sys.stdin

    if opts.output == "stdout":
        opts.output = sys.stdout

    if opts.veba_database is None:
        assert "VEBA_DATABASE" in os.environ, "Please set the following environment variable 'export VEBA_DATABASE=/path/to/veba_database' or provide path to --veba_database"
        opts.veba_database = os.environ["VEBA_DATABASE"]

    # Blacklist labels
    blacklist = set()
    if opts.blacklist:
        for label in opts.blacklist.strip().split(","):
            blacklist.add(label)

    # Threshold 
    assert 0 < opts.threshold <= 1, "--threshold must be (0,1])"

    # Classifications
    df_input = pd.read_csv(opts.input, sep="\t", index_col=0, header=None)
    mag_to_slc = df_input.iloc[:,0]
    mag_to_classification = df_input.iloc[:,1]

    include_weights = False
    if df_input.shape[1] == 3:
        mag_to_weight = df_input.iloc[:,2]
        include_weights = True

    assert not np.any(df_input.isnull()), "Cannot have any missing values in --input"

    # TaxDB
    taxdb = taxopy.TaxDb(
        nodes_dmp=os.path.join(opts.veba_database, "Classify", "NCBITaxonomy", "nodes.dmp"), 
        names_dmp=os.path.join(opts.veba_database, "Classify", "NCBITaxonomy", "names.dmp"),
        merged_dmp=os.path.join(opts.veba_database, "Classify", "NCBITaxonomy", "merged.dmp"),
    )

    # Organize SLCs
    slc_to_mags = defaultdict(list)
    slc_to_classifications = defaultdict(list)
    for id_mag, row in df_input.iloc[:,:2].iterrows():
        id_slc = row.iloc[0]
        classification = row.iloc[1]
        slc_to_mags[id_slc].append(id_mag)
        slc_to_classifications[id_slc].append(classification)
    slc_to_classifications = pd.Series(slc_to_classifications)
    
    # Get taxonomy ids
    slc_to_classification = dict()
    slc_to_homogeneity = dict()
    slc_to_taxonid = dict()
    slc_to_numberofgenomes = pd.Series(slc_to_mags).map(len)
    slc_to_weights = dict()

    for id_slc, classifications in slc_to_classifications.items():

        taxa = list()
        weights = None
        if include_weights:
            weights = list()
        for i, classification in enumerate(classifications):
            id_mag = slc_to_mags[id_slc][i]
            if include_weights:
                w = mag_to_weight[id_mag]
            if classification in blacklist:
                if opts.verbose:
                    print("[Skipping Blacklist] [MAG={}] [SLC={}] [Taxonomy={}]".format(id_mag, id_slc, classification), file=sys.stderr)
            else:
                for field in classification.split(opts.delimiter)[::-1]:
                    results = taxopy.taxid_from_name(field, taxdb)
                    if len(results) == 0:
                        if opts.verbose:
                            print("[No Results] [MAG={}] [SLC={}] [Taxonomy={}] [Field={}]".format(id_mag, id_slc, classification, field), file=sys.stderr)
                    else:
                        taxon = taxopy.Taxon(results[0], taxdb)
                        taxa.append(taxon)
                        if include_weights:
                            weights.append(w)
                        if opts.verbose:
                            print("[Match] [MAG={}] [SLC={}] [Taxonomy={}] [Field={}] [TaxID={}]".format(id_mag, id_slc, classification, field, taxon._taxid), file=sys.stderr)
                        break
        if len(taxa) == 0:
            # if opts.verbose:
            #     print("[Unclassified] [MAG={}] [SLC={}] [Taxonomy={}] [Field={}]".format(id_mag, id_slc, classification, field), file=sys.stderr)
            slc_to_classification[id_slc] = opts.unclassified_label
            slc_to_homogeneity[id_slc] = 1.0
            slc_to_taxonid[id_slc] = opts.unclassified_taxonid
        else:
            if len(taxa) == 1:
                taxon = taxa[0]
                slc_to_classification[id_slc] = taxon
                slc_to_homogeneity[id_slc] = 1.0
                slc_to_taxonid[id_slc] = taxon._taxid
            else:
                majority_vote = taxopy.find_majority_vote(taxa, taxdb, fraction=opts.threshold, weights=weights)
                results = taxopy.taxid_from_name(majority_vote.name, taxdb)
                taxon = taxopy.Taxon(results[0], taxdb)
                slc_to_classification[id_slc] = taxon
                slc_to_homogeneity[id_slc] = majority_vote.agreement
                slc_to_taxonid[id_slc] = taxon._taxid
        if include_weights:
            slc_to_weights[id_slc] = weights 
        else:
            slc_to_weights[id_slc] = np.nan

    df_output = pd.DataFrame(
        OrderedDict([ 
            ("consensus_classification", pd.Series(slc_to_classification)),
            ("homogeneity", pd.Series(slc_to_homogeneity)),
            ("number_of_unique_classification", slc_to_classifications.map(lambda x: len(set(x) - {"Unclassified"}))),
            ("number_of_genomes", slc_to_numberofgenomes),
            ("genomes", slc_to_mags),
            ("classifications", slc_to_classifications),
            ("weights", slc_to_weights),
            ("consensus_taxon_id", pd.Series(slc_to_taxonid)),

        ])
    ).sort_index()
    df_output.index.name = "id_genome_cluster"
    df_output.to_csv(opts.output, sep="\t")
            
if __name__ == "__main__":
    main()
