#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse
from collections import OrderedDict, defaultdict
import pandas as pd
import numpy as np


pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.11.3"

# RANK_TO_PREFIX="superkingdom:d__,phylum:p__,class:c__,order:o__,family:f__,genus:g__,species:s__"

RANK_PREFIXES="d__,p__,c__,o__,f__,g__,s__"

# Fill empty taxonomic levels for consensus classification
def fill_lower_taxonomy_levels(
    classifications:pd.Series,
    rank_prefixes:list,
    delimiter:str=";",
    ):
    
    rank_prefixes = list(rank_prefixes)
    number_of_taxonomic_levels = len(rank_prefixes)
    classifications_ = dict()
    for id_genome, classification in pd.Series(classifications).items():
        taxonomy = classification.split(delimiter)
        classifications_[id_genome] = delimiter.join(taxonomy + rank_prefixes[len(taxonomy):])
    return pd.Series(classifications_)[classifications.index]

# Get consensus classification
def get_consensus_classification(
    classification:pd.Series, 
    classification_weights:pd.Series,
    genome_to_genomecluster:pd.Series, 
    rank_prefixes:list,
    number_of_taxonomic_levels="infer",
    delimiter=";",
    leniency:float=1.382,
    ):
    # Assertions
    assert np.all(classification.notnull())
    assert np.all(classification_weights.notnull())
    assert np.all(genome_to_genomecluster.notnull())
    
    # Set and index overlap
    a = set(classification.index)
    b = set(classification_weights.index)
    c = set(genome_to_genomecluster.index)
    assert a == b, "`classification` and `classification_weights` must  have the same keys in the index"
    assert a <= c, "`classification` and `classification_weights` must be a subset (or equal) to the keys in `genome_to_genomecluster` index"
    index_genomes = pd.Index(sorted(a & b & c ))
    classification = classification[index_genomes]
    classification_weights = classification_weights[index_genomes]
    genome_to_genomecluster = genome_to_genomecluster[index_genomes]
    
    # Taxonomic levels
    taxonomic_levels = classification.map(lambda x: x.count(delimiter)).unique()
    if len(taxonomic_levels):
        assert len(taxonomic_levels) == 1, "Taxonomic levels in `classification` should all have the same number of delimiters" #! Might need to change this to allow for missing taxonomic levels
    else:
        number_of_taxonomic_levels = 1

    if number_of_taxonomic_levels == "infer":
        number_of_taxonomic_levels = taxonomic_levels[0] + 1

    # Scaling factors
    scaling_factors = np.arange(1, number_of_taxonomic_levels + 1) # d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Actinomycetales;f__Dermabacteraceae;g__Brachybacterium
    scaling_factors = np.power(scaling_factors, leniency)
    
    # Get container for scores [SLC -> Taxonomy -> Score]
    #
    # For example the following MAG: 
    # CLASSIFICATION=d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Corynebacterium;s__Corynebacterium aurimucosum_E
    # MSA_PERCENT=80.0
    #
    # Would be stored and appended for it's corresponding SLC:
    # d__Bacteria += 80.0
    # d__Bacteria;p__Actinobacteriota += 80.0
    # d__Bacteria;p__Actinobacteriota;c__Actinomycetia += 80.0
    # d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales += 80.0
    # d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae += 80.0
    # d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Corynebacterium += 80.0
    # d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Corynebacterium;s__Corynebacterium aurimucosum_E += 80.0
    genomecluster_taxa_scores = defaultdict(lambda: defaultdict(float))
    
    # Iterate through MAG, classification, and score
    df = pd.concat([genome_to_genomecluster.to_frame("id"), classification.to_frame("classification"), classification_weights.to_frame("weight")], axis=1)
    genomecluster_to_genomes = defaultdict(list)
    for id_genome, (id_genome_cluster, classification, w) in df.iterrows():
        genomecluster_to_genomes[id_genome_cluster].append(id_genome)
        # Split the taxonomy classification by levels
        levels = classification.split(delimiter)
        # Remove the empty taxonomy levels (e.g., g__Corynebacterium;s__ --> g__Corynebacterium)
        # levels = list(filter(lambda x:x not in rank_prefixes, levels))
        number_of_query_levels = len(levels)
        # Iterate through each level, scale score by the leniency weights, and add to running sum
        for i in range(1, number_of_query_levels + 1):
            taxon_at_level = levels[i-1]
            taxon_level_is_missing = taxon_at_level in rank_prefixes
            if taxon_level_is_missing:
                weighted_score = 0.0
                print("`{}` is missing taxonomic level `{}`".format(id_genome, taxon_at_level), file=sys.stderr)

            else:
                weighted_score = float(w) * scaling_factors[i-1]
            genomecluster_taxa_scores[id_genome_cluster][tuple(levels[:i])] += weighted_score
    genomecluster_to_genomes = pd.Series(genomecluster_to_genomes)

    # Build datafarme
    genomecluster_taxa_scores = pd.Series(genomecluster_taxa_scores)
    df_consensus_classification = pd.DataFrame(genomecluster_taxa_scores.map(lambda taxa_scores: sorted(taxa_scores.items(), key=lambda x:(x[1], len(x[0])), reverse=True)[0]).to_dict(), index=["consensus_classification", "score"]).T
    df_consensus_classification["consensus_classification"] = df_consensus_classification["consensus_classification"].map(";".join)
    df_consensus_classification["number_of_unique_classifications"] = df["classification"].groupby(genome_to_genomecluster).apply(lambda x: len(set(x)))
    df_consensus_classification["number_of_components"] = genomecluster_to_genomes.map(len) #df["classification"].groupby(genome_to_genomecluster).apply(len)
    df_consensus_classification["components"] = genomecluster_to_genomes
    df_consensus_classification["classifications"] = df["classification"].groupby(genome_to_genomecluster).apply(lambda x: list(x))
    df_consensus_classification["weights"] = df["weight"].groupby(genome_to_genomecluster).apply(lambda x: list(x))
    df_consensus_classification.index.name = "id"
    
    # Homogeneity
    slc_taxa_homogeneity = defaultdict(lambda: defaultdict(float))
    for id_genome_cluster, (classifications, weights) in df_consensus_classification[["classifications", "weights"]].iterrows():
        for (c, w) in zip(classifications, weights):
            slc_taxa_homogeneity[id_genome_cluster][c] += w
    df_consensus_classification["homogeneity"] = pd.DataFrame(slc_taxa_homogeneity).T.apply(lambda x: np.nanmax(x)/np.nansum(x), axis=1)
        
    fields = [
        "consensus_classification", 
        "homogeneity", 
        "number_of_unique_classifications",
        "number_of_components",
        "components",
        "classifications", 
        "weights", 
        "score",
    ]
    return df_consensus_classification.loc[:,fields]



def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <input> -o <output>".format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", default="stdin", type=str, help = "path/to/genome_to_classification.tsv  [id_genome]<tab>[id_genome_cluster]<tab>[classification]<tab>[weight]; No header. [Default: stdin]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "Output table with consensus classification [Default: stdout]")
    parser.add_argument("-l","--leniency", default=1.382, type=float, help = "Leniency parameter. Lower value means more conservative weighting. A value of 1 indiciates no weight bias. A value greater than 1 puts higher weight on higher level taxonomic assignments. A value less than 1 puts lower weights on higher level taxonomic assignments.  [Default: 1.382]")
    parser.add_argument("-r", "--rank_prefixes", type=str, default=RANK_PREFIXES, help = "Rank prefixes separated by , delimiter'\n[Default: {}]".format(RANK_PREFIXES))
    parser.add_argument("-d", "--delimiter", type=str, default=";", help = "Taxonomic delimiter [Default: ; ]")
    parser.add_argument("-s", "--simple", action="store_true", help = "Simple classification that does not use lineage information from --rank_prefixes")
    # parser.add_argument("--assert_resolved_taxonomy", action="store_true", help = "Do not allow missing taxonomic levels. (e.g., d__Eukaryota;p__;c__Pelagophyceae;o__Pelagomonadales;f__;g__Aureococcus;s__Aureococcus anophagefferens is missing phylum)")
    parser.add_argument("--remove_missing_classifications", action="store_true", help = "Remove all classifications and weights that are  null.  For viruses this could cause an error if this isn't selected.")
    parser.add_argument("-u", "--unclassified_label", default="Unclassified", type=str, help = "Unclassified label [Default: Unclassified]")
    parser.add_argument("-w", "--unclassified_weight", default=100.0,type=float, help = "Unclassified label weight [Default: 100.0]")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # I/O
    if opts.input == "stdin":
        opts.input = sys.stdin

    if opts.output == "stdout":
        opts.output = sys.stdout

    # Leniency 
    assert opts.leniency > 0, "--leniency must be > 0"
    # Format rank to lineage
    opts.rank_prefixes = opts.rank_prefixes.strip().split(",")

    # Classifications
    df_input = pd.read_csv(opts.input, sep="\t", index_col=0, header=None)
    genome_to_genomecluster = df_input.iloc[:,0]
    genome_to_classification = df_input.iloc[:,1].reindex(genome_to_genomecluster.index)
    genome_to_weights = df_input.iloc[:,2].reindex(genome_to_genomecluster.index)
    if opts.remove_missing_classifications:
        genome_to_weights = genome_to_weights.dropna()
        genome_to_classification = genome_to_classification[genome_to_weights.index]
    else:
        mask = genome_to_weights.isnull()
        genome_to_classification[mask] = ";".join(map(lambda x: f"{x}__{opts.unclassified_label}", opts.rank_prefixes))
        genome_to_weights[mask] = opts.unclassified_weight
        
        
    # Consensus classification
    df_consensus_classification = get_consensus_classification(
        classification=genome_to_classification, 
        classification_weights=genome_to_weights,
        genome_to_genomecluster=genome_to_genomecluster,
        rank_prefixes=opts.rank_prefixes,
        number_of_taxonomic_levels="infer",
        delimiter=opts.delimiter,
        leniency=opts.leniency,
    )

    if not opts.simple:
        # Fill empty taxonomy levels
        df_consensus_classification["consensus_classification"] = fill_lower_taxonomy_levels(
            classifications=df_consensus_classification["consensus_classification"],
            rank_prefixes=opts.rank_prefixes,
            delimiter=opts.delimiter,
        )

    df_consensus_classification.to_csv(opts.output, sep="\t")

if __name__ == "__main__":
    main()
