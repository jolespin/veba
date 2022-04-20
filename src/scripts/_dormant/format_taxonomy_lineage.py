#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, warnings
import pandas as pd
import numpy as np
from collections import OrderedDict

pd.options.display.max_colwidth = 100
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.02.02"

# Defaults
RANK_TO_PREFIX="superkingdom:d__,phylum:p__,class:c__,order:o__,family:f__,genus:g__,species:s__"

# Read taxonkit lineage
def read_taxonkit_lineage(path, ranks="infer"):
    df = pd.read_csv(path, sep="\t", index_col=0, header=None)
    df.index.name = "id_taxon"
    if ranks == "infer":
        ranks = False
        if df.shape[1] == 3:
            ranks = True
    if ranks:
        df.columns = ["lineage", "names", "ranks"]
    else:
        df.columns = ["lineage", "names"]
    return df

# Format taxonomy lineage
def format_taxonomy_lineage(
    lineage:pd.Series, 
    ranks:pd.Series,
    rank_to_prefix:OrderedDict,
    delimiter=";",
    partial_ok=False,
    ):
    # Check inputs
    assert set(lineage.index) == set(ranks.index), "`lineage` and `ranks` must have same keys in index"
    ranks = ranks[lineage.index]
    assert np.all(lineage.map(lambda x: len(x.split(delimiter))) == ranks.map(lambda x: len(x.split(delimiter)))), "Fields delimitted by `delimiter={}` do not match up between `lineage` and `ranks`".format(delimiter)
    prefixes = np.asarray(list(rank_to_prefix.values()))
    
    # Format lineage
    taxonid_to_lineage = dict()
    for id_taxon, (t, r) in tqdm(pd.DataFrame([lineage, ranks]).T.iterrows(), "Formatting taxonomy lineages", unit=" NCBItaxids"):
        
        # Split by delimiter
        rank_fields = np.asarray(r.split(delimiter))
        taxon_fields = np.asarray(t.split(delimiter))
        
        # Get positions of query ranks
        rank_indices = list()
        for i, rank in enumerate(rank_fields):
            if rank in rank_to_prefix:
                rank_indices.append(i)
                
        # Allow or don't allow for only a subset of ranks 
        all_ranks_included = len(rank_indices) == len(prefixes)
        msg = "[NCBItaxid: {}] Missing taxonomic ranks: {}".format(id_taxon, set(rank_to_prefix.keys()) - set(rank_fields))
        if not partial_ok:
            assert all_ranks_included, msg
        else:
            
            if not all_ranks_included:
                warnings.warn(msg)
        
        # Merge prefixes with query ranks
        taxonid_to_lineage[id_taxon] = delimiter.join(map(lambda item: item[0] + item[1], zip(prefixes, taxon_fields[rank_indices])))
    
    return pd.Series(taxonid_to_lineage, name="id_taxon")[lineage.index]


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
    parser.add_argument("-i","--taxonkit_lineage", default="stdin", type=str, help = "path/to/taxonkit_lineage.tsv with -R option [id_taxon]<tab>[taxon_level_identifiers]<tab>[taxon_level_ranks] [Default: stdin]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "Output table [id_taxon]<tab>[lineage] [Default: stdout]")
    parser.add_argument("-r","--rank_to_prefix", type=str, default=RANK_TO_PREFIX, help = "Rank to prefix with the following format with ranks from NCBI taxonomy: 'rank_1:prefix_a,rank_2:prefix_b'\n[Default: {}".format(RANK_TO_PREFIX))
    parser.add_argument("-d", "--delimiter", type=str, default=";", help = "Taxonomic delimiter [Default: ;]")
    parser.add_argument( "--partial_ok", action="store_true", help = "Allow for only a subset of ranks from --rank_to_prefix")
    parser.add_argument( "--include_header", action="store_true", help = "Include header on output [id_taxon]<tab>[lineage]")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # I/O
    if opts.taxonkit_lineage == "stdin":
        opts.taxonkit_lineage = sys.stdin

    if opts.output == "stdout":
        opts.output = sys.stdout

    # Format rank to lineage
    opts.rank_to_prefix = OrderedDict(map(lambda item: item.split(":"), opts.rank_to_prefix.strip().split(",")))

    # Read taxonkit lineage table
    df_taxonkit_lineage = read_taxonkit_lineage(opts.taxonkit_lineage, ranks=True)

    # Format lineage 
    taxonid_to_lineage = format_taxonomy_lineage(
        lineage=df_taxonkit_lineage["lineage"], 
        ranks=df_taxonkit_lineage["ranks"],
        rank_to_prefix=opts.rank_to_prefix,
        delimiter=opts.delimiter,
        partial_ok=bool(opts.partial_ok),
    )
    taxonid_to_lineage.to_frame("lineage").to_csv(opts.output, sep="\t", header=bool(opts.include_header))

if __name__ == "__main__":
    main()
