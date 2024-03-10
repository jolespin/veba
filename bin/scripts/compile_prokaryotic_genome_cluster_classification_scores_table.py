#!/usr/bin/env python
import sys, os, glob, argparse, warnings
from collections import OrderedDict
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.9.27"

def main(argv=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <gtdbtk_results.tsv> -c <clusters> -o <output>".format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    
    # I/O
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-i","--gtdbtk_results", type=str, required=True, help = "path/to/gtdbtk_results.tsv for when MASH_DB is used. [Required]")
    parser_io.add_argument("-c","--clusters", type=str, required=True, help = "path/to/clusters.tsv where [id_genome]<tab>[id_genome-cluster], No header. [Required]")
    parser_io.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")
    parser_io.add_argument("--prioritize", type=str,  default="msa", help = "Prioritize {msa, ani} [Default: msa]")
    parser_io.add_argument("--fill_missing_weight", type=float,  help = "Fill missing weight between [0, 100.0].  [Default is to throw error if value is missing]")
    parser_io.add_argument("--header", action="store_true", help = "Include header")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Checks
    assert opts.prioritize in {"msa", "ani"}, "--prioritize must be either msa or ani"
    if opts.fill_missing_weight:
        assert 0 <= opts.fill_missing_weight <= 100.0, "--fill_missing_weight must be between 0 and 100"

    # Output
    if opts.output == "stdout":
        opts.output = sys.stdout 
    
    # Inputs
    df_gtdbtk_results = pd.read_csv(opts.gtdbtk_results, sep="\t", index_col=0)
    genome_to_genomecluster = pd.read_csv(opts.clusters, sep="\t", index_col=0, header=None).iloc[:,0]

    # Check overlap between genomes
    A = set(df_gtdbtk_results.index)
    B = set(genome_to_genomecluster.index)
    C = A & B
    D = (A | B) - C
    assert A <= B, "Not all genomes from --gtdbtk_results are in --clusters.\n * Genomes only in --gtdbtk_results: {}\n".format(A - B)
    genome_to_genomecluster = genome_to_genomecluster.loc[df_gtdbtk_results.index]

    # Get classifications and weight for genomes
    genome_to_weight = dict()
    genome_to_classification = dict()
    for id_genome, row in df_gtdbtk_results.iterrows():
        classification, ani, msa = row[["classification", "fastani_ani", "msa_percent"]]
        assert pd.notnull(classification), "Missing classification for `{}`".format(id_genome)
        genome_to_classification[id_genome] = classification

        ani_not_null = pd.notnull(ani)
        msa_not_null = pd.notnull(msa)
        if all([ani_not_null, msa_not_null]):
            warnings.warn(f"`{id_genome}` has values for both `fastani_ani` and `msa_percent`. Prioritizing by `{opts.prioritize}`")
            if opts.prioritize == "msa":
                genome_to_weight[id_genome] = msa
            if opts.prioritize == "ani":
                genome_to_weight[id_genome] = ani 
        else:
            if msa_not_null:
                genome_to_weight[id_genome] = msa
            if ani_not_null:
                genome_to_weight[id_genome] = ani
    genome_to_weight = pd.Series(genome_to_weight)
    genome_to_classification = pd.Series(genome_to_classification)

    # Missing weights
    null_weights = genome_to_weight.isnull()
    if null_weights.sum() >= 1:
        warnings.warn("Missing weights for the following genomes: {}".format(null_weights.index[null_weights].tolist()))
        if not opts.fill_missing_weight:
            raise Exception("Please ensure all genomes have a weight attribute for msa_percent or fastani_ani.  The alternative is to fill missing values with --fill_missing_weight")
        else:
            genome_to_weight[null_weights] = opts.fill_missing_weight
    
    df_output = genome_to_genomecluster.to_frame("id_genome-cluster")
    df_output["classification"] = genome_to_classification
    df_output["weight"] = genome_to_weight
    df_output.index.name = "id_genome"

    df_output.to_csv(opts.output, sep="\t", header=bool(opts.header))



    
if __name__ == "__main__":
    main()
    
                

