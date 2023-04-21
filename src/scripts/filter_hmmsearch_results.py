#!/usr/bin/env python
import sys, os, glob, argparse, gzip
from collections import defaultdict
import pandas as pd
from soothsayer_utils import pv, get_file_object, read_hmmer, assert_acceptable_arguments

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.4.18"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <hmmsearch_tblout> -o <proteins> -o <output_directory> -s '|--|'".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-i","--hmmsearch_tblout", type=str, help = "path/to/hmmsearch_tblout.tsv", required=True)
    parser.add_argument("-s","--scores_cutoff", type=str, help = "path/to/scores_cutoff.tsv. No header. [id_hmm]<tab>[score]")
    parser.add_argument("-o","--output", type=str,  default="stdout", help = "Output file")
    parser.add_argument("-f", "--hmm_marker_field", default="accession", type=str, help="HMM reference type (accession, name) [Default: accession]")
    parser.add_argument("--synopsis",  type=str, help="path/to/synopsis.tsv [id_protein]<tab>[markers]<tab>[scores]<tab>[e-values]")
    parser.add_argument("--region",  default="full_sequence", type=str, help="{full_sequence, best_domain} [Default: full_sequence]")
    parser.add_argument("--identifiers_only", action="store_true", help="Output query identifiers only")

    # parser.add_argument("--use_hmmsearch_header", action="store_true", help="Use the original hmmsearch header instead of the multiindex version")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert_acceptable_arguments(opts.hmm_marker_field, {"accession", "name"})
    assert_acceptable_arguments(opts.region, {"full_sequence", "best_domain"})

    if opts.output == "stdout":
        opts.output = sys.stdout

    # Read HMMER
    df_tblout = read_hmmer(opts.hmmsearch_tblout, program="hmmsearch", format="tblout", add_header_as_index=False)
    if opts.scores_cutoff:
        hmm_to_cutoff = pd.read_csv(opts.scores_cutoff, sep="\t", index_col=0, header=None).iloc[:,0]
        hmmsearch_filtered = list()
        for i, row in df_tblout.iterrows():
            id_hmm = row[("identifier", {"name":"query_name", "accession":"query_accession"}[opts.hmm_marker_field])]
            score = float(row[(opts.region, "score")])
            if score >= hmm_to_cutoff[id_hmm]:
                hmmsearch_filtered.append(row)
        df_tblout = pd.DataFrame(hmmsearch_filtered)

    if opts.identifiers_only:
        if opts.output == sys.stdout:
            for id_hmm in sorted(set(df_tblout[("identifier", "target_name")])):
                print(id_hmm, file=sys.stdout)
        else:
            with open(opts.output, "w") as f_out:
                for id_hmm in sorted(set(df_tblout[("identifier", "target_name")])):
                    print(id_hmm, file=f_out)
    else:
        df_tblout.to_csv(opts.output, sep="\t", index=None)

    if opts.synopsis:
        synopsis = defaultdict(lambda: defaultdict(list))
        for i, (id_protein, id_hmm, score, e) in pv(df_tblout.loc[:,[("identifier", "target_name"), ("identifier", "query_{}".format(opts.hmm_marker_field)), (opts.region, "score"), (opts.region, "e-value")]].iterrows(), "Compiling synopsis"):
            synopsis[id_protein]["id_hmms"].append(id_hmm)
            synopsis[id_protein]["scores"].append(float(score))
            synopsis[id_protein]["e-values"].append(float(e))
        df_synopsis = pd.DataFrame(synopsis).T
        df_synopsis.index.name = "id_protein"
        df_synopsis.to_csv(opts.synopsis, sep="\t")


    
if __name__ == "__main__":
    main()
    
                

