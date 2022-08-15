#!/usr/bin/env python
import sys, os, glob, argparse, gzip
from collections import defaultdict
import pandas as pd
from soothsayer_utils import get_file_object, read_hmmer, assert_acceptable_arguments

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.06.16"

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
    parser.add_argument("--identifiers_only", action="store_true", help="Output query identifiers only")
    # parser.add_argument("--use_hmmsearch_header", action="store_true", help="Use the original hmmsearch header instead of the multiindex version")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert_acceptable_arguments(opts.hmm_marker_field, {"accession", "name"})

    if opts.output == "stdout":
        f_out = sys.stdout
    else:
        f_out = open(opts.output, "w")

    # Read HMMER

    # header_lines=None
    # if opts.use_hmmsearch_header:
    #     header_lines = list()
    #     with get_file_object(opts.hmmsearch_tblout, mode="read") as f:
    #         for line in f:
    #             line = line.strip()
    #             if line.startswith("#"):
    #                 header_lines.append(line)
    #             else:
    #                 break

    df_tblout = read_hmmer(opts.hmmsearch_tblout, program="hmmsearch", format="tblout", add_header_as_index=False)
    if not opts.scores_cutoff:
        df_tblout.to_csv(f_out, sep="\t", index=None)
    else:
        hmm_to_cutoff = pd.read_csv(opts.scores_cutoff, sep="\t", index_col=0, header=None).iloc[:,0]

        hmmsearch_filtered = list()
        for i, row in df_tblout.iterrows():
            id_hmm = row[("identifier", {"name":"query_name", "accession":"query_accession"}[opts.hmm_marker_field])]
            score = float(row[("full_sequence", "score")])
            if score >= hmm_to_cutoff[id_hmm]:
                hmmsearch_filtered.append(row)
        df_tblout_filtered = pd.DataFrame(hmmsearch_filtered)

        if opts.identifiers_only:
            for id_hmm in sorted(set(df_tblout_filtered[("identifier", "target_name")])):
                print(id_hmm, file=f_out)
        else:
            df_tblout_filtered.to_csv(f_out, sep="\t", index=None)

    if f_out != sys.stdout:
        f_out.close()
    
if __name__ == "__main__":
    main()
    
                

