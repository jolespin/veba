#!/usr/bin/env python
import sys, os, glob, argparse 
import pandas as pd

# DATABASE_CHECKV="/usr/local/scratch/CORE/jespinoz/db/checkv/checkv-db-v1.0/"

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.03.08"


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <checkv_results> -c <mag_to_clusters.tsv> -o <output_filepath>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--checkv_results", type=str, required=True, help = "path/to/checkv/checkv_results.filtered.tsv")
    parser.add_argument("-c","--clusters", type=str,  help = "path/to/clusters.tsv, Format: [id_mag]<tab>[id_cluster], No header [Optional]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")
    parser.add_argument("--checkv_database", type=str, default=None, required=True, help="CheckV | path/to/checkv_database (e.g. --arg 1 )")
    parser.add_argument("--header", type=int, default=1, help="Include header in output {0=No, 1=Yes) [Default: 1]")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # I/O

    # Load CheckV results
    df_checkv_results = pd.read_csv(opts.checkv_results, sep="\t", index_col=0)

    # Load CheckV database
    df_checkvdb_circular = pd.read_csv(os.path.join(opts.checkv_database, "genome_db", "checkv_circular.tsv"), sep="\t", index_col=0)
    df_checkvdb_genbank = pd.read_csv(os.path.join(opts.checkv_database, "genome_db", "checkv_genbank.tsv"), sep="\t", index_col=0)

    # Output
    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Get best hit and score
    mag_to_aai_top_hits = df_checkv_results["aai_top_hit"].dropna()
    mag_to_aai_score = 100 - df_checkv_results["aai_error"][mag_to_aai_top_hits.index]

    # Get genbank info
    mag_to_name__gb = mag_to_aai_top_hits[mag_to_aai_top_hits.map(lambda x: not x.startswith("DTR"))].map(lambda x: df_checkvdb_genbank.loc[x, "ncbi_name"])
    mag_to_vog__gb = mag_to_aai_top_hits[mag_to_aai_top_hits.map(lambda x: not x.startswith("DTR"))].map(lambda x: df_checkvdb_genbank.loc[x, "vog_clade"])
    mag_to_lineage__gb = mag_to_aai_top_hits[mag_to_aai_top_hits.map(lambda x: not x.startswith("DTR"))].map(lambda x: df_checkvdb_genbank.loc[x, "lineage"])

    # MAG -> Classification
    mag_to_classification = pd.concat([
        mag_to_aai_top_hits[mag_to_aai_top_hits.map(lambda x: x.startswith("DTR"))].map(lambda x: df_checkvdb_circular.loc[x, "lineage"]),
        mag_to_vog__gb,
    ])

    # MAG -> Source
    mag_to_source = pd.concat([
        mag_to_aai_top_hits[mag_to_aai_top_hits.map(lambda x: x.startswith("DTR"))].map(lambda x: df_checkvdb_circular.loc[x, "habitat"]),
        mag_to_aai_top_hits[mag_to_aai_top_hits.map(lambda x: not x.startswith("DTR"))].map(lambda x: df_checkvdb_genbank.loc[x, "isolation"]),

    ])

    # Output
    df_mag_classifications = pd.concat([
        mag_to_classification.to_frame("vog_clade"),
        mag_to_aai_score.to_frame("aai_score"),
        mag_to_aai_top_hits.to_frame("aai_top_hit"),
        mag_to_source.to_frame("source"),
        mag_to_name__gb.to_frame("ncbi_name"),
        mag_to_lineage__gb.to_frame("lineage"),
    ],
    axis=1)
    df_mag_classifications.index.name = "id_genome"

    # Add clusters if provided
    if opts.clusters:
        mag_to_cluster = pd.read_csv(opts.clusters, sep="\t", index_col=0, header=None).iloc[:,0]
        df_mag_classifications.insert(loc=0, column="id_genome_cluster", value=mag_to_cluster[df_mag_classifications.index])

    # Output
    df_mag_classifications.to_csv(opts.output, sep="\t", header=bool(opts.header))



if __name__ == "__main__":
    main()
    
                

