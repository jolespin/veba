#!/usr/bin/env python
import sys, os, glob, argparse 
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.03.24"

def main(argv=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -a <gtdbtk.ar122.summary.tsv> -b <gtdbtk.bac120.summary.tsv> -o <output.tsv>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-a","--archaea_summary", required=True, type=str, help = "path/to/gtdbtk.ar122.summary.tsv")
    parser.add_argument("-b","--bacteria_summary", required=True, type=str, help = "path/to/gtdbtk.bac120.summary.tsv")
    parser.add_argument("-c","--clusters", type=str, help = "path/to/clusters.tsv, Format: [id_mag]<tab>[id_cluster], No header. [Optional]")
    parser.add_argument("-d","--domain", type=str, default="domain", help = "Domain column label [Default: domain]")
    parser.add_argument("-o","--output", type=str,  default="stdout", help = "Output merged multiple sequence alignment [Default: stdout]")

    # Options
    opts = parser.parse_args(argv)
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Open output file
    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Read GTDB-Tk
    dataframes = list()
    if os.path.exists(opts.archaea_summary):
        df_ar122 = pd.read_csv(opts.archaea_summary, sep="\t", index_col=0)
        # Add domain
        df_ar122.insert(loc=0, column=opts.domain, value="Archaea")
        dataframes.append(df_ar122)

    if os.path.exists(opts.bacteria_summary):
        df_bac120 = pd.read_csv(opts.bacteria_summary, sep="\t", index_col=0)
        # Add domain
        df_bac120.insert(loc=0, column=opts.domain, value="Bacteria")
        dataframes.append(df_bac120)

    # Merge
    df_output = pd.concat(dataframes, axis=0)

    if opts.clusters:
        if opts.clusters != "None":
            mag_to_cluster = pd.read_csv(opts.clusters, sep="\t", index_col=0, header=None).iloc[:,0]
            gtdbtk_mags = set(df_output.index)
            cluster_mags = set(mag_to_cluster.index)
            assert gtdbtk_mags <= cluster_mags, "All genomes from --archaea_summary and --bacteria_summary should be in --clusters file:\n{}".format(gtdbtk_mags - cluster_mags)
            df_output.insert(loc=0, column="id_cluster", value=mag_to_cluster)

    # Output
    df_output.to_csv(opts.output, sep="\t")



if __name__ == "__main__":
    main()
    
                

