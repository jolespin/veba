#!/usr/bin/env python
import sys, os, glob, argparse 
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.02.23"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <scaffold_to_mag> -c <mag_to_clusters.tsv> -o <output_filepath>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--scaffolds_to_bins", type=str, required=True,  help = "path/to/scaffolds_to_bins.tsv, Format: [id_scaffold]<tab>[id_bin], No header")
    parser.add_argument("-c","--clusters", type=str, required=True, help = "path/to/clusters.tsv, Format: [id_mag]<tab>[id_cluster], No header")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")
    parser.add_argument("--header", action="store_true", help="Specify if header should be in output")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Organize identifiers
    scaffold_to_mag = pd.read_csv(opts.scaffolds_to_bins, sep="\t", index_col=0, header=None).iloc[:,0]
    mag_to_cluster = pd.read_csv(opts.clusters, sep="\t", index_col=0, header=None).iloc[:,0]
    scaffold_to_cluster = scaffold_to_mag.map(lambda id_mag: mag_to_cluster[id_mag])

    # Create output table
    df_output = pd.concat( [
        scaffold_to_mag.to_frame("id_genome"),
        scaffold_to_cluster.to_frame("id_cluster"),
    ]
    , axis=1)
    df_output.index.name = "id_scaffold"

    if opts.output == "stdout":
        opts.output = sys.stdout 
    
    df_output.to_csv(opts.output, sep="\t", header=bool(opts.header))


    

if __name__ == "__main__":
    main()
    
                

