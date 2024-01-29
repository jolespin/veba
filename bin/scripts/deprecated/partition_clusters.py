#!/usr/bin/env python
import sys, os, glob, argparse 
from collections import OrderedDict
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.08.12"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -c <clusters> -a <proteins> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-c","--clusters", type=str, help = "path/to/clusters.tsv")
    parser.add_argument("-a","--proteins", type=str, help = "path/to/proteins.tsv")
    parser.add_argument("-o","--output_directory", type=str, default="clusters", help = "Output directory [Default: clusters]")
    # parser.add_argument("-p", "--cluster_prefix", type=str, default="SLC", help="Cluster prefix [Default: 'SLC")
    # parser.add_argument("-s", "--cluster_suffix", type=str, default="", help="Cluster suffix [Default: '")
    parser.add_argument("-x", "--extension", type=str, default="faa", help="File extension for proteins [Default: 'faa")
    parser.add_argument("--clone_singletons", action="store_true", help="Clone singletons (i.e. a dummy genome) to do paralog analysis instead of ortholog analysis")
    parser.add_argument("--clone_label", type=str, default="___clone", help="Cluster suffix [Default: '___clone")
    parser.add_argument("--copy", action="store_true", help="Copy instead of symlink")



    # parser.add_argument("--sep", type=str, default="\t",  help = "Seperator [Default: '\t'")
    # parser.add_argument("--scaffold_column_name", type=str, default="Scaffold", help="Scaffold column name [Default: Scaffold")
    # parser.add_argument("--bin_column_name", type=str, default="Bin", help="Bin column name [Default: Bin")
    # parser.add_argument("--column_order", type=str, default="scaffold,bin", help="Column order.  Specify either 'scaffold,bin' or 'bin,scaffold' [Default:scaffold,bin")
    # parser.add_argument("--bin_prefix", type=str, default="", help="Bin prefix [Default: '")
    # parser.add_argument("--header", action="store_true", help="Specify if header should be in output")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Parse
    print(" * Reading clusters table:", opts.clusters, file=sys.stderr)
    opts.clusters = pd.read_csv(opts.clusters, sep="\t", index_col=0, header=None).iloc[:,0].map(str)
    print(" * Reading proteins table:", opts.proteins, file=sys.stderr)
    opts.proteins = pd.read_csv(opts.proteins, sep="\t", index_col=0, header=None).iloc[:,0]
    assert set(opts.clusters.index) == set(opts.clusters.index), "--clusters and --proteins must have the same genomes in the first column"
    opts.proteins = opts.proteins.loc[opts.clusters.index]


    # # Add cluster prefix and suffix
    # if bool(opts.cluster_prefix):
    #     opts.clusters = opts.clusters.map(lambda x: "{}{}".format(opts.cluster_prefix, x))
    # if bool(opts.cluster_suffix):
    #     opts.clusters = opts.clusters.map(lambda x: "{}{}".format(x, opts.cluster_suffix))

    singleton_clusters = opts.clusters.value_counts()[lambda x: x == 1].index

    print(" * Making cluster subdirectories:", opts.output_directory, file=sys.stderr)
    os.makedirs(opts.output_directory, exist_ok=True)

    # Make directories
    clusters_to_proteins = OrderedDict()
    for id_cluster in opts.clusters.unique():
        path = os.path.join(opts.output_directory, id_cluster)
        os.makedirs(path, exist_ok=True)
        clusters_to_proteins[id_cluster] = path
    clusters_to_proteins = pd.Series(clusters_to_proteins) 
    clusters_to_proteins.to_frame().to_csv(os.path.join(opts.output_directory, "clusters_to_proteins.tsv"), sep="\t", header=None)
    

    
    if opts.copy:
        print(" * Copying proteins to subdirectories", file=sys.stderr)
        from shutil import copyfile
        for id_mag, src in opts.proteins.items():
            id_cluster = opts.clusters[id_mag]
            src = os.path.realpath(src)
            dst = os.path.join(opts.output_directory, id_cluster, "{}.{}".format(id_mag, opts.extension))
            copyfile(src, dst)

            # Singletons 
            conditions = [ 
                id_cluster in singleton_clusters,
                bool(opts.clone_singletons),
            ]
            if all(conditions):
                dst = os.path.join(opts.output_directory, id_cluster, "{}{}.{}".format(id_mag, opts.clone_label, opts.extension))
                copyfile(src, dst)
        
    else:
        print(" * Symlinking proteins to subdirectories", file=sys.stderr)
        for id_mag, src in opts.proteins.items():
            id_cluster = opts.clusters[id_mag]
            src = os.path.realpath(src)
            dst = os.path.join(opts.output_directory, id_cluster, "{}.{}".format(id_mag, opts.extension))
            if os.path.exists(dst):
                os.remove(dst)
            os.symlink(src, dst)

            # Singletons 
            conditions = [ 
                id_cluster in singleton_clusters,
                bool(opts.clone_singletons),
            ]
            if all(conditions):
                dst = os.path.join(opts.output_directory, id_cluster, "{}{}.{}".format(id_mag, opts.clone_label, opts.extension))
                if os.path.exists(dst):
                    os.remove(dst)
                os.symlink(src, dst)


if __name__ == "__main__":
    main()
    
                

