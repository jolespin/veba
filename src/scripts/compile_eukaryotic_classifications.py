#!/usr/bin/env python
import sys, os, glob, argparse, pickle, warnings, gzip
import pandas as pd
import numpy as np
from collections import OrderedDict
from tqdm import tqdm 


# DATABASE_EUKARYOTIC="/usr/local/scratch/CORE/jespinoz/db/veba/v1.0/Classify/Eukaryotic/"

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.12.07"


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <metaeuk_identifier_mapping> -s <scaffolds_to_bins.tsv> -c <mag_to_clusters.tsv> -o <output_filepath>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--metaeuk_identifier_mapping", type=str, required=True, help = "path/to/identifier_mapping.metaeuk.tsv")
    parser.add_argument("-s","--scaffolds_to_bins", type=str, required=True, help = "path/to/scaffolds_to_bins.tsv")
    parser.add_argument("-c","--clusters", type=str,  help = "path/to/clusters.tsv, Format: [id_mag]<tab>[id_cluster], No header [Optional]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")
    parser.add_argument("--eukaryotic_database", type=str, default=None, required=True, help="path/to/eukaryotic_database (e.g. --arg 1 )")
    # parser.add_argument("--veba_database", type=str, default=None, help=f"VEBA database location.  [Default: $VEBA_DATABASE environment variable]")
    parser.add_argument("--header", type=int, default=1, help="Include header in output {0=No, 1=Yes) [Default: 1]")
    parser.add_argument("--debug", action="store_true")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Database
    # if opts.veba_database is None:
    #     assert "VEBA_DATABASE" in os.environ, "Please set the following environment variable 'export VEBA_DATABASE=/path/to/veba_database' or provide path to --veba_database"
    # else:
    #     opts.veba_database = os.environ["VEBA_DATABASE"]
    # opts.eukaryotic_database = os.path.join(opts.veba_database, "Classify", "Microeukaryotic")

    # I/O
    # Scaffolds -> Bins
    fp = opts.scaffolds_to_bins
    print("* Reading scaffolds to bins table {}".format(fp), file=sys.stderr)
    scaffold_to_bin = pd.read_csv(fp, sep="\t", index_col=0, header=None).iloc[:,0]
    if opts.debug:
        print(fp, file=sys.stderr)
        scaffold_to_bin.head().to_csv(sys.stderr, sep="\t", header=None)
        print("\n", file=sys.stderr)
        
    # SourceID -> Taxonomy
    fp = os.path.join(opts.eukaryotic_database,"source_taxonomy.tsv.gz")
    print("* Reading source taxonomy table {}".format(fp), file=sys.stderr)
    df_source_taxonomy = pd.read_csv(fp, sep="\t", index_col=0)
    df_source_taxonomy.index = df_source_taxonomy.index.map(str)
    if opts.debug:
        print(fp, file=sys.stderr)
        df_source_taxonomy.head().to_csv(sys.stderr, sep="\t")
        print("\n", file=sys.stderr)

    # VEBA -> SourceID
    fp = os.path.join(opts.eukaryotic_database,"target_to_source.dict.pkl.gz")
    print("* Reading target to source mapping {}".format(fp), file=sys.stderr)
    with gzip.open(fp, "rb") as f:
        target_to_source = pickle.load(f)
    #target_to_source = pd.read_csv(fp, sep="\t", index_col=0, dtype=str, usecols=["id_veba", "id_source"], squeeze=True)#.iloc[:,0]
    if opts.debug:
        print(fp, file=sys.stderr)
        print(list(target_to_source.items())[:5], sep="\n", file=sys.stderr)
        print("\n", file=sys.stderr)

    # MetaEuk headers parsed
    fp = opts.metaeuk_identifier_mapping
    print("* Reading MetaEuk identifier mapping table {}".format(fp), file=sys.stderr)
    df_metaeuk = pd.read_csv(fp, sep="\t", index_col=0)
    if opts.debug:
        print(fp, file=sys.stderr)
        df_metaeuk.head().to_csv(sys.stderr, sep="\t")
        print("\n", file=sys.stderr)

    orf_to_bitscore = df_metaeuk["bitscore"].map(float)
    orf_to_scaffold = df_metaeuk["C_acc"].map(str)
    orf_to_mag = orf_to_scaffold.map(lambda id_scaffold: scaffold_to_bin[id_scaffold])

    orf_to_target = df_metaeuk["T_acc"]
    orf_to_source = orf_to_target.map(lambda id_target: target_to_source.get(id_target,np.nan))
    if np.any(pd.isnull(orf_to_source)):
        warnings.warn("The following gene - target identifiers are not in the database file: {}".format(
            os.path.join(opts.eukaryotic_database,"target_to_source.dict.pkl.gz"), 
            ),
        )
        orf_to_target[orf_to_source[orf_to_source.isnull()].index].to_frame().to_csv(sys.stderr, sep="\t", header=None)
        orf_to_source = orf_to_source.dropna()

    # Lineage
    orf_to_lineage = OrderedDict()

    missing_lineage = list()
    for id_orf, id_source in tqdm(orf_to_source.items(), desc="Retrieving lineage", unit = " ORFs"):
        if id_source in df_source_taxonomy.index:
            lineage = df_source_taxonomy.loc[id_source] #  class   order   family  genus   species
            lineage = ";".join(map(lambda items: "".join(items), zip(["c__", "o__", "f__", "g__", "s__"], lineage)))
            orf_to_lineage[id_orf] = lineage
        else:
            missing_lineage.append(id_source)
    orf_to_lineage = pd.Series(orf_to_lineage)

    if len(missing_lineage):
        warnings.warn("The following source identifiers are not in the database file: {}\n{}`".format(
            os.path.join(opts.eukaryotic_database,"source_taxonomy.tsv.gz"), 
            "\n".join(set(missing_lineage)),
            ),
        )

    # Output
    # ["id_orf", "id_mag", "bitscore", "lineage"]
    df_orf_classifications = pd.concat([
        orf_to_scaffold.to_frame("id_scaffold"),
        orf_to_mag.to_frame("id_mag"),
        orf_to_target.to_frame("id_target"),
        orf_to_source.to_frame("id_source"),
        orf_to_lineage.to_frame("lineage"),
        orf_to_bitscore.to_frame("bitscore"),
    ],
    axis=1)
    df_orf_classifications.index.name = "id_gene"


    # Add clusters if provided
    if opts.clusters:
        if opts.clusters != "None": # Hack for when called internally 
            mag_to_cluster = pd.read_csv(opts.clusters, sep="\t", index_col=0, header=None).iloc[:,0]
            orf_to_cluster = orf_to_mag.map(lambda id_orf: mag_to_cluster[id_orf])
            df_orf_classifications.insert(loc=2, column="id_cluster", value=orf_to_cluster)

    # Output
    if opts.output == "stdout":
        opts.output = sys.stdout
    df_orf_classifications = df_orf_classifications.dropna(how="any", axis=0)
    df_orf_classifications.to_csv(opts.output, sep="\t", header=bool(opts.header))

  


    

if __name__ == "__main__":
    main()
    
                

