#!/usr/bin/env python
import sys, os, glob, argparse, pickle, warnings, gzip
import pandas as pd
import numpy as np
from collections import OrderedDict
from tqdm import tqdm 

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.12.14"


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <metaeuk_identifier_mapping> -s <scaffolds_to_bins.tsv> -c <genome_to_clusters.tsv> -o <output_filepath>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--metaeuk_identifier_mapping", type=str, required=True, help = "path/to/identifier_mapping.metaeuk.tsv")
    parser.add_argument("-s","--scaffolds_to_bins", type=str, required=False, help = "path/to/scaffolds_to_bins.tsv")
    # parser.add_argument("-g","--genes_to_contigs", type=str, required=False, help = "path/to/genes_to_contigs.tsv cannot be used with --scaffolds_to_bins")
    parser.add_argument("-c","--clusters", type=str,  help = "path/to/clusters.tsv, Format: [id_genome]<tab>[id_cluster], No header [Optional]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/gene-source_lineage.tsv [Default: stdout]")
    parser.add_argument("-d", "--eukaryotic_database", type=str, default=None, required=True, help="path/to/eukaryotic_database directory (e.g. --arg 1 )")
    # parser.add_argument("--veba_database", type=str, default=None, help=f"VEBA database location.  [Default: $VEBA_DATABASE environment variable]")
    parser.add_argument("--header", type=int, default=1, help="Include header in output {0=No, 1=Yes) [Default: 1]")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--remove_genes_with_missing_values", action="store_true")
    parser.add_argument("--use_original_metaeuk_gene_identifiers", action="store_true")

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

        
    # SourceID -> Taxonomy
    fp = os.path.join(opts.eukaryotic_database,"source_taxonomy.tsv.gz")
    print("* Reading source taxonomy table {}".format(fp), file=sys.stderr)
    df_source_taxonomy = pd.read_csv(fp, sep="\t", index_col=0)
    df_source_taxonomy.index = df_source_taxonomy.index.map(str)
    df_source_taxonomy = pd.DataFrame(df_source_taxonomy.to_dict()) # Hack for duplicate entries that will be resolved in MicroEuk_v3.1
    
    if opts.debug:
        print(fp, file=sys.stderr)
        df_source_taxonomy.head().to_csv(sys.stderr, sep="\t")
        print("\n", file=sys.stderr)

    # VEBA -> SourceID
    fp = os.path.join(opts.eukaryotic_database,"target_to_source.dict.pkl.gz")
    print("* Reading target to source mapping {} (Note: This one takes a little longer to load...)".format(fp), file=sys.stderr)
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

    gene_to_bitscore = df_metaeuk["bitscore"].map(float)
    gene_to_scaffold = df_metaeuk["C_acc"].map(str)
    gene_to_genome = pd.Series([np.nan]*df_metaeuk.shape[0], index=df_metaeuk.index) 
    gene_to_target = df_metaeuk["T_acc"]
    gene_to_source = gene_to_target.map(lambda id_target: target_to_source.get(id_target,np.nan))

    if opts.scaffolds_to_bins:
        # Scaffolds -> Bins
        fp = opts.scaffolds_to_bins
        print("* Reading scaffolds to bins table {}".format(fp), file=sys.stderr)
        scaffold_to_bin = pd.read_csv(fp, sep="\t", index_col=0, header=None).iloc[:,0]
        if opts.debug:
            print(fp, file=sys.stderr)
            scaffold_to_bin.head().to_csv(sys.stderr, sep="\t", header=None)
            print("\n", file=sys.stderr)
        gene_to_genome = gene_to_scaffold.map(lambda id_scaffold: scaffold_to_bin[id_scaffold])

    if np.any(pd.isnull(gene_to_source)):
        warnings.warn("The following gene - target identifiers are not in the database file: {}".format(
            os.path.join(opts.eukaryotic_database,"target_to_source.dict.pkl.gz"), 
            ),
        )
        gene_to_target[gene_to_source[gene_to_source.isnull()].index].to_frame().to_csv(sys.stderr, sep="\t", header=None)
        gene_to_source = gene_to_source.dropna()

    # Lineage
    gene_to_lineage = OrderedDict()

    missing_lineage = list()
    for id_gene, id_source in tqdm(gene_to_source.items(), desc="Retrieving lineage", unit = " genes"):
        if id_source in df_source_taxonomy.index:
            lineage = df_source_taxonomy.loc[id_source, ["class", "order", "family", "genus", "species"]] #  class   order   family  genus   species
            lineage = lineage.fillna("")
            lineage = ";".join(map(lambda items: "".join(items), zip(["c__", "o__", "f__", "g__", "s__"], lineage)))
            gene_to_lineage[id_gene] = lineage
        else:
            missing_lineage.append(id_source)
    gene_to_lineage = pd.Series(gene_to_lineage)

    if len(missing_lineage):
        warnings.warn("The following source identifiers are not in the database file: {}\n{}`".format(
            os.path.join(opts.eukaryotic_database,"source_taxonomy.tsv.gz"), 
            "\n".join(set(missing_lineage)),
            ),
        )

    # Output
    df_gene_classifications = pd.DataFrame({
        "id_scaffold":gene_to_scaffold,
        "id_genome":gene_to_genome,
        "id_target":gene_to_target,
        "id_source":gene_to_source,
        "lineage":gene_to_lineage,
        "bitscore":gene_to_bitscore,
    })
    df_gene_classifications.index.name = "id_gene"


    # df_gene_classifications = pd.concat([
    #     gene_to_scaffold.to_frame("id_scaffold"),
    #     gene_to_genome.to_frame("id_genome"),
    #     gene_to_target.to_frame("id_target"),
    #     gene_to_source.to_frame("id_source"),
    #     gene_to_lineage.to_frame("lineage"),
    #     gene_to_bitscore.to_frame("bitscore"),
    # ],
    # axis=1)
    # df_gene_classifications.index.name = "id_gene"


    # Add clusters if provided
    if opts.clusters:
        if opts.clusters != "None": # Hack for when called internally 
            genome_to_cluster = pd.read_csv(opts.clusters, sep="\t", index_col=0, header=None).iloc[:,0]
            gene_to_cluster = gene_to_genome.map(lambda id_gene: genome_to_cluster[id_gene])
            df_gene_classifications.insert(loc=2, column="id_cluster", value=gene_to_cluster)

    # Output
    if opts.output == "stdout":
        opts.output = sys.stdout
    if opts.remove_genes_with_missing_values:
        df_gene_classifications = df_gene_classifications.dropna(how="any", axis=0)

    if not opts.use_original_metaeuk_gene_identifiers:
        metaeuk_to_gene = df_metaeuk["gene_id"].to_dict()
        df_gene_classifications.index = df_gene_classifications.index.map(lambda x: metaeuk_to_gene[x])

    df_gene_classifications.to_csv(opts.output, sep="\t", header=bool(opts.header))

  


    

if __name__ == "__main__":
    main()
    
                

