#!/usr/bin/env python
import sys, os, glob, argparse 
from collections import defaultdict
import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.10.11"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <classify_directory> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--classify_directory", type=str, default="veba_output/classify", help = "path/to/classify_directory.  Assumes directory structure is [classification_directory]/[domain]/output/taxonomy.tsv/[taxonomy.clusters.tsv] [Default: veba_output/classify]")
    parser.add_argument("-o","--output_directory", type=str, help = "path/to/output_directory [Default: veba_output/classify]", default="veba_output/classify")
    parser.add_argument("-d","--no_domain", action="store_true", help = "Don't domain in output")
    parser.add_argument("--no_header", action="store_true", help = "Don't write header")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    os.makedirs(opts.output_directory, exist_ok=True)
    
    # Get cluster taxonomy
    df_taxonomy_clusters = None
    cluster_dataframes = list()
    cluster_taxonomy = glob.glob(os.path.join(opts.classify_directory, "*", "output", "taxonomy.clusters.tsv"))
    if cluster_taxonomy:
        print("* Compiling taxonomy from the following files:", *cluster_taxonomy, sep="\n    ", file=sys.stdout)
        for fp in cluster_taxonomy:
            id_domain = fp.split("/")[-3]
            df = pd.read_csv(fp, sep="\t", index_col=0)
            if not opts.no_domain:
                df.insert(loc=0, column="domain", value=id_domain)
            cluster_dataframes.append(df)
        df_taxonomy_clusters = pd.concat(cluster_dataframes, axis=0)
        df_taxonomy_clusters.index.name = "id_genome_cluster"
        print("\tWriting genome cluster taxonomy output for {} clusters: {}".format(df_taxonomy_clusters.shape[0], os.path.join(opts.output_directory, "taxonomy_classifications.clusters.tsv")), file=sys.stdout)
        df_taxonomy_clusters.to_csv(os.path.join(opts.output_directory, "taxonomy_classifications.clusters.tsv"), sep="\t",  header=not bool(opts.no_header))

    else:
        print("Could not find any cluster taxonomy tables in the following directory: {}".format(opts.classify_directory), file=sys.stdout)

      # Get genome taxonomy
    df_taxonomy_genomes = None
    genome_to_data = defaultdict(dict)
    genome_taxonomy = glob.glob(os.path.join(opts.classify_directory, "*", "output", "taxonomy.tsv"))
    if genome_taxonomy:
        print("* Compiling taxonomy from the following files:", *genome_taxonomy, sep="\n    ", file=sys.stdout)
        for fp in genome_taxonomy:
            id_domain = fp.split("/")[-3]
            df = pd.read_csv(fp, sep="\t", index_col=0)
            if id_domain.lower() in {"viral", "virus"}:
                for id_genome, taxonomy in df["lineage"].items():
                    genome_to_data[id_genome] = {"domain":id_domain, "taxonomy_classification":taxonomy}
            if id_domain.lower() in {"prokaryotic", "prokaryotes", "prokarya", "bacteria", "archaea", "bacterial","archael", "prok", "proks"}:
                for id_genome, taxonomy in df["classification"].items():
                    genome_to_data[id_genome] = {"domain":id_domain, "taxonomy_classification":taxonomy}
            if id_domain.lower() in {"eukaryotic", "eukaryotes", "eukarya", "microeukaryotes", "microeukarya", "microeukaryotic", "protists", "protista", "euk", "euks"}:
                for id_genome, taxonomy in df["consensus_classification"].items():
                    genome_to_data[id_genome] = {"domain":id_domain, "taxonomy_classification":taxonomy}
        df_taxonomy_genomes = pd.DataFrame(genome_to_data).T
        df_taxonomy_genomes = df_taxonomy_genomes.loc[:,["domain", "taxonomy_classification"]]
        if opts.no_domain:
            df_taxonomy_genomes = df_taxonomy_genomes.drop(["domain"], axis=1)


        df_taxonomy_genomes.index.name = "id_genome"
        # if df_taxonomy_clusters is not None:
        #     mag_to_slc = dict()
        #     for id_cluster, genomes in df_taxonomy_clusters["genomes"].items():
        #         genomes = eval(genomes)
        #         for id_genome in genomes:
        #             mag_to_slc[id_genome] = id_cluster
        #     mag_to_slc = pd.Series(mag_to_slc)
        #     df_taxonomy_genomes.insert(loc=1, column="id_genome_cluster", value=mag_to_slc)
        print("\tWriting genome taxonomy output for {} genomes: {}".format(df_taxonomy_genomes.shape[0], os.path.join(opts.output_directory, "taxonomy_classifications.tsv")), file=sys.stdout)
        df_taxonomy_genomes.to_csv(os.path.join(opts.output_directory, "taxonomy_classifications.tsv"), sep="\t",  header=not bool(opts.no_header))

    else:
        print("Could not find any taxonomy tables in the following directory: {}".format(opts.classify_directory), file=sys.stdout)


if __name__ == "__main__":
    main()
    
                

