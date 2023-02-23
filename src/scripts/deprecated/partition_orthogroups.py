#!/usr/bin/env python
import sys, os, glob, argparse 
from collections import OrderedDict
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.04.01"


def restore_modified_identifier(x):
    fields = x.split("_")
    number_of_fields = len(fields)

    identifier_ok = False
    if len(fields) == 7: # Prodigal
        identifier_ok = True
        return x
    if len(fields) == 10: # Modified MetaEuk
        identifier_ok = True
        return "{}_{}_{}_{}_{}_{}_{}:{}({})".format(*fields)
    assert identifier_ok, "Could not infer identifier type for {}. Should either be either prodigal format: NODE_1460_length_4773_cov_9.241628 or modified MetaEuk: NODE_1460_length_4773_cov_9.241628_2314:2934(-)'".format(x)


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
    parser.add_argument("-i","--orthogroups", type=str, help = "path/to/orthogroups.tsv")
    parser.add_argument("-c","--clusters", type=str, help = "path/to/clusters.tsv")
    parser.add_argument("-a","--proteins", type=str, help = "path/to/proteins.tsv")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "Output table [Default: stdout]")
    parser.add_argument("-x", "--extension", type=str, default="faa", help="File extension for proteins [Default: 'faa")
    parser.add_argument("--clone_label", type=str, default="___clone", help="Cluster suffix [Default: '___clone")
    parser.add_argument("--proteins_to_contigs", type=str,  help="Assumes prodigal format for gene calls.  If not, then specify here.  No header. [id_protein]<tab>[id_contig]")
    # parser.add_argument("--metaeuk", type=str,  default="prodigal", help="gene_calling_algorithm {prodigal,metaeuk, infer}")

    parser.add_argument("--sep", type=str, default="_",  help = "Seperator for cluster and orthogroups. e.g. SLC0_OG0000000 with _ as the sep [Default: '_'")
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
    print(" * Reading orthogroups table:", opts.orthogroups, file=sys.stderr)
    opts.orthogroups = pd.read_csv(opts.orthogroups, sep="\t", index_col=0, header=None).iloc[:,0]
    opts.orthogroups.index = opts.orthogroups.index.map(str)

    print(" * Reading clusters table:", opts.clusters, file=sys.stderr)
    opts.clusters = pd.read_csv(opts.clusters, sep="\t", index_col=0, header=None).iloc[:,0].map(str)
    mag_to_cluster =  opts.clusters

    print(" * Reading proteins table:", opts.proteins, file=sys.stderr)
    opts.proteins = pd.read_csv(opts.proteins, sep="\t", index_col=0, header=None).iloc[:,0]
    assert set(opts.clusters.index) == set(opts.clusters.index), "--clusters and --proteins must have the same genomes in the first column"
    opts.proteins = opts.proteins.loc[opts.clusters.index]


    # Protein to MAG
    protein_to_mag = dict()
    for id_mag, path in tqdm(opts.proteins.items(), " * Parsing protein files", unit=" MAG"):
        with open(path, "r") as f_protein:
            for header, seq in SimpleFastaParser(f_protein):
                id_protein = header.split(" ")[0]
                # protein_to_contig[id_protein] 
                protein_to_mag[id_protein] = id_mag
    protein_to_mag = pd.Series(protein_to_mag)

    if opts.proteins_to_contigs:
        protein_to_contig = dict()
        with open(opts.proteins_to_contigs, "r") as f:
            for line in f.readlines():
                line = line.strip()
                if line:
                    id_protein, id_contig = line.split("\t")
                    protein_to_contig[id_protein] = id_contig
        protein_to_contig= pd.Series(protein_to_contig)
    else:
        protein_to_contig = pd.Series(
            data=protein_to_mag.index.map(lambda x: "_".join(x.split("_")[:-1])), 
            index=protein_to_mag.index,
        )

    orthogroup_to_cluster = dict() 
    protein_to_orthogroup = dict() 


    for id_cluster, path in tqdm(opts.orthogroups.items(), " * Parsing orthogroup files", unit=" orthogroup"):
        df_orthogroups = pd.read_csv(path, sep="\t", index_col=0)
        df_orthogroups.index = df_orthogroups.index.map(lambda x: "{}{}{}".format(id_cluster, opts.sep, x))
        df_orthogroups = df_orthogroups.applymap(lambda x: set() if pd.isnull(x) else set(x.split(", ")))
        for id_orthogroup, row in df_orthogroups.iterrows():
            orthogroup_to_cluster[id_orthogroup] = id_cluster
            for id_mag, proteins in row.items():
                if not id_mag.endswith(opts.clone_label):
                    for id_protein in proteins:
                        id_protein = restore_modified_identifier(id_protein) # Hack to handle modified MetaEuk labels
                        protein_to_orthogroup[id_protein] = id_orthogroup


    protein_to_orthogroup = pd.Series(protein_to_orthogroup)
    orthogroup_to_cluster = pd.Series(orthogroup_to_cluster)

    # Output table
    df_output = pd.concat([ 
        protein_to_contig.to_frame("id_contig"),
        protein_to_mag.to_frame("id_genome"),
    ], axis=1)

    df_output["id_cluster"] = protein_to_mag.map(lambda id_mag: mag_to_cluster[id_mag]).reindex(protein_to_mag.index)
    df_output["id_orthogroup"] = protein_to_orthogroup.reindex(protein_to_mag.index)
    df_output.index.name = "id_orf"

    if opts.output == "stdout":
        opts.output = sys.stdout
    df_output.to_csv(opts.output, sep="\t")



if __name__ == "__main__":
    main()
    
                

