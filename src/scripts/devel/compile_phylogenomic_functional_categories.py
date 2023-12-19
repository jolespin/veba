#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, pickle
from collections import defaultdict
# import numpy as np
import pandas as pd
from tqdm import tqdm 
import ensemble_networkx as enx

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.10.23"


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <annotation_results.tsv[.gz]> -l genome -o <output_table>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--annotation_results",  type=str, default = "stdin", help = "path/to/annotation.tsv from annotate.py [Default: stdin]")
    parser.add_argument("-X","--counts",  type=str, required = True, help = "path/to/X_orfs.tsv[.gz] from mapping.py at the ORF/gene/protein level.  Rows=Samples, Columns=Genes")
    parser.add_argument("-g","--genes",  type=str, help = "path/to/genes.ffn[.gz] fasta used for scaling-factors")
    parser.add_argument("-o","--output_directory", type=str, default="phylogenomic_functional_categories", help = "path/to/output_directory [Default: phylogenomic_functional_categories]")
    parser.add_argument("-l","--level", type=str, default="genome_cluster", help = "level {genome, genome_cluster} [Default: genome_cluster]")
    parser.add_argument("--minimum_count", type=float, default=1.0, help = "Minimum count to include gene [Default: 1 ]")
    parser.add_argument("--veba_database", type=str,  help = "VEBA Database [Default: $VEBA_DATABASE environment variable]")

    # parser.add_argument("-p", "--include_protein_identifiers", action="store_true", help = "Write protein identifiers")
    # parser.add_argument("--header", action="store_true", help = "Write header")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    assert opts.level in {"genome", "genome_cluster"}, "--level must be either {genome, genome_cluster}"

    if opts.level == "genome":
        level_field = ("Identifiers", "id_genome")
    if opts.level == "genome_cluster":
        level_field = ("Identifiers", "id_genome_cluster")

    if not opts.veba_database:
        opts.veba_database = os.environ["VEBA_DATABASE"]

    os.makedirs(opts.output_directory, exist_ok=True)
    # os.makedirs(os.path.join(opts.output_directory, opts.level), exist_ok=True)

    # Read annotations
    if opts.annotation_results == "stdin":
        opts.annotation_results = sys.stdin 
    df_annotations = pd.read_csv(opts.annotation_results, sep="\t", index_col=0, header=[0,1])
    protein_to_organism  = df_annotations[level_field]

    # KEGG Database
    delimiters = [",","_","-","+"]

    # Load MicrobeAnnotator KEGG dictionaries
    module_to_kos__unprocessed = defaultdict(set)
    for fp in glob.glob(os.path.join(opts.veba_database, , "*.pkl")):
        with open(fp, "rb") as f:
            d = pickle.load(f)
            
        for id_module, v1 in d.items():
            if isinstance(v1, list):
                try:
                    module_to_kos__unprocessed[id_module].update(v1)
                except TypeError:
                    for v2 in v1:
                        module_to_kos__unprocessed[id_module].update(v2)
            else:
                for k2, v2 in v1.items():
                    if isinstance(v2, list):
                        try:
                            module_to_kos__unprocessed[id_module].update(v2)
                        except TypeError:
                            for v3 in v2:
                                module_to_kos__unprocessed[id_module].update(v3)

    # Flatten the KEGG orthologs
    module_to_kos = dict()
    for id_module, kos_unprocessed in module_to_kos__unprocessed.items():
        kos_processed = set()
        for id_ko in kos:
            composite=False
            for sep in delimiters:
                if sep in id_ko:
                    id_ko = id_ko.replace(sep,";")
                    composite = True
            if composite:
                kos_composite = set(map(str.strip, filter(bool, id_ko.split(";"))))
                kos_processed.update(kos_composite)
            else:
                kos_processed.add(id_ko)
        module_to_kos[id_module] = kos_processed

    # Read counts
    X_counts = pd.read_csv(opts.counts, sep="\t", index_col=0)

    # Organisms
    organisms = df_annotations[level_field].unique()

    # Organizing KOs
    organism_to_kos = defaultdict(set)
    protein_to_kos = dict()
    kos_global = list()
    for id_protein, (id_organism, ko_ids) in tqdm(df_annotations.loc[:,[level_field, ("KOFAM", "ids")]].iterrows(), "Compiling KO identifiers", total=df_annotations.shape[0]):
        ko_ids = eval(ko_ids)
        if len(ko_ids):
            ko_ids = set(ko_ids)
            protein_to_kos[id_protein] = ko_ids
            organism_to_kos[id_organism].update(ko_ids)
            for id_ko in ko_ids:
                kos_global.append([id_protein, id_organism, id_ko])
    df_kos_global = pd.DataFrame(kos_global, columns=["id_protein", level_field[1], "id_kegg-ortholog"])
    del kos_global
    df_kos_global.to_csv(os.path.join(opts.output_directory, "kos.{}s.tsv".format(opts.level)), sep="\t", index=False)

    # Sample -> Organisms -> KOs
    sample_to_organism_to_kos = defaultdict(lambda: defaultdict(set))
    for id_sample, row in X_counts.iterrows():
        for id_protein, count in tqdm(row.items(), total=X_counts.shape[1]):
            if id_protein in protein_to_kos:
                if count >= opts.minimum_count:
                    id_organism = protein_to_organism[id_protein]
                    kos = protein_to_kos[id_protein]
                    sample_to_organism_to_kos[id_sample][id_organism].update(kos)

                    







if __name__ == "__main__":
    main()
