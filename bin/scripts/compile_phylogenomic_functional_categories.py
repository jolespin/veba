#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, pickle, time
from collections import defaultdict
# import numpy as np
import pandas as pd
from tqdm import tqdm 
import numpy as np
from ensemble_networkx import CategoricalEngineeredFeature
from soothsayer_utils import read_fasta, format_duration
from genopype import *

script_directory  =  os.path.dirname(os.path.abspath( __file__ ))

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.2.5"


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
    parser.add_argument("-f","--fasta",  type=str, help = "path/to/sequences.faa[.gz] fasta used for scaling-factors.  If --level is genome, this should be proteins from ORFs.  If --level is genome_cluster, this should be the representative gene sequence using the SSPC label as the identifier.")
    parser.add_argument("-o","--output_directory", type=str, default="veba_output/phylogenomic_functional_categories", help = "path/to/output_directory [Default: veba_output/phylogenomic_functional_categories]")
    parser.add_argument("-l","--level", type=str, default="genome_cluster", help = "level {genome, genome_cluster} [Default: genome_cluster]")
    parser.add_argument("-n","--name", type=str,  help = "Name of CategoricalEngineerFeature object (e.g., TestDataset-Metagenomics or CohortX-Metatranscriptomics)")
    parser.add_argument("--minimum_count", type=float, default=1.0, help = "Minimum count to include gene [Default: 1 ]")
    parser.add_argument("--multiplier",  type=int, choices={1,3}, default=3, help = "Multiplier for scaling factors.  If genes in nucleotide-space are used then this should 1 and genes in protein-space are used this should be 3. [Default: 3]")
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
    os.makedirs(os.path.join(opts.output_directory, "intermediate"), exist_ok=True)

    with open(os.path.join(opts.output_directory, "phylogenomic_functional_categories.log"), "w") as f_log:
        # os.makedirs(os.path.join(opts.output_directory, opts.level), exist_ok=True)

        # Read annotations
        if opts.annotation_results == "stdin":
            opts.annotation_results = sys.stdin 
        df_annotations = pd.read_csv(opts.annotation_results, sep="\t", index_col=0, header=[0,1])

        # KEGG Database
        delimiters = [",","_","-","+"]

        # Load MicrobeAnnotator KEGG dictionaries
        module_to_kos__unprocessed = defaultdict(set)
        for fp in glob.glob(os.path.join(opts.veba_database, "Annotate", "MicrobeAnnotator-KEGG", "*.pkl")):
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
            for id_ko in kos_unprocessed:
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

        # KOS to Module
        ko_to_modules = defaultdict(set)
        for id_module, kos in module_to_kos.items():
            for id_ko in kos:
                ko_to_modules[id_ko].add(id_module)

        # Read counts
        X_counts = pd.read_csv(opts.counts, sep="\t", index_col=0)

        # Organisms
        organisms = df_annotations[level_field].unique()

        # Organizing KOs
        organism_to_kos = defaultdict(set)
        protein_to_kos = dict()
        protein_to_modules = defaultdict(set)
        kos_global = list()
        for id_protein, (id_organism, ko_ids) in tqdm(df_annotations.loc[:,[level_field, ("KOFAM", "ids")]].iterrows(), "Compiling KO identifiers", total=df_annotations.shape[0]):
            ko_ids = eval(ko_ids)
            if len(ko_ids):
                ko_ids = set(ko_ids)
                protein_to_kos[id_protein] = ko_ids
                organism_to_kos[id_organism].update(ko_ids)
                for id_ko in ko_ids:
                    kos_global.append([id_protein, id_organism, id_ko])
                    if id_ko in ko_to_modules:
                        protein_to_modules[id_protein].update(ko_to_modules[id_ko])
                    else:
                        print("[{}][{}]".format(id_organism, id_protein),  "No KEGG module information for {} in {}".format(id_ko, os.path.join(opts.veba_database, "Annotate", "MicrobeAnnotator-KEGG")), file=f_log)
        

        df_kos_global = pd.DataFrame(kos_global, columns=[{"genome":"id_protein","genome_cluster":"id_protein_cluster"}[opts.level], level_field[1], "id_kegg-ortholog"])
        del kos_global
        df_kos_global.to_csv(os.path.join(opts.output_directory, "intermediate", "kos.{}s.tsv".format(opts.level)), sep="\t", index=False)

        A = set(protein_to_modules.keys())
        B = set(X_counts.columns)
        # assert A <= B, "All the genes/proteins in --annotation_results must be present in --counts.  There are {} unique to --annotation_results and {} unique to --counts with {} shared.".format(len(A - B), len(B-A), len(A & B))
        print("There are {} unique to --annotation_results and {} unique to --counts with {} shared.".format(len(A - B), len(B-A), len(A & B)), file=f_log)


        # Index of proteins with KEGG module info
        protein_to_modules = pd.Series(protein_to_modules)
        protein_to_organism  = df_annotations[level_field]

        # Sample -> Organisms -> KOs
        sample_to_organism_to_kos = defaultdict(lambda: defaultdict(set))
        for id_sample, row in X_counts.iterrows():
            for id_protein, count in tqdm(row.items(), total=X_counts.shape[1]):
                if id_protein in protein_to_kos:
                    if count >= opts.minimum_count:
                        id_organism = protein_to_organism[id_protein]
                        kos = protein_to_kos[id_protein]
                        sample_to_organism_to_kos[id_sample][id_organism].update(kos)

        # Sample-specific KO tables
        os.makedirs(os.path.join(opts.output_directory, "intermediate", "sample-specific_ko_tables"), exist_ok=True)
        for id_sample, organism_to_kos in sample_to_organism_to_kos.items():
            with open(os.path.join(opts.output_directory, "intermediate", "sample-specific_ko_tables", "{}.tsv".format(id_sample)), "w") as f_kos:
                for id_organism, kos in organism_to_kos.items():
                    for id_ko in kos:
                        print(id_organism, id_ko, sep="\t", file=f_kos)

        scaling_factors = None
        if opts.fasta:
            gene_to_length = read_fasta(opts.fasta, description=False, verbose=True)
            gene_to_length = gene_to_length.map(len)
            assert set(protein_to_modules.keys()) <= set(gene_to_length.index), "All the genes/proteins in --annotation_results must be present in --genes"
            scaling_factors = gene_to_length.loc[protein_to_modules.index]
            scaling_factors = scaling_factors * opts.multiplier

        # Feature engineering
        protein_to_organism = protein_to_organism.loc[protein_to_modules.index]

        cef = CategoricalEngineeredFeature(
            initial_feature_type="protein_cluster" if opts.level == "genome_cluster" else "protein",
            engineered_feature_type="pgfc",
            observation_type="id_sample",
            unit_type="normalized_counts" if opts.fasta else "counts",
            name=opts.name,
        )

        cef.add_category(
            name_category="Taxonomy",
            mapping=protein_to_organism,
        )
        cef.add_category(
            name_category="Functional",
            mapping=protein_to_modules,
        )
        cef.compile(scaling_factors=scaling_factors)
        cef.to_file(os.path.join(opts.output_directory, "categorical_engineered_feature.pgfc.{}s.pkl".format(opts.level)))    
        cef.synopsis_.to_csv(os.path.join(opts.output_directory, "categorical_engineered_feature.pgfc.{}s.synopsis.tsv.gz".format(opts.level)), sep="\t")

        X_transformed = cef.fit_transform(
            X=X_counts, 
            aggregate_fn=np.sum,
        )
        X_transformed = X_transformed.loc[X_counts.index]
        X_transformed.to_csv(os.path.join(opts.output_directory, "categorical_engineered_feature.pgfc.{}s.transformed_counts.tsv.gz".format(opts.level)), sep="\t")



        # Module completion ratios
        with open(os.path.join(opts.output_directory, "intermediate", "commands.sh"), "w") as f_cmds:
            # Directories
            t0 = time.time()
            directories = dict()
            directories["output"] = create_directory(os.path.join(opts.output_directory, "intermediate", "output"))
            directories["log"] = create_directory(os.path.join(opts.output_directory, "intermediate", "log"))
            # directories["tmp"] = create_directory(os.path.join(opts.output_directory, "intermediate", "tmp"))
            directories["checkpoints"] = create_directory(os.path.join(opts.output_directory, "intermediate", "checkpoints"))

            for fp in glob.glob(os.path.join(opts.output_directory, "intermediate", "sample-specific_ko_tables", "*.tsv")):
                id_sample = fp.split("/")[-1][:-4]
                name = "mcr__{}".format(id_sample)
                description = "[Module Completion Ratios] [id_sample={}]".format(id_sample)
                cmd = Command([
                    os.path.join(script_directory, "module_completion_ratios.py"),
                    "-i {}".format(fp),
                    "-o {}".format(os.path.join(directories["output"], "{}.tsv".format(id_sample))),
                    "-d {}".format(os.path.join(opts.veba_database, "Annotate", "MicrobeAnnotator-KEGG")),
                    ], 
                    name=name, 
                    f_cmds=f_cmds,
                    )

                # Run command
                cmd.run(
                    checkpoint_message_notexists="[Running ({})] | {}".format(format_duration(t0), description),
                    checkpoint_message_exists="[Loading Checkpoint ({})] | {}".format(format_duration(t0), description),
                    write_stdout=os.path.join(directories["log"], "{}.o".format(name)),
                    write_stderr=os.path.join(directories["log"], "{}.e".format(name)),
                    write_returncode=os.path.join(directories["log"], "{}.returncode".format(name)),
                    checkpoint=os.path.join(directories["checkpoints"], name),
                    )
                if hasattr(cmd, "returncode_"):
                    if cmd.returncode_ != 0:
                        print("[Error] | {}".format(description), file=sys.stdout)
                        print("Check the following files:\ncat {}".format(os.path.join(directories["log"], "{}.*".format(name))), file=sys.stdout)
                        sys.exit(cmd.returncode_)

        # Concatenate MCRs
        sample_to_mcrs = dict()
        for fp in glob.glob(os.path.join(opts.output_directory, "intermediate", "output", "*.tsv")):
            id_sample = fp.split("/")[-1][:-4]
            df = pd.read_csv(fp, sep="\t", index_col=1)
            df = df.drop(["pathway_group", "module_name"], axis=1)
            mcrs = df.stack()
            mcrs.index = mcrs.index.map(lambda x: x[::-1])
            sample_to_mcrs[id_sample] = mcrs
        df_mcr = pd.DataFrame(sample_to_mcrs).T
        df_mcr.columns.names = ["Taxonomy", "Functional"]
        df_mcr.index.name = "id_sample"
        df_mcr = df_mcr.reindex(index=X_transformed.index, columns=X_transformed.columns).fillna(0.0)
        df_mcr.to_csv(os.path.join(opts.output_directory, "categorical_engineered_feature.pgfc.{}s.module_completion_ratios.tsv.gz".format(opts.level)), sep="\t")




if __name__ == "__main__":
    main()
