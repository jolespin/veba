#!/usr/bin/env python
import sys
import os
import glob
import argparse
from collections import defaultdict
import pandas as pd
from tqdm import tqdm
from pyexeggutor import (
    build_logger,
    write_pickle,
    open_file_writer,
    copy_file,
)

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2025.3.3"

def main(argv=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <path/to/veba_output/essentials/> -o <path/to/output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--veba_essentials", type=str,  help = "path/to/veba_output/essentials/")
    parser.add_argument("-o","--output_directory", type=str, help = "path/to/output_directory/")

    # Options
    opts = parser.parse_args(argv)
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Copy genomes
    quality_dataframes = list()
    taxonomy_series = list()
    protein_annotations = defaultdict(dict)

    logger = build_logger("Reformat VEBA for NAL Standards")
    os.makedirs(opts.output_directory, exist_ok=True)
    os.makedirs(os.path.join(opts.output_directory, "Genomes"), exist_ok=True)
    os.makedirs(os.path.join(opts.output_directory, "Metadata"), exist_ok=True)
    os.makedirs(os.path.join(opts.output_directory, "Traits"), exist_ok=True)
    
    for organism_type in ["eukaryotic", "prokaryotic", "viral"]:
        genome_directory = os.path.join(opts.veba_essentials, "genomes", organism_type)
        if os.path.exists(genome_directory):
            source_files = list(filter(os.path.isfile, glob.glob(os.path.join(genome_directory, "*"))))
            if source_files:
                f = open_file_writer(os.path.join(opts.output_directory, f"organisms.{organism_type}.list"), "w")
                destination_directory = os.path.join(opts.output_directory, "Genomes", organism_type.capitalize())
                os.makedirs(destination_directory, exist_ok=True)
                # =======
                # Genomes 
                # =======
                logger.info(f"Processing genomes: {organism_type}")
                for src in source_files:
                    logger.info(f"Copying: {src} to {destination_directory}")
                    copy_file(src, destination_directory)
                    if src.endswith(".fa.gz"):
                        id_genome = os.path.basename(src).rsplit(".", maxsplit=2)
                        print(id_genome, file=f)
                logger.info(f"Processing quality: {organism_type}")

                # =======
                # Quality
                # =======
                if organism_type == "eukaryotic":
                    df = pd.read_csv(f"{opts.veba_essentials}/quality/quality.{organism_type}.tsv.gz", sep="\t", index_col=0, header=[0,1])["specific"]
                    busco_version = None
                    if "Complete" in df.columns:
                        busco_version = "BUSCO_v5.4.x"
                        df = df[["Complete", "Multi copy", "dataset_name"]]
                        df.insert(2, "busco_version", busco_version)

                    elif "Complete percentage" in df.columns:
                        busco_version = "BUSCO_v5.6.x"
                        df = df[["Complete percentage", "Multi copy percentage", "dataset_name"]]
                        df.insert(2, "busco_version", busco_version)
                    df.columns = range(df.shape[1])
                    df.to_csv(os.path.join(opts.output_directory,"Metadata", f"quality.{organism_type}.tsv.gz"), sep="\t", header=None)
                    quality_dataframes.append(df)
                    
                if organism_type == "prokaryotic":
                    df = pd.read_csv(f"{opts.veba_essentials}/quality/quality.{organism_type}.tsv.gz", sep="\t", index_col=0)
                    if "Origin" in df.columns:
                        df = df[["Completeness", "Contamination", "Origin"]]
                        df.insert(2, "method", "Binette_v1.0.5")
                    elif "Completeness_Model_Used" in df.columns:
                        df = df[["Completeness", "Contamination", "Completeness_Model_Used"]]
                        df.insert(2, "method", "CheckM2_v1.0.1")
                    df.columns = range(df.shape[1])
                    df.to_csv(os.path.join(opts.output_directory,"Metadata", f"quality.{organism_type}.tsv.gz"), sep="\t", header=None)
                    quality_dataframes.append(df)
                    
                if organism_type == "viral":
                    df = pd.read_csv(f"{opts.veba_essentials}/quality/quality.{organism_type}.tsv.gz", sep="\t", index_col=0)
                    info = df[["checkv_quality", "miuvig_quality","completeness_method"]].apply(lambda x: ", ".join(map(lambda y: "=".join(y), x.items())), axis=1)
                    df = df[["completeness", "contamination"]]
                    df["checkv_version"] = "CheckV_v1.0.1"
                    df["info"] = info
                    df.columns = range(df.shape[1])
                    df.to_csv(os.path.join(opts.output_directory,"Metadata", f"quality.{organism_type}.tsv.gz"), sep="\t", header=None)
                    quality_dataframes.append(df)
                    
            # =======
            # Taxonomy
            # =======
            taxonomy_filepath = os.path.join(opts.veba_essentials, "taxonomy_classification", organism_type, "taxonomy.tsv.gz")
            if os.path.exists(taxonomy_filepath):
                logger.info(f"Processing taxonomy: {organism_type}")

                if organism_type == "eukaryotic":
                    header = "# VEBA_v2.x.x"
                if organism_type == "prokaryotic":
                    header = "# GTDB_r220"
                if organism_type == "viral":
                    header = "# geNomad_v1.5.x"

                data = pd.read_csv(taxonomy_filepath, sep="\t", index_col=0).iloc[:,0]
                with open_file_writer(os.path.join(opts.output_directory,"Metadata", f"taxonomy.{organism_type}.tsv.gz"), "w") as f:
                    print(f"{header}", file=f)
                    data.to_csv(f, sep="\t", header=None, mode="a")
                taxonomy_series.append(data)
    df_quality = pd.concat(quality_dataframes)
    df_quality.to_csv(os.path.join(opts.output_directory,"Metadata", f"quality.tsv.gz"), sep="\t", header=None)
            
    genome_to_taxonomy = pd.concat(taxonomy_series)
    genome_to_taxonomy.to_frame().to_csv(os.path.join(opts.output_directory,"Metadata", f"taxonomy.tsv.gz"), sep="\t", header=None)

    # ===========
    # Annotations
    # ===========
    protein_annotations = defaultdict(dict)
    df = pd.read_csv(os.path.join(opts.veba_essentials, "annotation", "annotations.proteins.tsv.gz"), sep="\t", index_col=0, header=[0,1])
    mask = df[[("KOfam", "number_of_hits"), ("Enzymes","number_of_hits"),("Pfam","number_of_hits"), ('NCBIfam-AMR', "number_of_hits")]].sum(axis=1) > 0
    for id_protein, row in tqdm(df.loc[mask,["KOfam", "Enzymes", "Pfam", 'NCBIfam-AMR']].iterrows(), total=mask.sum(), desc="Parsing annotations", unit=" proteins"):
        if row["KOfam"]["number_of_hits"] > 0:
            protein_annotations[id_protein]["KOfam"] = eval(row["KOfam"]["ids"])
        else:
            protein_annotations[id_protein]["KOfam"] = set()

        if row["Enzymes"]["number_of_hits"] > 0:
            protein_annotations[id_protein]["EC"] = eval(row["Enzymes"]["ids"])
        else:
            protein_annotations[id_protein]["EC"] = set()

        if row["Pfam"]["number_of_hits"] > 0:
            protein_annotations[id_protein]["Pfam"] = eval(row["Pfam"]["ids"])
        else:
            protein_annotations[id_protein]["Pfam"] = set()

        if row["NCBIfam-AMR"]["number_of_hits"] > 0:
            protein_annotations[id_protein]["NCBIfam-AMR"] = eval(row["NCBIfam-AMR"]["ids"])
        else:
            protein_annotations[id_protein]["NCBIfam-AMR"] = set()

    df_annot = pd.DataFrame(protein_annotations).T
    df_annot.index.name = "id_protein"
    df_annot.to_csv(os.path.join(opts.output_directory,"Metadata", "protein_annotations.tsv.gz"), sep="\t")
    write_pickle(protein_annotations, os.path.join(opts.output_directory,"Metadata", "protein_annotations.pkl.gz"))


    df_identifier_mapping = pd.read_csv(os.path.join(opts.veba_essentials,"clustering", "identifier_mapping.proteins.tsv.gz"), sep="\t", index_col=0)
    protein_to_genome = df_identifier_mapping["id_genome"]
    del df_identifier_mapping
    df_identifier_mapping = protein_to_genome.to_frame()
    df_identifier_mapping.to_csv(os.path.join(opts.output_directory,"Metadata", "identifier_mapping.lite.tsv.gz"), sep="\t", header=None)
    df_identifier_mapping.insert(0, "id_contig", df_identifier_mapping.index.map(lambda x: x.rsplit("_", maxsplit=1)[0]))
    df_identifier_mapping.to_csv(os.path.join(opts.output_directory,"Metadata", "identifier_mapping.tsv.gz"), sep="\t", header=None)

    genomic_traits = defaultdict(dict)
    for id_protein, annotations in tqdm(protein_annotations.items(), desc="Building traits matrix", unit=" proteins"):
        kofam = annotations["KOfam"]
        if kofam:
            id_genome = protein_to_genome[id_protein]
            for id_ko in kofam:
                genomic_traits[id_genome][id_ko] = 1

    df_genomic_traits = pd.DataFrame(genomic_traits).T.fillna(0).astype(int)
    df_genomic_traits.to_csv(os.path.join(opts.output_directory,"Traits", "genomic_traits.kofam.bool-int.tsv.gz"), sep="\t")


if __name__ == "__main__":
    main()
    
                

