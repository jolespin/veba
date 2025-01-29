#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, shutil, gzip
from collections import defaultdict
import pandas as pd
from tqdm import tqdm 
from pyexeggutor import (
    copy_file,
    gzip_file,
    build_logger,
)

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2025.1.24"
    
def merge_quality(organism_type, organism_specific_binning_directory):
    """
    Merges quality assessment files based on the organism type and directory structure.

    Parameters:
        organism_type (str): Type of organism ("prokaryotic", "viral", or "eukaryotic").
        organism_specific_binning_directory (str): Path to the organism-specific binning directory.

    Returns:
        pd.DataFrame: Merged DataFrame containing all quality assessment results.
    """
    dataframes = []

    if organism_type == "prokaryotic":
        file_pattern = os.path.join(organism_specific_binning_directory, "*", "output", "checkm2_results.filtered.tsv")
        header = 0

    elif organism_type == "viral":
        file_pattern = os.path.join(organism_specific_binning_directory, "*", "output", "checkv_results.filtered.tsv")
        header = 0

    elif organism_type == "eukaryotic":
        file_pattern = os.path.join(organism_specific_binning_directory, "*", "output", "busco_results.filtered.tsv")
        header = [0, 1]

    else:
        raise ValueError(f"Unsupported organism type: {organism_type}")

    for filepath in tqdm(glob.glob(file_pattern), desc=f"Processing {organism_type} files"):
        id_sample = filepath.split("/")[-3]
        df = pd.read_csv(filepath, sep="\t", index_col=0, header=header)
        if organism_type == "eukaryotic":
            df.insert(0, ("identifiers", "id_sample"), id_sample)
        else:
            df.insert(0, "id_sample", id_sample)
        dataframes.append(df)

    # Concatenate all DataFrames
    if dataframes:
        return pd.concat(dataframes, axis=0)
    else:  
        return pd.DataFrame()
    
def merge_identifier_mapping(veba_directory):
    dataframes = list()
    for filepath in tqdm(glob.glob(os.path.join(veba_directory, "binning","*", "*", "output", "genomes", "identifier_mapping.tsv"))):
        df = pd.read_csv(filepath, sep="\t", index_col=0, header=None)
        dataframes.append(df)
    # Concatenate all DataFrames
    if dataframes:
        return pd.concat(dataframes, axis=0)
    else:  
        return pd.DataFrame()
    
def merge_identified_mapping_metaeuk(eukaryotic_binning_directory):
    dataframes = list()
    for filepath in tqdm(glob.glob(os.path.join(eukaryotic_binning_directory, "*", "output", "genomes", "identifier_mapping.metaeuk.tsv"))):
        df = pd.read_csv(filepath, sep="\t", index_col=0)
        dataframes.append(df)
    # Concatenate all DataFrames
    if dataframes:
        return pd.concat(dataframes, axis=0)
    else:  
        return pd.DataFrame()
    

def merge_plasmid_summary(viral_binning_directory):
    dataframes = list()
    for filepath in tqdm(glob.glob(os.path.join(viral_binning_directory, "*", "output", "plasmid_summary.tsv"))):
        df = pd.read_csv(filepath, sep="\t", index_col=0)
        id_sample = filepath.split("/")[-3]
        df.insert(0, "id_sample", id_sample)
        dataframes.append(df)
    # Concatenate all DataFrames
    if dataframes:
        return pd.concat(dataframes, axis=0)
    else:  
        return pd.DataFrame()

def determine_preprocess_type(preprocess_directory):
    """
    Determines the preprocessing type based on file patterns.

    Parameters:
        patterns (dict): A dictionary containing file patterns.

    Returns:
        tuple: A tuple indicating the type of data and preprocessing (e.g., ("paired", "cleaned")).

    Raises:
        ValueError: If file pattern constraints are violated.
    """
    patterns = dict(
    cleaned_r1_filepaths = glob.glob(os.path.join(preprocess_directory, "*", "output", "cleaned_1.fastq.gz")),
    cleaned_r2_filepaths = glob.glob(os.path.join(preprocess_directory, "*", "output", "cleaned_2.fastq.gz")),
    trimmed_r1_filepaths = glob.glob(os.path.join(preprocess_directory, "*", "output", "trimmed_1.fastq.gz")),
    trimmed_r2_filepaths = glob.glob(os.path.join(preprocess_directory, "*", "output", "trimmed_2.fastq.gz")),
    cleaned_long_filepaths = glob.glob(os.path.join(preprocess_directory, "*", "output", "cleaned.fastq.gz")),
    trimmed_long_filepaths = glob.glob(os.path.join(preprocess_directory, "*", "output", "trimmed.fastq.gz")),
    )
    
    cleaned_r1 = patterns.get("cleaned_r1_filepaths", [])
    cleaned_r2 = patterns.get("cleaned_r2_filepaths", [])
    trimmed_r1 = patterns.get("trimmed_r1_filepaths", [])
    trimmed_r2 = patterns.get("trimmed_r2_filepaths", [])
    cleaned_long = patterns.get("cleaned_long_filepaths", [])
    trimmed_long = patterns.get("trimmed_long_filepaths", [])

    if cleaned_r1:
        if not cleaned_r2 or len(cleaned_r1) != len(cleaned_r2):
            raise ValueError("cleaned_r1_filepaths and cleaned_r2_filepaths must both be non-empty and have the same length.")
        if trimmed_r1 or trimmed_r2 or cleaned_long or trimmed_long:
            raise ValueError("Other file patterns must be empty when using cleaned_r1_filepaths.")
        return ("paired", "cleaned")

    if trimmed_r1:
        if not trimmed_r2 or len(trimmed_r1) != len(trimmed_r2):
            raise ValueError("trimmed_r1_filepaths and trimmed_r2_filepaths must both be non-empty and have the same length.")
        if cleaned_r1 or cleaned_r2 or cleaned_long or trimmed_long:
            raise ValueError("Other file patterns must be empty when using trimmed_r1_filepaths.")
        return ("paired", "trimmed")

    if cleaned_long:
        if cleaned_r1 or cleaned_r2 or trimmed_r1 or trimmed_r2 or trimmed_long:
            raise ValueError("Other file patterns must be empty when using cleaned_long_filepaths.")
        return ("long", "cleaned")

    if trimmed_long:
        if cleaned_r1 or cleaned_r2 or trimmed_r1 or trimmed_r2 or cleaned_long:
            raise ValueError("Other file patterns must be empty when using trimmed_long_filepaths.")
        return ("long", "trimmed")

    raise ValueError("No valid preprocessing type determined from the provided patterns.")


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <veba_directory> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i","--veba_directory", type=str, help = "path/to/veba_output/", required=True)
    parser_io.add_argument("-o","--output_directory", type=str,  default="veba_output/essentials/", help = "Output directory [Default: veba_output/essentials/]")
    
    parser_preprocess = parser.add_argument_group('Preprocess arguments')
    parser_preprocess.add_argument("-q", "--include_fastq", action="store_true", help="Include trimmed or cleaned fastq files")

    parser_assembly = parser.add_argument_group('Assembly arguments')
    parser_assembly.add_argument("-c", "--include_scaffolds", action="store_true", help="Include scaffolds assembly fasta")
    parser_assembly.add_argument("-b", "--include_bam", action="store_true", help="Inlcude mapped.sorted.bam file")
    
    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    # logger
    logger = build_logger("Gather VEBA essentials")
    
    # Output
    os.makedirs(opts.output_directory, exist_ok=True)

    # Gather essential files
    # ======================
    # Binning
    # -------
    parent_directory = os.path.join(opts.veba_directory, "binning")
    if os.path.exists(parent_directory):
        logger.info(f"Gathering binning results: {parent_directory}")
        # Identifier mapping
        df_identifier_mapping = merge_identifier_mapping(opts.veba_directory)
        if not df_identifier_mapping.empty:
            destination_filepath = os.path.join(opts.output_directory, "genomes", "identifier_mapping.tsv.gz")
            logger.info(f"Merging identifier mapping table: {destination_filepath}")
            os.makedirs(os.path.join(opts.output_directory, "genomes"), exist_ok=True)
            df_identifier_mapping.to_csv(destination_filepath, sep="\t", header=None)
        else:
            logger.warning("No identifier mapping files found")
        
        # Genomes
        for source_directory in glob.glob(os.path.join(opts.veba_directory, "binning", "*")):
            source_directory = os.path.normpath(source_directory)
            if os.path.isdir(source_directory):
                # if source_directory.endswith("/"):
                #     source_directory = source_directory[:-1]
                organism_specific_binning_directory = source_directory
                organism_type = organism_specific_binning_directory.split("/")[-1]
                logger.info(f"Gathering genome files [{organism_type}]: {source_directory}")
                destination_directory = os.path.join(opts.output_directory, "genomes", organism_type)
                os.makedirs(destination_directory, exist_ok=True)
                for source_filepath in tqdm(glob.glob(os.path.join(organism_specific_binning_directory, "*", "output", "genomes", "*")), organism_specific_binning_directory):
                    if not source_filepath.endswith((".tsv.gz", ".tsv")):
                        if source_filepath.endswith((".gz", ".fai")):
                            copy_file(source_filepath, destination_directory, gzip=False, logger=logger)
                        else:
                            copy_file(source_filepath, destination_directory, gzip=True, logger=logger)
                    

                # Quality assessment                 
                df_merged = merge_quality(organism_type, organism_specific_binning_directory)
                if not df_merged.empty:
                    destination_filepath = os.path.join(opts.output_directory, "quality", f"quality.{organism_type}.tsv.gz")
                    logger.info(f"Merging quality assessment for {organism_type}: {destination_filepath}")
                    os.makedirs(os.path.join(opts.output_directory, "quality"), exist_ok=True)
                    df_merged.to_csv(destination_filepath, sep="\t")
                else:
                    logger.warning(f"No quality assessment files found for {organism_type}")
                    
                # Plasmids
                if organism_type == "viral":
                    df_plasmids = merge_plasmid_summary(organism_specific_binning_directory)
                    if not df_plasmids.empty:
                        destination_filepath = os.path.join(opts.output_directory, "genomes", "plasmid_summary.tsv.gz")
                        logger.info(f"Merging plasmid summary: {destination_filepath}")
                        df_plasmids.to_csv(destination_filepath, sep="\t")
                    else:
                        logger.warning("No plasmid summary files found")
                        
                # MetaEuk
                if organism_type == "eukaryotic":
                    df_identifier_mapping_metaeuk = merge_identified_mapping_metaeuk(organism_specific_binning_directory)
                    if not df_identifier_mapping_metaeuk.empty:
                        destination_filepath = os.path.join(destination_directory, "identifier_mapping.metaeuk.tsv.gz")
                        logger.info(f"Merging MetaEuk identifier mapping: {destination_filepath}")
                        df_identifier_mapping_metaeuk.to_csv(destination_filepath, sep="\t")
                    else:
                        logger.warning("Identifier mapping files found for MetaEuk")
                        
    # Annotations
    # -----------
    source_directory = os.path.join(opts.veba_directory, "annotation")
    if os.path.exists(source_directory):
        logger.info(f"Gathering annotation results: {source_directory}")
        destination_directory = os.path.join(opts.output_directory, "annotation")
        os.makedirs(destination_directory, exist_ok=True)
        for source_filepath in tqdm(glob.glob(os.path.join(opts.veba_directory, "annotation","output", "*")), source_directory):
            if source_filepath.endswith(".gz"):
                copy_file(source_filepath, destination_directory, gzip=False, logger=logger)
            else:
                copy_file(source_filepath, destination_directory, gzip=True, logger=logger)
    # Clustering
    # ----------
    source_directory = os.path.join(opts.veba_directory, "cluster")
    if os.path.exists(source_directory):
        logger.info(f"Gathering clustering results: {source_directory}")
        destination_directory = os.path.join(opts.output_directory, "clustering")
        os.makedirs(destination_directory, exist_ok=True)
        #! Need to gzip files
        for source_filepath in tqdm(glob.glob(os.path.join(opts.veba_directory, "cluster","output", "global", "*")), source_directory):
            if os.path.isdir(source_filepath):
                if source_filepath.endswith("pangenome_core_sequences"):
                    destination_subdirectory = os.path.join(opts.output_directory, "clustering", "pangenome_core_sequences")
                    os.makedirs(destination_subdirectory, exist_ok=True)
                    for source_subfilepath in tqdm(glob.glob(os.path.join(source_filepath, "*")), source_filepath):
                        copy_file(source_subfilepath, destination_subdirectory, gzip=True, logger=logger)
                else:
                    copy_file(source_filepath, destination_directory, gzip=False, logger=logger)
            elif source_filepath.endswith(".gz"):
                copy_file(source_filepath, destination_directory, gzip=False, logger=logger)     
            else:
                copy_file(source_filepath, destination_directory, gzip=True, logger=logger)
            
    # Biosynthetic
    # ------------
    source_directory = os.path.join(opts.veba_directory, "biosynthetic")
    if os.path.exists(source_directory):
        logger.info(f"Gathering biosynthetic results: {source_directory}")
        destination_directory = os.path.join(opts.output_directory, "biosynthetic")
        os.makedirs(destination_directory, exist_ok=True)
        for source_filepath in tqdm(glob.glob(os.path.join(opts.veba_directory, "cluster","output", "global", "*")), source_directory):
            copy_file(source_filepath, destination_directory, logger=logger)
            
    # Phylogeny
    # ---------
    source_directory = os.path.join(opts.veba_directory, "phylogeny")
    if os.path.exists(source_directory):
        logger.info(f"Gathering phylogeny results: {source_directory}")
        destination_directory = os.path.join(opts.output_directory, "phylogeny")
        os.makedirs(destination_directory, exist_ok=True)
        for source_filepath in tqdm(glob.glob(os.path.join(opts.veba_directory, "phylogeny","output", "*")), source_directory):
            copy_file(source_filepath, destination_directory, logger=logger)
            
    # Classify
    # --------
    parent_directory = os.path.join(opts.veba_directory, "classify")
    if os.path.exists(parent_directory):
        logger.info(f"Gathering classify results: {parent_directory}")
        
        # Genomes
        for source_directory in glob.glob(os.path.join(opts.veba_directory, "classify", "*")):
            source_directory = os.path.normpath(source_directory)
            if os.path.isdir(source_directory):
                # if source_directory.endswith("/"):
                #     source_directory = source_directory[:-1]
                organism_specific_classify_directory = source_directory
                organism_type = organism_specific_classify_directory.split("/")[-1]
                logger.info(f"Gathering classification files [{organism_type}]: {source_directory}")
                destination_directory = os.path.join(opts.output_directory, "taxonomy_classification", organism_type)
                os.makedirs(destination_directory, exist_ok=True)
                for source_filepath in tqdm(glob.glob(os.path.join(organism_specific_classify_directory, "output",  "*")), organism_specific_classify_directory):
                    if source_filepath.endswith((".gz", ".html")):
                        copy_file(source_filepath, destination_directory, gzip=False, logger=logger)
                    else:
                        copy_file(source_filepath, destination_directory, gzip=True, logger=logger)
                
    # Assembly
    # --------           
    if opts.include_scaffolds:
        source_directory = os.path.join(opts.veba_directory, "assembly")
        if os.path.exists(source_directory):
            logger.info(f"Gathering assembly results [scaffolds]: {source_directory}")
            destination_directory = os.path.join(opts.output_directory, "assembly")
            os.makedirs(destination_directory, exist_ok=True)
            for source_filepath in tqdm(glob.glob(os.path.join(opts.veba_directory, "assembly","*", "output", "scaffolds.fasta")), source_directory):
                id_sample = source_filepath.split("/")[-3]
                destination_filepath = os.path.join(destination_directory, f"{id_sample}.fasta.gz")
                gzip_file(source_filepath, destination_filepath, logger=logger)
                    
    # BAM
    # --------           
    if opts.include_bam:
        source_directory = os.path.join(opts.veba_directory, "assembly")
        if os.path.exists(source_directory):
            logger.info(f"Gathering assembly results [BAM]: {source_directory}")
            destination_directory = os.path.join(opts.output_directory, "assembly")
            os.makedirs(destination_directory, exist_ok=True)
            for source_filepath in tqdm(glob.glob(os.path.join(opts.veba_directory, "assembly","*", "output", "mapped.sorted.bam")), source_directory):
                id_sample = source_filepath.split("/")[-3]
                destination_filepath = os.path.join(destination_directory, f"{id_sample}.bam")
                copy_file(source_filepath, destination_filepath, logger=logger)
                
    # Preprocess
    # --------           
    if opts.include_fastq:
        source_directory = os.path.join(opts.veba_directory, "preprocess")
        if os.path.exists(source_directory):
            logger.info(f"Gathering preprocess results [Fastq]: {source_directory}")
            destination_directory = os.path.join(opts.output_directory, "preprocess")
            os.makedirs(destination_directory, exist_ok=True)
            sequencing_type, preprocess_type = determine_preprocess_type(source_directory)
            
            # Cleaned 
            if sequencing_type == "paired":
                if preprocess_type == "cleaned":            
                    for source_filepath in tqdm(glob.glob(os.path.join(opts.veba_directory, "preprocess","*", "output", "cleaned_1.fastq.gz")), source_directory):
                        id_sample = source_filepath.split("/")[-3]
                        destination_filepath = os.path.join(destination_directory, f"{id_sample}_1.fastq.gz")
                        copy_file(source_filepath, destination_filepath, logger=logger)
                        
                    for source_filepath in tqdm(glob.glob(os.path.join(opts.veba_directory, "preprocess","*", "output", "cleaned_2.fastq.gz")), source_directory):
                        id_sample = source_filepath.split("/")[-3]
                        destination_filepath = os.path.join(destination_directory, f"{id_sample}_2.fastq.gz")
                        copy_file(source_filepath, destination_filepath, logger=logger)
                elif preprocess_type == "trimmed":            
                    for source_filepath in tqdm(glob.glob(os.path.join(opts.veba_directory, "preprocess","*", "output", "trimmed_1.fastq.gz")), source_directory):
                        id_sample = source_filepath.split("/")[-3]
                        destination_filepath = os.path.join(destination_directory, f"{id_sample}_1.fastq.gz")
                        copy_file(source_filepath, destination_filepath, logger=logger)
                        
                    for source_filepath in tqdm(glob.glob(os.path.join(opts.veba_directory, "preprocess","*", "output", "trimmed_2.fastq.gz")), source_directory):
                        id_sample = source_filepath.split("/")[-3]
                        destination_filepath = os.path.join(destination_directory, f"{id_sample}_2.fastq.gz")
                        copy_file(source_filepath, destination_filepath, logger=logger)
            elif sequencing_type == "long":
                if preprocess_type == "cleaned":            
                    for source_filepath in tqdm(glob.glob(os.path.join(opts.veba_directory, "preprocess","*", "output", "cleaned.fastq.gz")), source_directory):
                        id_sample = source_filepath.split("/")[-3]
                        destination_filepath = os.path.join(destination_directory, f"{id_sample}.fastq.gz")
                        copy_file(source_filepath, destination_filepath, logger=logger)
                        
                elif preprocess_type == "trimmed":            
                    for source_filepath in tqdm(glob.glob(os.path.join(opts.veba_directory, "preprocess","*", "output", "trimmed.fastq.gz")), source_directory):
                        id_sample = source_filepath.split("/")[-3]
                        destination_filepath = os.path.join(destination_directory, f"{id_sample}.fastq.gz")
                        copy_file(source_filepath, destination_filepath, logger=logger)
                        

if __name__ == "__main__":
    main()

