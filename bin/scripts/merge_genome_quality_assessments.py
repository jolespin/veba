#!/usr/bin/env python
import sys, os, glob, argparse 
from collections import defaultdict
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.11.9"

def get_prokaryotic_description(x, fields=["Completeness_Model_Used", "Additional_Notes"]):
    output = list()
    for k, v in x.reindex(fields).items():
        if pd.notnull(v):
            if k == "Completeness_Model_Used":
                k = "checkm2_model"
            if k == "Additional_Notes":
                k = "warnings"
            label = "[{}={}]".format(k,v)
            output.append(label)
    return "".join(output)

def get_viral_description(x, fields=["checkv_quality", "miuvig_quality", "aai_confidence", "hmm_completeness_upper", "virus_score"]):
    output = list()
    for k, v in x.reindex(fields).items():
        if pd.notnull(v):
            label = "[{}={}]".format(k,v)
            output.append(label)
    return "".join(output)

def get_eukaryotic_description(x, fields=["one_line_summary", "dataset_name"]):
    output = list()
    for k, v in x.reindex(fields).items():
        if pd.notnull(v):
            if k == "one_line_summary":
                k = "busco_notation"
            if k == "dataset_name":
                k = "busco_marker_set"
            label = "[{}={}]".format(k,v)
            output.append(label)
    return "".join(output)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <binning_directory> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--binning_directory", type=str, default="veba_output/binning", help = "path/to/binning_directory.  Assumes directory structure is [binning_directory]/[domain]/[sample]/output/[quality].tsv [Default: veba_output/binning]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/merge_genome_quality.tsv [Default: stdout]")
    parser.add_argument("-V","--viral_subdirectory_name", type=str, default="viral", help = "Viral subdirectory label in --binning_directory [Default: viral]")
    parser.add_argument("-P","--prokaryotic_subdirectory_name", type=str, default="prokaryotic", help = "Prokaryotic subdirectory label in --binning_directory [Default: prokaryotic]")
    parser.add_argument("-E","--eukaryotic_subdirectory_name", type=str, default="eukaryotic", help = "Eukaryotic subdirectory label in --binning_directory [Default: eukaryotic]")
    parser.add_argument("--no_header", action="store_true", help = "Don't write header")



    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Output
    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Concatenated    
    genome_quality_dataframes = list()

    # Prokaryotic
    prokaryotic_genome_quality_files = glob.glob(os.path.join(opts.binning_directory, opts.prokaryotic_subdirectory_name, "*", "output", "checkm2_results.filtered.tsv"))
    if prokaryotic_genome_quality_files:
        print("* Compiling prokaryotic genome quality from following files:", *prokaryotic_genome_quality_files, sep="\n    ", file=sys.stdout)

        for fp in prokaryotic_genome_quality_files:
            id_domain = fp.split("/")[-4]
            id_sample = fp.split("/")[-3]
            df_checkm2 = pd.read_csv(fp, sep="\t", index_col=0)
            df_unified = df_checkm2.loc[:,["Completeness", "Contamination"]]
            df_unified["description"] = df_checkm2.apply(get_prokaryotic_description, axis=1)
            df_unified.columns = ["completeness", "contamination", "description"]
            df_unified.insert(0, "domain", id_domain)
            df_unified.insert(1, "id_sample", id_sample)
            df_unified.index.name = "id_genome"
            genome_quality_dataframes.append(df_unified)

    else:
        print("Could not find any prokaryotic genome assessment tables from CheckM2 in the following directory: {}".format(opts.binning_directory), file=sys.stdout)

    # Viral
    viral_genome_quality_files = glob.glob(os.path.join(opts.binning_directory, opts.viral_subdirectory_name, "*", "output", "checkv_results.filtered.tsv"))
    if viral_genome_quality_files:
        print("* Compiling viral genome quality from following files:", *viral_genome_quality_files, sep="\n    ", file=sys.stdout)

        for fp in viral_genome_quality_files:
            id_domain = fp.split("/")[-4]
            id_sample = fp.split("/")[-3]
            df_checkv = pd.read_csv(fp, sep="\t", index_col=0)
            df_unified = df_checkv.loc[:,["completeness", "contamination"]]
            df_unified["description"] = df_checkv.apply(get_viral_description, axis=1)
            df_unified.columns = ["completeness", "contamination", "description"]
            df_unified.insert(0, "domain", id_domain)
            df_unified.insert(1, "id_sample", id_sample)
            df_unified.index.name = "id_genome"
            genome_quality_dataframes.append(df_unified)

    else:
        print("Could not find any viral genome assessment tables from CheckV in the following directory: {}".format(opts.binning_directory), file=sys.stdout)

    # Eukaryotic
    eukaryotic_genome_quality_files = glob.glob(os.path.join(opts.binning_directory, opts.eukaryotic_subdirectory_name, "*", "output", "busco_results.filtered.tsv"))
    if eukaryotic_genome_quality_files:
        print("* Compiling eukaryotic genome quality from following files:", *eukaryotic_genome_quality_files, sep="\n    ", file=sys.stdout)

        for fp in eukaryotic_genome_quality_files:
            id_domain = fp.split("/")[-4]
            id_sample = fp.split("/")[-3]
            df_busco = pd.read_csv(fp, sep="\t", index_col=0, header=[0,1])
            use_generic = True
            if "specific" in df_busco.columns.get_level_values(0):
                if not df_busco[("specific", "Complete")].dropna().empty:
                    use_generic = False 
            if use_generic:
                df_busco = df_busco["generic"]
            else:
                df_busco = df_busco["specific"]
            
            df_unified = df_busco.loc[:,["Complete", "Multi copy"]]
            df_unified["description"] = df_busco.apply(get_eukaryotic_description, axis=1)

            df_unified.columns = ["completeness", "contamination", "description"]
            df_unified.insert(0, "domain", id_domain)
            df_unified.insert(1, "id_sample", id_sample)
            df_unified.index.name = "id_genome"
            genome_quality_dataframes.append(df_unified)

    else:
        print("Could not find any eukaryotic genome assessment tables from BUSCO in the following directory: {}".format(opts.binning_directory), file=sys.stdout)

    assert len(genome_quality_dataframes) > 0, "No genome quality assessment tables were found in the following directory: {}".format(opts.binning_directory)
    df_genome_quality = pd.concat(genome_quality_dataframes, axis=0)
    df_genome_quality.to_csv(opts.output, sep="\t", header=not bool(opts.no_header))


if __name__ == "__main__":
    main()
    
                

