#!/usr/bin/env python
import sys, os, glob, argparse 
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.7.5"

def read_busco_json(path):
    json_data = pd.read_json(path, typ="series")[["results", "lineage_dataset"]]
    results = pd.Series(json_data["results"]).loc[["one_line_summary", "Complete", "Single copy", "Multi copy", "Fragmented", "Missing", "n_markers"]]
    lineage_dataset = pd.Series(json_data["lineage_dataset"]).loc[["name", "creation_date", "number_of_buscos", "number_of_species"]]
    lineage_dataset.index = ["dataset_name", "creation_date", "number_of_busco_markers", "number_of_species"] 
    return pd.concat([results, lineage_dataset])

def main(argv=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -a <gtdbtk.ar122.summary.tsv> -b <gtdbtk.bac120.summary.tsv> -o <output.tsv>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--busco_directory", required=True, type=str, help = "path/to/gtdbtk.ar122.summary.tsv")
    parser.add_argument("-j","--json_output", type=str, help = "path/to/merged_busco.json")
    parser.add_argument("-o","--output", type=str,  default="stdout", help = "Output merged multiple sequence alignment [Default: stdout]")

    # Options
    opts = parser.parse_args(argv)
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Open output file
    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Get generic json files
    generic_output = dict()
    for fp in glob.glob(os.path.join(opts.busco_directory, "*", "short_summary.generic.*.json")):
        id_genome = fp.split("/")[-2]
        generic_output[id_genome] = read_busco_json(fp)
    df_generic = pd.DataFrame(generic_output).T
    df_generic.columns = df_generic.columns.map(lambda x: ("generic", x))

    # Get specfic json files
    specific_output = dict()
    for fp in glob.glob(os.path.join(opts.busco_directory, "*", "short_summary.specific.*.json")):
        id_genome = fp.split("/")[-2]
        specific_output[id_genome] = read_busco_json(fp)
    df_specific = pd.DataFrame(specific_output).T
    df_specific.columns = df_specific.columns.map(lambda x: ("specific", x))


    # Concatenate tables
    dataframes = list() 
    if not df_generic.empty:
        dataframes.append(df_generic)

    if not df_specific.empty:
        dataframes.append(df_specific)

    if dataframes:
        df_output = pd.concat(dataframes, axis=1)
    else:
        print("Neither generic nor specific json were available. This likely means no markers were found via BUSCO. Please check the BUSCO log files", file=sys.stderr)
        sys.exit(1)
    if not df_output.empty:
        df_output.index.name = "id_genome"
    
    # Output 
    if opts.json_output:
        df_output.T.to_json(opts.json_output)

    df_output.to_csv(opts.output, sep="\t")





    



if __name__ == "__main__":
    main()
    
                

