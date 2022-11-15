#!/usr/bin/env python
import sys, os, glob, argparse 
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.03.10"

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
        data = pd.read_json(fp, typ="series")
        data["dataset_name"] = data["dataset"].split("/")[-1]
        generic_output[id_genome] = data
    df_generic = pd.DataFrame(generic_output).T
    df_generic.columns = df_generic.columns.map(lambda x: ("generic", x))

    # Get specfic json files
    specific_output = dict()
    for fp in glob.glob(os.path.join(opts.busco_directory, "*", "short_summary.specific.*.json")):
        id_genome = fp.split("/")[-2]
        data = pd.read_json(fp, typ="series")
        data["dataset_name"] = data["dataset"].split("/")[-1]
        specific_output[id_genome] = data
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
        df_output.to_json(opts.json_output)

    df_output.to_csv(opts.output, sep="\t")

    # Run Info
    # runinfo_columns = ["input_file", "mode"]

    # df_info = pd.DataFrame()
    # if not df_generic.empty:
    #     df_info = df_generic[runinfo_columns]

    # df_info.columns = df_info.columns.map(lambda x: ("run_info", x))


    # Remove run info from generic
    # if not df_generic.empty:
        # df_generic = df_generic.drop(runinfo_columns, axis=1)

        # Add multiindex to generic





    



if __name__ == "__main__":
    main()
    
                

