#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse
import subprocess
import pandas as pd
import numpy as np
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.02.02"




def main(argv=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <input> -o <output_directory>".format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-i","--input", required=True, type=str, help = "Tab-seperated table with the following format: [id_protein]<tab>[id_protein_cluster]<tab>[annotation]")
    parser.add_argument("-o","--output_directory", type=str, default="consensus_orthogroup_annotations", help = "Output directory [Default: consensus_orthogroup_annotations]")
    parser.add_argument("--similarity_threshold", type=float, default=0.8, help = "Threshold for similarity analysis [Default: 0.8]")
    parser.add_argument("--retain_hypothetical", type=int, default=1, help = "Consider hypothetical annotations as being annotations [Default: 1]")
    parser.add_argument("--retain_unannotated", type=int, default=1, help = "Consider unannotations (i.e., blank functions) in the scording system [Default: 1]")
    parser.add_argument("--unannotated_weight", type=float, default=0.382, help = "Weight for unannotations (i.e., blank functions) in the scording system? [Default: 0.382]")
    parser.add_argument("--representative_threshold", type=float, default=0.618, help = "Score to consider a function representative [Default: 0.618]")
    parser.add_argument("--unifunc", type=str, default="CONDA_PREFIX", help = "path/to/unifunc [Default: CONDA_PREFIX]")
    parser.add_argument("--orthogroup_label", type=str, default="id_protein_cluster", help = "orthogroup label [Default: id_protein_cluster]")

    # Options
    opts = parser.parse_args(argv)

    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # if opts.input == "stdin": # Can't handle it yet
    #     opts.input = sys.stdin
    os.makedirs(opts.output_directory, exist_ok=True)

    if opts.unifunc == "CONDA_PREFIX":
        opts.unifunc = os.path.join(os.environ["CONDA_PREFIX"], "bin", "unifunc")
    os.environ["unifunc"] = opts.unifunc

    # Arguments
    args = [
        os.environ["unifunc"],
        "cluster_function",
        "-o {}".format(opts.output_directory),
        "-st {}".format(opts.similarity_threshold),
        "-uw {}".format(opts.unannotated_weight),
        "-rt {}".format(opts.representative_threshold),
        "-owr",
    ]

    if bool(opts.retain_hypothetical):
        args.append("-kh")
    if not bool(opts.retain_unannotated):
        args.append("-ra")
    if opts.input != sys.stdin:
        args.append("-i {}".format(opts.input))

    df_input = pd.read_csv(opts.input, sep="\t", header=None, index_col=0)
    print("Input: {}".format(opts.input), file=sys.stdout)
    print(" * Number of ORFs: {}".format(df_input.shape[0]), file=sys.stdout)
    print(" * Number of orthogroups: {}".format(df_input.iloc[:,0].nunique()), file=sys.stdout)
    print(" * Number of unique annotations: {}".format(df_input.iloc[:,1].nunique()), file=sys.stdout)
    print("", file=sys.stdout)

    # Run unifunc
    cmd = " ".join(args)
    print("Running UniFunc Command:\n{}".format(cmd), file=sys.stdout)
    # if opts.input == sys.stdin:
    #     process = subprocess.Popen([cmd], stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    #     stdout, stderr = process.communicate(input=sys.stdin.buffer.read())
    # else:
    process = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()

    returncode = process.returncode

    assert returncode == 0, stderr

    print("", file=sys.stdout)

    with open(os.path.join(opts.output_directory, "rep_func_cluster.tsv"), "r") as f:
        orthogroup_annotations = list()
        for line in tqdm(f, "Parsing annotations", unit=" Orthogroups"):
            line = line.strip()
            fields = line.split("\t")
            # id_orthogroup = fields[0]
            if len(fields) == 1:
                fields += [np.nan, np.nan]
            
            orthogroup_annotations.append(fields)
        df_output = pd.DataFrame(orthogroup_annotations, columns=[opts.orthogroup_label, "score", "annotation"])
        df_output = df_output.sort_values([opts.orthogroup_label, "score"], ascending=[True, False]).set_index(opts.orthogroup_label)
        df_output = df_output[~df_output.index.duplicated(keep='first')]
        df_output.to_csv(os.path.join(opts.output_directory, "consensus_annotations.tsv"), sep="\t")
        
        print("", file=sys.stdout)
        print("Output: {}".format(os.path.join(opts.output_directory, "consensus_annotations.tsv")), file=sys.stdout)
        print(" * Number of annotated orthogroups: {}".format(df_output.dropna(how="any").shape[0]), file=sys.stdout)
        print(" * Number of unique annotations: {}".format(df_output.dropna(how="any")["annotation"].nunique()), file=sys.stdout)


if __name__ == "__main__":
    main()
