#!/usr/bin/env python
import sys, os, argparse
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.9.15"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <components> -s <synopsis> -d <diamond> -o <output.tsv[.gz]>".format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-d","--diamond", type=str, required=True, help = "path/to/homology.tsv.gz")
    parser_io.add_argument("-c","--components", type=str, required=True, help = "path/to/identifier_mapping.components.tsv.gz")
    parser_io.add_argument("-b","--bgcs", type=str, required=True, help = "path/to/identifier_mapping.bgcs.tsv.gz")
    parser_io.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")

    parser_thresholds = parser.add_argument_group('Novelty score threshold arguments')
    parser.add_argument("--pident", type=float, default=0.0, help = "pident lower bound [float:0 ≤ x < 100] [Default: 0]")
    parser.add_argument("--qcovhsp", type=float, default=0.0, help = "qcovhsp lower bound [float:0 ≤ x < 100] [Default: 0]")
    parser.add_argument("--scovhsp", type=float, default=0.0, help = "scovhsp lower bound [float:0 ≤ x < 100] [Default: 0]")
    parser.add_argument("--evalue", type=float, default=1e-3, help = "e-value lower bound [float:0 < x < 1] [Default: 1e-3]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Output
    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Load tables
    df_components = pd.read_csv(opts.components, sep="\t", index_col=None)
    df_bgcs = pd.read_csv(opts.bgcs, sep="\t", index_col=0)
    df_diamond = pd.read_csv(opts.diamond, sep="\t", index_col=0, header=[0,1])


    # Gene -> BGC
    component_to_bgc = dict(zip(df_components["component_id"], df_components["bgc_id"]))

    # MiBIG
    df_mibig = df_diamond["MiBIG"].dropna(how="all", axis=0)
    df_mibig = df_mibig.query("pident >= {}".format(opts.pident))
    df_mibig = df_mibig.query("qcovhsp >= {}".format(opts.qcovhsp))
    df_mibig = df_mibig.query("scovhsp >= {}".format(opts.scovhsp))
    df_mibig = df_mibig.query("evalue <= {}".format(opts.evalue))

    bgc_to_mibignhits = pd.Series([0]*df_bgcs.shape[0], df_bgcs.index)
    bgc_to_mibignhits.update(df_mibig.index.map(lambda x: component_to_bgc[x]).value_counts())
    df_bgcs["number_of_mibig_hits"] = bgc_to_mibignhits
    df_bgcs["novelty_score"] = 1 - df_bgcs["number_of_mibig_hits"]/df_bgcs["number_of_genes"]

    # VFDB
    df_vfdb = df_diamond["VFDB"].dropna(how="all", axis=0)
    df_vfdb = df_vfdb.query("pident >= {}".format(opts.pident))
    df_vfdb = df_vfdb.query("qcovhsp >= {}".format(opts.qcovhsp))
    df_vfdb = df_vfdb.query("scovhsp >= {}".format(opts.scovhsp))
    df_vfdb = df_vfdb.query("evalue <= {}".format(opts.evalue))

    bgc_to_vfdbnhits = pd.Series([0]*df_bgcs.shape[0], df_bgcs.index)
    bgc_to_vfdbnhits.update(df_vfdb.index.map(lambda x: component_to_bgc[x]).value_counts())
    df_bgcs["number_of_vfdb_hits"] = bgc_to_mibignhits
    df_bgcs["virulence_ratio"] = df_bgcs["number_of_vfdb_hits"]/df_bgcs["number_of_genes"]

    # Output
    df_bgcs.to_csv(opts.output, sep="\t")

if __name__ == "__main__":
    main()
    
                
