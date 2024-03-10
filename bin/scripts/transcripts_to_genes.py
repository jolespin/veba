#!/usr/bin/env python
import sys, os, glob, argparse, gzip
from collections import OrderedDict
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.2.20"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <binning_directory> -x <extension> --sep '\t' > <output.tsv>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--transcript_fasta", type=str,default="stdin", help = "path/to/transcripts.fasta. Assumes rnaSPAdes format.")
    parser.add_argument("-d", "--sep", type=str, default="\t",  help = "Seperator [Default: <tab>")
    parser.add_argument("-t", "--transcript_column_name", type=str, default="id_transcript", help="Transcript column name [Default: id_transcript]")
    parser.add_argument("-g", "--gene_column_name", type=str, default="id_gene", help="Bin column name [Default: id_gene]")
    parser.add_argument("-c", "--column_order", type=str, default="transcript,gene", help="Column order.  Specify either 'transcript,gene' or 'gene,transcript' [Default:transcript,gene]")
    parser.add_argument("-p", "--gene_prefix", type=str,  default="g", help="Bin prefix. Default is to not have a prefix.")
    parser.add_argument("-P", "--preserve_order", action="store_true", help="Preserve ordering")
    parser.add_argument("--header", action="store_true", help="Specify if header should be in output")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Parse
    assert opts.column_order in {"transcript,gene", "gene,transcript"}, "Must choose either 'transcript,gene' or 'gene,transcript' for --column_order"

    transcript_to_gene = OrderedDict()
    if opts.transcript_fasta.endswith(".gz"):
        f = gzip.open(opts.transcript_fasta, "rt")
    else:
        f = open(opts.transcript_fasta, "r")
    for line in f:
        if line.startswith(">"):
            id_transcript = line.strip()[1:].split(" ")[0]
            gene_number = id_transcript.split("_")[-2][1:]
            id_gene = "{}{}".format(opts.gene_prefix, gene_number)
            transcript_to_gene[id_transcript] = id_gene 
    f.close()
        
    df = pd.Series(transcript_to_gene).to_frame(opts.gene_column_name)

    if not opts.preserve_order:
        # index = df.index
        # if opts.column_order == "gene,transcript":
        #     length_of_prefix = len(opts.gene_prefix)
        #     index = sorted(df.index, key=lambda x: int(x[length_of_prefix:]))
        # if opts.column_order == "transcript,gene":
        def f(x):
            id_gene, id_isoform = x.split("_")[-2:]
            return (int(id_gene[1:]), int(id_isoform[1:]))
        index = sorted(df.index, key=f)
        df = df.loc[index,:]



    df.index.name = opts.transcript_column_name 
    df = df.reset_index()

    if opts.column_order == "gene,transcript":
        df = df.iloc[:,[1,0]]
        
    df.to_csv(sys.stdout, sep=opts.sep, header=opts.header, index=None)

if __name__ == "__main__":
    main()
    
                

