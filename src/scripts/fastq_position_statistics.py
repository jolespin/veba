#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, gzip
import numpy as np
import pandas as pd
from tqdm import tqdm
# import pyfastx
from Bio.SeqIO.QualityIO import FastqGeneralIterator

pd.options.display.max_colwidth = 100
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.10.24"

def statistics(fp, phred):
    # Build table
    quality_table = list()
    failed_reads = list()
    
    # for id,seq,quality in tqdm(pyfastx.Fastq(fp, build_index=False), desc="Reading fastq filepath: {}".format(fp), unit=" read"):
    file_open_function = {True:gzip.open, False:open}[fp.endswith(".gz")]
    with file_open_function(fp, "rt") as f:
        for id,seq,quality in tqdm(FastqGeneralIterator(f), desc="Reading fastq filepath: {}".format(fp), unit=" read"):
            try:
                quality = np.asarray(list(map(lambda q: ord(q) - phred, quality)))#.astype(int)
                quality_table.append(quality)
            except UnicodeDecodeError as e:
                failed_reads.append(id)


    quality_table = np.stack(quality_table)

    df_minmeanmax = pd.DataFrame([
        pd.Series(np.min(quality_table, axis=0), name="min"),
        pd.Series(np.mean(quality_table, axis=0), name="mean"),
        pd.Series(np.max(quality_table, axis=0), name="max"),
    ]).T
    df_quantiles = pd.DataFrame(np.quantile(quality_table, q=[0.25,0.5,0.75], axis=0).T, columns = ["q=0.25", "q=0.5", "q=0.75"])
    df_output = pd.concat([df_minmeanmax, df_quantiles], axis=1)
    df_output.index = df_output.index.values + 1
    df_output.index.name = "position"

    return df_output

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = """
    Single fastq file and output using -o argument: {} -o output.tsv file.fastq.gz\n
    Multiple fastq files, only include basename, and only include a certain field per file: {} -b -f q=0.25 -o output.tsv file_1.fq file_2.fq,...file_n.fq]\n
    Multiple fastq files and output to stdout: {}  file_1.fastq file_2.fq,...file_n.fq.gz > output.tsv
    """.format(__program__, __program__, __program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("fastq_filepath",  type=str, nargs="+", help = "path/to/fastq[.gz] file(s)")
    parser.add_argument("-o","--output", default="stdout", type=str, help = "Output filepath [Default: stdout]")
    parser.add_argument("-f","--field",  type=str, help = "Field when batch [min, mean, max, q=0.25, q=0.5, q=0.75]")
    parser.add_argument("-b","--basename",  action="store_true", help = "Output basename for multiple fastq files")
    parser.add_argument("--phred",  type=int, default=33, help = "Phred offset [Default: 33]")

    # parser.add_argument("-r","--retain_index",  action="store_true", help = "Keep fastq index created by `pyfastx`")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Checks 
    assert opts.phred in {33,64}, "--phred must be either 33 or 64. The following is invalid: {}".format(opts.phred)

    # Build stats table for single fastq file
    if len(opts.fastq_filepath) == 1:
        df_output = statistics(opts.fastq_filepath[0], phred=opts.phred)
    else:
        # Build stats table for multiple fastq files but just pull out one column
        if opts.field:
            assert opts.field in {"min", "mean", "max", "q=0.25", "q=0.5", "q=0.75"}, "--field value is not supported.  Please choose between [min, mean, max, q=0.25, q=0.5, q=0.75] or none at all."
            field_table = dict()
            for fp in opts.fastq_filepath:
                name = fp
                if opts.basename:
                    name = fp.split("/")[-1]
                field_table[name] = statistics(fp, phred=opts.phred)[opts.field]
            df_output = pd.DataFrame(field_table)

        # Build stats table for multiple fastq files but include all columns (results in multiindex)
        else:
            dataframes = list()
            for fp in opts.fastq_filepath:
                name = fp
                if opts.basename:
                    name = fp.split("/")[-1]
                df = statistics(fp, phred=opts.phred)
                df.columns = df.columns.map(lambda x: (name, x))
                dataframes.append(df)
            df_output = pd.concat(dataframes, axis=1)
        
    # Output table
    if opts.output == "stdout":
        opts.output = sys.stdout 
        
    df_output.to_csv(opts.output, sep="\t")

if __name__ == "__main__":
    main()
