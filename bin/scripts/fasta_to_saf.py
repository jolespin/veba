#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser


# Soothsayer Ecosystem
from soothsayer_utils import get_file_object, pv

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.04.04"


def fasta_to_saf(path, compression="infer"):
    """
    #     GeneID	Chr	Start	End	Strand
    # http://bioinf.wehi.edu.au/featureCounts/

    # Useful:
    import re
    record_id = "lcl|NC_018632.1_cds_WP_039228897.1_1 [gene=dnaA] [locus_tag=MASE_RS00005] [protein=chromosomal replication initiator protein DnaA] [protein_id=WP_039228897.1] [location=410..2065] [gbkey=CDS]"
    re.search("\[locus_tag=(\w+)\]", record_id).group(1)
    # 'MASE_RS00005'

    """


    saf_data = list()

    if path == "stdin":
        f = sys.stdin
    else:
        f = get_file_object(path, mode="read", compression=compression, verbose=False)
        
    for id_record, seq in pv(SimpleFastaParser(f), "Reading sequences [{}]".format(path)):
        id_record = id_record.split(" ")[0]
        fields = [
            id_record, 
            id_record, 
            1, 
            len(seq),
            "+",
        ]
        saf_data.append(fields)
    if f is not sys.stdin:
        f.close()
    return pd.DataFrame(saf_data, columns=["GeneID", "Chr", "Start", "End", "Strand"])

#
def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <fasta_file> --compression infer > <saf_file>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--fasta_file", default="stdin", type=str, help = "path/to/fasta.fa[.gz] [Default: stdin")
    parser.add_argument("-c","--compression", type=str, default="infer", help = "Compression type {None, gzip, bz2, infer} [Default: infer")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Fasta -> SAF
    df_saf = fasta_to_saf(opts.fasta_file, compression=opts.compression)
    df_saf.to_csv(sys.stdout, sep="\t", index=None)

if __name__ == "__main__":
    main()
