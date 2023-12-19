#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, hashlib
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser


# Soothsayer Ecosystem
from soothsayer_utils import get_file_object, pv, assert_acceptable_arguments


pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.12.13"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__

    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <record_table> -o [id_sample]/concatenated.fa".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", default="stdin", type=str, help = "[id_sample]<tab>[path/to/reference.fa] [Default: stdin]")
    parser.add_argument("-o","--output_directory", default="concatenated_output", type=str, help = "Concatenated output directory [Default: concatenated_output]")
    parser.add_argument("-m", "--minimum_contig_length", type=int, default=1, help="Minimum contig length. [Default: 1] ")
    parser.add_argument("-x", "--extension", type=str, default="fa", help="Concatenated fasta output extension [Default: fa] ")
    parser.add_argument("-b", "--basename", type=str, default="concatenated", help="Concatenated fasta output extension [Default: concatenated] ")
    parser.add_argument("-M", "--mode", type=str, default="infer", help="Concatenate all references with global and build index or build index for each reference {global, local, infer}")
    parser.add_argument("--no_sort", action="store_true", help = "Don't sort the grouped filepaths")
    parser.add_argument("--no_subdirectory", action="store_true", help = "Don't create a nested directory structure")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Input 
    f_in = opts.input
    if opts.input == "stdin":
        f_in = sys.stdin 

    # Output
    os.makedirs(opts.output_directory, exist_ok=True)

    # Read table
    df = pd.read_csv(f_in, sep="\t", index_col=0, header=None)
    m = df.shape[1]

    if opts.mode == "infer":
        assert_acceptable_arguments(m, {0,1})
        opts.mode = {0:"global", 1:"local"}[m]
    
    if opts.mode == "global":
        assert m == 0, "There should be only one column if mode='global'"
    if opts.mode == "local":
        assert m == 1, "There should be two columns if mode='local'"
    
    assert_acceptable_arguments(opts.mode, {"local", "global"})

    if opts.mode == "local":
        # GroupBy and parse
        for id_sample, filepaths in df.groupby(df.index):
            # Get filepaths
            filepaths = list(filepaths.values.ravel())
            if not opts.no_sort:
                filepaths = sorted(filepaths)
                
            # Create output files (and directories)
            if opts.no_subdirectory:
                f_out = get_file_object(
                    path=os.path.join(opts.output_directory, "{}.{}".format(id_sample, opts.extension)),
                    mode="write", 
                    safe_mode=False, 
                    verbose=False,
                )

                saf_filepath = os.path.join(opts.output_directory, "{}.saf".format(id_sample))

                f_duplicates = get_file_object(
                    path=os.path.join(opts.output_directory, "{}.duplicates_removed.list".format(id_sample)),
                    mode="write", 
                    safe_mode=False, 
                    verbose=False,
                )
            else:
                os.makedirs(os.path.join(opts.output_directory, id_sample), exist_ok=True)

                f_out = get_file_object(
                    path=os.path.join(opts.output_directory, id_sample, "{}.{}".format(opts.basename, opts.extension)),
                    mode="write", 
                    safe_mode=False, 
                    verbose=False,
                )

                saf_filepath = os.path.join(opts.output_directory, id_sample, "{}.saf".format(opts.basename))

                f_duplicates = get_file_object(
                    path=os.path.join(opts.output_directory, id_sample, "{}.duplicates_removed.list".format(opts.basename)),
                    mode="write", 
                    safe_mode=False, 
                    verbose=False,
                )
            
            # Read input fasta, filter out short sequences, and write to concatenated file
            sequence_hashes = set()
            saf_data = list()
            for fp in pv(filepaths, description=id_sample, unit= " files"):
                f_query = get_file_object(fp, mode="read", verbose=False)
                for id, seq in SimpleFastaParser(f_query):
                    if len(seq) >= opts.minimum_contig_length:
                        id_hash = hashlib.md5(seq.upper().encode()).hexdigest()
                        id_record = id.split(" ")[0]
                        if id_hash not in sequence_hashes:
                            print(">{}\n{}".format(id, seq), file=f_out)
                            fields = [
                                id_record, 
                                id_record, 
                                1, 
                                len(seq),
                                "+",
                            ]
                            saf_data.append(fields)
                            sequence_hashes.add(id_hash)
                        else:
                            print(id_record, file=f_duplicates)

                f_query.close()

            f_out.close()
            f_duplicates.close()

            df_saf = pd.DataFrame(saf_data, columns=["GeneID", "Chr", "Start", "End", "Strand"])
            df_saf.to_csv(saf_filepath, sep="\t", index=None)

    if opts.mode == "global":
        # GroupBy and parse
        filepaths = df.index 
        if not opts.no_sort:
            filepaths = sorted(filepaths)

        # Create output files (and directories)
        f_out = get_file_object(
            path=os.path.join(opts.output_directory, "{}.{}".format(opts.basename, opts.extension)),
            mode="write", 
            safe_mode=False, 
            verbose=False,
        )

        saf_filepath = os.path.join(opts.output_directory, "{}.saf".format(opts.basename))

        f_duplicates = get_file_object(
            path=os.path.join(opts.output_directory, "{}.duplicates_removed.list".format(opts.basename)),
            mode="write", 
            safe_mode=False, 
            verbose=False,
        )

        # Read input fasta, filter out short sequences, and write to concatenated file
        sequence_hashes = set()
        saf_data = list()
        for fp in pv(filepaths, unit= " files"):
            f_query = get_file_object(fp, mode="read", verbose=False)
            for id, seq in SimpleFastaParser(f_query):
                if len(seq) >= opts.minimum_contig_length:
                    id_hash = hashlib.md5(seq.upper().encode()).hexdigest()
                    id_record = id.split(" ")[0]
                    if id_hash not in sequence_hashes:
                        print(">{}\n{}".format(id, seq), file=f_out)
                        fields = [
                            id_record, 
                            id_record, 
                            1, 
                            len(seq),
                            "+",
                        ]
                        saf_data.append(fields)
                    else:
                        print(id_record, file=f_duplicates)
            f_query.close()

        f_out.close()
        f_duplicates.close()

        df_saf = pd.DataFrame(saf_data, columns=["GeneID", "Chr", "Start", "End", "Strand"])
        df_saf.to_csv(saf_filepath, sep="\t", index=None)




if __name__ == "__main__":
    main()
