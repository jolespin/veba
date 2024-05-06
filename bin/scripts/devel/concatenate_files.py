#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, shutil, gzip, warnings
from collections import Counter
from tqdm import tqdm 

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.4.30"


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
    parser.add_argument("-o","--output", type=str,  default="stdout", help = "Output filepath [Default: stdout]")
    parser.add_argument("-u", "--no_duplicate_files", action="store_true", help = "Assert there are no duplicates")
    parser.add_argument("-m", "--mixed_decompressed_and_compressed_types", action="store_true", help = "Mixed gzip and decompressed files")

    parser_input = parser.add_mutually_exclusive_group(required=True)
    parser_input.add_argument("-f", "--filepaths", type=str, nargs="+", help="Filepaths of files to concatenate")
    parser_input.add_argument("-g", "--glob_pattern", type=str, help="Glob pattern (e.g., path/to/*.tsv)")
    parser_input.add_argument("-l", "--list", type=str, help="Filepath of file containing filepaths of files to be concatenated with each filepath on a separate line")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Output
    if opts.output == "stdout":
        f_output = sys.stdout
    else:
        f_output = open(opts.output, "w")

    # --filepaths
    if opts.filepaths:
        filepaths = list(filter(bool, map(lambda fp: fp.strip(), opts.filepaths)))
        assert len(filepaths) > 0, "Cleaning for --filepaths '{}' returns 0 files".format(opts.filepaths)


    # --glob_pattern
    if opts.glob_pattern:
        filepaths = glob.glob(opts.glob_pattern)
        assert len(filepaths) > 0, "Glob pattern: '{}' returns 0 files".format(opts.glob_pattern)

    # --list
    if opts.list:
        filepaths = list()
        if opts.list.endswith(".gz"):
            f_list = gzip.open(opts.list, "rt")
        else:
            f_list = open(opts.list, "r")
        for line in f_list:
            line = line.strip()
            if line:
                filepaths.append(line)
        f_list.close()

    if opts.no_duplicate_files:
        c = Counter(filepaths)
        most_common_filepath = c.most_common()[0]
        assert most_common_filepath[1] < 2, "The following filepath {} is represented {} times.  To allow this, do not select --no_duplicate_files".format(*most_common_filepath)


    # Copy filepaths into new file
    # Mixed file types
    if opts.mixed_decompressed_and_compressed_types:
        for fp in tqdm(filepaths, total=len(filepaths), desc="Concatenting files"):
            if fp.endswith(".gz"):
                f = gzip.open(fp, "rt")
            else:
                f = open(fp, "r")
            shutil.copyfileobj(f, f_output)
            f.close()

    # Not mixed file types
    else:
        number_of_gzipped_filepaths = 0
        for fp in filepaths:
            if fp.endswith(".gz"):
                number_of_gzipped_filepaths += 1
        if all([
            number_of_gzipped_filepaths > 0,
            number_of_gzipped_filepaths != len(filepaths),
            ]):
            warnings.warn("Mixed compressed and decompressed files exist which will result in a corrupted file.  If this is not desired, please select --mixed_decompressed_and_compressed_types")
        for fp in tqdm(filepaths, total=len(filepaths), desc="Concatenting files"):
            with open(fp,'rb') as f:
                shutil.copyfileobj(f, f_output)

    # Output
    if f_output != sys.stdout:
        f_output.close()
        
if __name__ == "__main__":
    main()

# https://stackoverflow.com/a/27077437/678572
# with open('output_file.txt','wb') as wfd:
#     for f in ['seg1.txt','seg2.txt','seg3.txt']:
#         with open(f,'rb') as fd:
#             shutil.copyfileobj(fd, wfd)