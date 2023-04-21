#!/usr/bin/env python
import sys, os, argparse, gzip
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.4.17"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__

    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} file.fasta".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("input", type=str, nargs="*", help = "Fasta files or stdin")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    if not opts.input:
        identifiers = set()
        duplicates = set()
        for line in tqdm(sys.stdin, "stdin"):
            if line.startswith(">"):
                id = line[1:].split(" ")[0].strip()
                if id not in identifiers:
                    identifiers.add(id)
                else:
                    duplicates.add(id)
        if duplicates:
            print("# Duplicates:", *sorted(duplicates), file=sys.stdout, sep="\n", end=None)
            sys.exit(1)
        else:
            sys.exit(0)
    else:
        files_with_duplicates = set()
        for fp in opts.input:
            identifiers = set()
            duplicates = set()
            f = {True:gzip.open(fp, "rt"), False:open(fp, "r")}[fp.endswith(".gz")]
            for line in tqdm(f, fp):
                if line.startswith(">"):
                    id = line[1:].split(" ")[0]
                    if id not in identifiers:
                        identifiers.add(id)
                    else:
                        duplicates.add(id)
            if duplicates:
                files_with_duplicates.add(fp)
                print(f"[Fail] {fp}", file=sys.stdout)
            else:
                print(f"[Pass] {fp}", file=sys.stdout)

        if files_with_duplicates:
            sys.exit(1)
        else:
            sys.exit(0)
                     
   
if __name__ == "__main__":
    main()
    
                

