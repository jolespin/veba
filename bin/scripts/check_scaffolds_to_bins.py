#!/usr/bin/env python
import sys, os, argparse

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.08.20"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__

    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "S2B_TMP = $({} -i <scaffolds,to,bins> -n <names,of,programs>)".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("-i","--scaffolds_to_bins", required=True, type=str, help = "Comma-separated list of scaffolds_to_bins.tsv files")
    parser.add_argument("-n","--names", required=True, type=str, help = "Comma-separated list of names for the --scaffolds_to_bins")
    parser.add_argument("--sep", default="\t", type=str, help = "Separator used for scaffolds to bins [Default:tab'")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename


    opts.scaffolds_to_bins = opts.scaffolds_to_bins.strip().split(",")
    opts.names = opts.names.strip().split(",")
    a = len(opts.scaffolds_to_bins)
    b = len(opts.names)
    assert a == b, "Number of --scaffolds_to_bins ({}) does not match number of --names ({})".format(a,b)

    output__scaffolds_to_bins = list() 
    output__names = list() 
    for fp, name in zip(opts.scaffolds_to_bins, opts.names):
        keep = False
        if os.path.exists(fp):
            with open(fp, "r") as f:
                for line in f.readlines():
                    line = line.strip()
                    fields = list(filter(bool, line.split(opts.sep)))
                    if len(fields) == 2:
                        keep = True 
                        break

        if keep:
            output__scaffolds_to_bins.append(fp)
            output__names.append(name)
        else:
            print("Discarding empty file: {}".format(fp), file=sys.stderr)
    print(",".join(output__scaffolds_to_bins), ",".join(output__names), sep=" ")
                     



   

if __name__ == "__main__":
    main()
    
                

