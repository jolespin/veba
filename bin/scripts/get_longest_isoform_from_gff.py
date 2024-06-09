#!/usr/bin/env python
import sys, os, argparse
from tqdm import tqdm
from collections import defaultdict

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.9.18"

def get_longest_interval(intervals):
    """
    Implementation credits to @yep-d from the following source: 
    * https://stackoverflow.com/questions/77116236/find-longest-interval-between-overlapping-intervals-in-python
    """
    #get sorted list of intervals by starting point
    aux = sorted(intervals, key=lambda x: x[0])
    new_group = [aux[0]]
    output = []
    for interval in aux[1:]:
        if interval[0] > new_group[-1][-1]:
        #if the interval does not overlap with the last interval we start a new group
            output.append(
                #find interval of maximum length
                max(new_group, key=lambda interval: interval[1] - interval[0])
            )
            #start new group
            new_group = [interval]
        else:
            #otherwise just append it onto the same group if it overlaps
            new_group.append(interval)
    #take care of the last group, doing the same as above
    output.append(max(new_group, key=lambda interval: interval[1] - interval[0]))
    return output


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
    parser.add_argument("-i","--input_gff", type=str, default="stdin", help = "path/to/input.gff [Default: stdin]")
    parser.add_argument("-o","--output_gff", type=str, default="stdout", help = "path/to/output.gff [Default: stdout]")
    parser.add_argument("-f","--field", type=str, default = "CDS", help = "Field to dereplicate [Default: CDS]")
    parser.add_argument("--no_comment", action = "store_true", help = "Don't include commented fields in output")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Input
    if opts.input_gff == "stdin":
        f_in = sys.stdin 
    else:
        f_in = open(opts.input_gff, "r")

    # Output
    if opts.output_gff == "stdout":
        f_out = sys.stdout 
    else:
        f_out = open(opts.output_gff, "w")

    # Parse input
    contigstrand_to_locations_to_record = defaultdict(dict)
    if opts.no_comment:
        for line in tqdm(f_in, "Reading input GFF: {}".format(opts.input_gff)):
            line = line.strip()
            if opts.field in line:
                fields = line.split("\t")
                id_contig = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                contigstrand_to_locations_to_record[(id_contig, strand)][(start,end)] = line

    else:
        for line in tqdm(f_in, "Reading input GFF: {}".format(opts.input_gff)):
            line = line.strip()
            if line.startswith("#"):
                print(line, file=f_out)
            if opts.field in line:
                fields = line.split("\t")
                id_contig = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                contigstrand_to_locations_to_record[(id_contig, strand)][(start,end)] = line
    if f_in != sys.stdin:
        f_in.close()

    # Dereplicate output
    records = list()
    for (id_contig, strand), d1 in contigstrand_to_locations_to_record.items():
        longest_intervals = get_longest_interval(d1.keys())
        for interval in longest_intervals:
            start, end = interval
            data = [id_contig, start, end, d1[interval]]
            records.append(data)
    records = sorted(records, key=lambda x: (x[0], x[1], x[2]))

    # Write output
    for (id_contig, start, end, record) in tqdm(records, "Writing output GFF with longest CDS: {}".format(opts.output_gff)):
        print(record, file=f_out)
    if f_out != sys.stdout:
        f_out.close()

if __name__ == "__main__":
    main()
    
                
