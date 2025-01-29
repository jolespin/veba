#!/usr/bin/env python
import sys,os, argparse, shutil
from datetime import datetime
from tqdm import tqdm
from pyexeggutor import (
    build_logger,
    open_file_reader, 
    open_file_writer, 
)

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "v2024.11.8"

DEFAULT_ATTRIBUTES = [
# General
"ID",
"Parent",
# Gene
"gene_id",
"geneId",
# Transcript
"transcript_id",
"transcriptId",
# Exon
"exon_id",
"exonId",    
# Contig
"contig_id",
"contigId",
# mRNA
"mRNA_id",
"mRNAId",
# Protein
"protein_id",
"proteinId",
]
DEFAULT_ATTRIBUTES = ",".join(DEFAULT_ATTRIBUTES)

now = datetime.now()
timestamp = now.strftime("%Y-%m-%d %H:%M:%S")

def main(args=None):
    # Options
    # =======
    # Path info
    python_executable = sys.executable
    bin_directory = "/".join(python_executable.split("/")[:-1])
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, sys.version.split(" ")[0], python_executable, script_filename)
    usage = f"{__program__} "
    epilog = "Copyright 2024 Josh L. Espinoza (jol.espinoz@gmail.com)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-i", "--input", type=str,  required=True, help = f"path/to/input.gff[.gz]")
    parser_io.add_argument("-p", "--prefix", type=str,  required=True, help = f"Genome name used for prefix followed by some separator [Recommended: Genome identifier (e.g., Monili1) followed by separator (e.g., Monili2__)]")
    parser_io.add_argument("-o", "--output", type=str, help = f"path/to/output.gff[.gz].  [Default: original_path.MODIFIED.gff[.gz]]")
    parser_io.add_argument("-a", "--attributes", type=str, default=DEFAULT_ATTRIBUTES,help = f"Attributes to prepend label separated by commans [Default: {DEFAULT_ATTRIBUTES}]")
    parser_io.add_argument( "--no_contig_prepend", action="store_true", help = f"Do not prepend to contig identifier in first column")
    parser_io.add_argument("--overwrite", action="store_true", help = f"Overwrite the GFF file")
    
    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    # Logger
    logger = build_logger(f"Prepend label to GFF identifiers [{opts.script_filename}]: {opts.prefix}")
    # reset_logger(logger)

    # Commands
    logger.info(f"Command: {sys.argv}")
    
    # Output
    if not opts.output:
        if opts.input.endswith(".gz"):
            fileprefix = ".".join(opts.input.split(".")[:-2])
            opts.output = fileprefix + ".MODIFIED" + ".gff.gz"
        else:
            fileprefix = ".".join(opts.input.split(".")[:-1])
            opts.output = fileprefix + ".MODIFIED" + ".gff"
        logger.info(f"--output was not provided.  Writing (temporary) output: {opts.output}")
        
    # Attributes
    opts.attributes = set(filter(bool, map(str.strip, opts.attributes.split(","))))

    msg = f"Reading input {opts.input} and writing output {opts.output} GFF files"
    logger.info(msg)

    with (
        open_file_reader(opts.input) as f_in,
        open_file_writer(opts.output) as f_out,
        ):
        print("# Modified on {} by {} {}".format(timestamp, __program__, __version__), file=f_out)
        print("# Command: {}".format(sys.argv), file=f_out)

        for i, line in tqdm(enumerate(f_in), desc=msg):
            line = line.strip()
            if line:
                if line.startswith("#"):
                    print(line, file=f_out)
                else:
                    try:
                        id_contig, source, feature, start, end, score, strand, frame, attributes_unformatted = line.split("\t")
                        if not opts.no_contig_prepend:
                            id_contig = "{}{}".format(opts.prefix, id_contig)
                                     
                        attributes = list()
                        for attribute in attributes_unformatted.split(";"):
                            attribute = attribute.strip()
                            if attribute:
                                k, v = attribute.split("=")
                                if k in opts.attributes:
                                    attributes.append("{}={}{}".format(k, opts.prefix, v))
                                else:
                                    attributes.append(attribute)
                        print("\t".join([id_contig, source, feature, start, end, score, strand, frame, ";".join(attributes)]), file=f_out)
                        
                        
                    except ValueError as e:
                        msg = f"{e}\n[Skipping] Please check formatting for the following record: {line}"
                        logger.warning(msg)
    if opts.overwrite:
        msg = f"Overwriting input {opts.input} with temporary output {opts.output} GFF file"
        logger.info(msg)
        shutil.copyfile(opts.output, opts.input)
        os.remove(opts.output)
                        
if __name__ == "__main__":
    main()
    
    

    

