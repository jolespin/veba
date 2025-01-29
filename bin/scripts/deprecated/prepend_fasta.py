#!/usr/bin/env python
import sys,os, argparse, shutil
from tqdm import tqdm
import pyfastx
from pyexeggutor import (
    build_logger,
    open_file_writer, 
    fasta_writer,
)

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "v2024.11.7"

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
    parser_io.add_argument("-i", "--input", type=str,  required=True, help = f"path/to/input.fasta[.gz]")
    parser_io.add_argument("-n", "--prefix", type=str,  required=True, help = f"Genome name used for prefix followed by some separator [Recommended: Genome identifier (e.g., Monili1) followed by separator (e.g., Monili2__)]")
    parser_io.add_argument("-o", "--output", type=str, help = f"path/to/output.fasta[.gz].  [Default: original_path.MODIFIED.fasta[.gz]]")
    parser_io.add_argument("-w", "--wrap", type=int, default=1000, help = f"Line width for wrapping fasta lines.  If 0, then no wrapping is used. [Default: 1000]")
    parser_io.add_argument("--overwrite", action="store_true", help = f"Overwrite the fasta file")
    
    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    # Logger
    logger = build_logger(f"Prepend label to fasta identifiers [{opts.script_filename}]: {opts.prefix}")
    # reset_logger(logger)

    # Commands
    logger.info(f"Command: {sys.argv}")
    
    # Output
    if not opts.output:
        fields = opts.input.split(".")
        if opts.input.endswith(".gz"):
            ext = fields[-2]
            fileprefix = ".".join(fields[:-2])
            opts.output = fileprefix + ".MODIFIED" + f".{ext}.gz"
        else:
            ext = fields[-1]
            fileprefix = ".".join(fields[:-1])
            opts.output = fileprefix + ".MODIFIED" + f".{ext}"
        logger.info(f"--output was not provided.  Writing (temporary) output: {opts.output}")
    
    msg = f"Reading input {opts.input} and writing output {opts.output} fasta files"
    logger.info(msg)

    with (
        open_file_writer(opts.output) as f_out,
        ):
        for id, seq in tqdm(pyfastx.Fasta(opts.input, build_index=False), desc=msg):
            header = f"{opts.prefix}{id}"
            fasta_writer(header, seq, file=f_out, wrap=opts.wrap)

    if opts.overwrite:
        msg = f"Overwriting input {opts.input} with temporary output {opts.output} fasta file"
        logger.info(msg)
        shutil.copyfile(opts.output, opts.input)
        os.remove(opts.output)                        
if __name__ == "__main__":
    main()
    
    

    

