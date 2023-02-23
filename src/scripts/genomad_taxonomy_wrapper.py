#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, shutil
from collections import OrderedDict, defaultdict

import pandas as pd
import numpy as np

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.2.17"


# MetaEuk
def get_genomad_taxonomy_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [
        "cat",
        opts.fasta,
        "|",
        os.environ["seqkit"],
        "replace",
        '-p "\s.+"',
        ">",
        os.path.join(directories["tmp"], "input.fasta"),

        "&&",

        os.environ["genomad"],
        "annotate",
        "--cleanup",
        "--threads {}".format(opts.n_jobs),
        "--splits {}".format(opts.splits),
        "--evalue {}".format(opts.mmseqs2_evalue),
        "--sensitivity {}".format(opts.sensitivity),
    ]

    if opts.use_minimal_database_for_taxonomy:
        cmd += [
            "--use-minimal-db",
        ]

    cmd += [ 
        os.path.join(directories["tmp"], "input.fasta"), # Input
        os.path.join(output_directory, "intermediate"), # Output
        os.path.join(opts.veba_database, "Classify", "geNomad"),

            "&&",

        "cp -rf", 
        os.path.join(output_directory, "intermediate", "input_annotate", "input_taxonomy.tsv"), 
        os.path.join(output_directory, "viral_taxonomy.tsv"), 

            "&&",

        "rm -rf",
        os.path.join(directories["tmp"], "input.fasta"),

        ]


    return cmd


# ============
# Run Pipeline
# ============



def create_pipeline(opts, directories, f_cmds):

    # .................................................................
    # Primordial
    # .................................................................
    # Commands file
    pipeline = ExecutablePipeline(name=__program__, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])


    # ==========
    # Binning
    # ==========
    step = 1

    program = "genomad_taxonomy"

    program_label = "{}__{}".format(step, program)

    # Add to directories


    # Info
    description = "Classifying viral contigs with geNomad"
    
    # i/o
    input_filepaths = [
        opts.fasta, 
        ]

    output_directory = directories["output"]
    output_filepaths = [ 
        os.path.join(output_directory, "viral_taxonomy.tsv"),
    ]

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_genomad_taxonomy_cmd(**params)


    # Add step to pipeline
    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
                log_prefix=program_label,
    )


    return pipeline



# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = set([ 

    ])
    

    required_executables={
        "seqkit",
        "genomad",
     } 
  

    required_executables |= accessory_scripts

    if opts.path_config == "CONDA_PREFIX":
        executables = dict()
        for name in sorted(required_executables):
            if name not in accessory_scripts:
                executables[name] = os.path.join(os.environ["CONDA_PREFIX"], "bin", name)
    else:
        if opts.path_config is None:
            opts.path_config = os.path.join(opts.script_directory, "veba_config.tsv")
        opts.path_config = format_path(opts.path_config)
        assert os.path.exists(opts.path_config), "config file does not exist.  Have you created one in the following directory?\n{}\nIf not, either create one, check this filepath:{}, or give the path to a proper config file using --path_config".format(opts.script_directory, opts.path_config)
        assert os.stat(opts.path_config).st_size > 1, "config file seems to be empty.  Please add 'name' and 'executable' columns for the following program names: {}".format(required_executables)
        df_config = pd.read_csv(opts.path_config, sep="\t")
        assert {"name", "executable"} <= set(df_config.columns), "config must have `name` and `executable` columns.  Please adjust file: {}".format(opts.path_config)
        df_config = df_config.loc[:,["name", "executable"]].dropna(how="any", axis=0).applymap(str)
        # Get executable paths
        executables = OrderedDict(zip(df_config["name"], df_config["executable"]))
        assert required_executables <= set(list(executables.keys())), "config must have the required executables for this run.  Please adjust file: {}\nIn particular, add info for the following: {}".format(opts.path_config, required_executables - set(list(executables.keys())))

    # Display

    for name in sorted(accessory_scripts):
        if name.endswith(".py"):
            executables[name] = "python " + os.path.join(opts.script_directory, name)
        else: 
            executables[name] = os.path.join(opts.script_directory, name)


    print(format_header( "Adding executables to path from the following source: {}".format(opts.path_config), "-"), file=sys.stdout)
    for name, executable in executables.items():
        if name in required_executables:
            print(name, executable, sep = " --> ", file=sys.stdout)
            os.environ[name] = executable.strip()
    print("", file=sys.stdout)


# Configure parameters
def configure_parameters(opts, directories):
    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -f <scaffolds.fasta> -i <scaffolds_to_bins.tsv>  -o <output_directory> --veba_database <path>".format(__program__)

    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-f","--fasta", type=str, required=True, help = "path/to/scaffolds.fasta")
    parser_io.add_argument("-o","--output_directory", type=str, default="genomad_taxonomy_output", help = "path/to/project_directory [Default: genomad_taxonomy_output]")
    parser_io.add_argument("-i","--scaffolds_to_bins", type=str, required=False,  help = "path/to/scaffolds_to_bins.tsv, [Optional] Format: [id_scaffold]<tab>[id_bin], No header")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    # parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Databases
    parser_databases = parser.add_argument_group('Database arguments')
    parser_databases.add_argument("--veba_database", type=str,  help=f"VEBA database location.  [Default: $VEBA_DATABASE environment variable]")

    # geNomad
    parser_genomad = parser.add_argument_group('geNomad arguments')
    parser_genomad.add_argument("--sensitivity", type=float, default=4.0, help = "MMseqs2 marker search sensitivity. Higher values will annotate more proteins, but the search will be slower and consume more memory. [Default: 4.0; x ≥ 0.0]")
    parser_genomad.add_argument("--splits", type=int, default=0, help = "Split the data for the MMseqs2 search. Higher values will reduce memory usage, but will make the search slower. If the MMseqs2 search is failing, try to increase the number of splits. Also used for VirFinder. [Default: 0; x ≥ 0]")
    parser_genomad.add_argument("--mmseqs2_evalue", type=float, default=1e-3, help = "Maximum accepted E-value in the MMseqs2 search. Used by genomad annotate when VirFinder is used as binning algorithm [Default: 1e-3]")
    parser_genomad.add_argument("--use_minimal_database_for_taxonomy", action="store_true", help="Use a smaller marker database to annotate proteins. This will make execution faster but sensitivity will be reduced.")

    # Unclassified
    parser_unclassified = parser.add_argument_group('Unclassified arguments')
    parser_unclassified.add_argument("-R", "--remove_unclassified", action="store_true", help = "Remove unclassified contigs")
    parser_unclassified.add_argument("-l", "--unclassified_lineage_label", type=str, default="Unclassified", help = "Unclassified lineage label [Default: Unclassified]")
    parser_unclassified.add_argument("-t", "--unclassified_taxid", type=int, default=-1, help = "Unclassified taxon id [Default: -1]")
    parser_unclassified.add_argument("-s", "--unclassified_score", type=float, default=1.0, help = "Unclassified agreement score [Default: 1.0]")

    # Options
    opts = parser.parse_args()

    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Threads
    if opts.n_jobs == -1:
        from multiprocessing import cpu_count 
        opts.n_jobs = cpu_count()
    assert opts.n_jobs >= 1, "--n_jobs must be ≥ 1.  To select all available threads, use -1."

    # Database
    if opts.veba_database is None:
        assert "VEBA_DATABASE" in os.environ, "Please set the following environment variable 'export VEBA_DATABASE=/path/to/veba_database' or provide path to --veba_database"
        opts.veba_database = os.environ["VEBA_DATABASE"]

    # Directories
    directories = dict()
    directories["output"] = create_directory(opts.output_directory)
    directories["intermediate"] = create_directory(os.path.join(directories["output"], "intermediate"))
    directories["log"] = create_directory(os.path.join(directories["output"], "log"))
    directories["tmp"] = create_directory(os.path.join(directories["output"],"intermediate", "tmp"))
    directories["checkpoints"] = create_directory(os.path.join(directories["output"],"intermediate", "checkpoints"))
    os.environ["TMPDIR"] = directories["tmp"]



    # Info
    print(format_header(__program__, "="), file=sys.stdout)
    print(format_header("Configuration:", "-"), file=sys.stdout)
    print("Python version:", sys.version.replace("\n"," "), file=sys.stdout)
    print("Python path:", sys.executable, file=sys.stdout) #sys.path[2]
    print("Script version:", __version__, file=sys.stdout)
    print("Moment:", get_timestamp(), file=sys.stdout)
    print("Directory:", os.getcwd(), file=sys.stdout)
    print("Commands:", list(filter(bool,sys.argv)),  sep="\n", file=sys.stdout)
    configure_parameters(opts, directories)
    sys.stdout.flush()


    # Run pipeline
    with open(os.path.join(directories["output"], "commands.sh"), "w") as f_cmds:
        pipeline = create_pipeline(
                    opts=opts,
                    directories=directories,
                    f_cmds=f_cmds,
        )
        pipeline.compile()
        pipeline.execute(restart_from_checkpoint=opts.restart_from_checkpoint)

    # shutil.rmtree(directories["intermediate"])

    update_table = False
    df_genomad_taxonomy = pd.read_csv(os.path.join(directories["output"], "viral_taxonomy.tsv"), sep="\t", index_col=0)

    if not opts.remove_unclassified:

        identifiers = list()
        with get_file_object(opts.fasta, mode="read") as f:
            for line in pv(f.readlines(), description="Reading {}".format(opts.fasta)):
                header = line.strip()
                if header.startswith(">"):
                    id_contig = header.split(" ")[0][1:]
                    identifiers.append(id_contig)

        print(" * Filling in missing taxonomy information", file=sys.stdout)
        df_genomad_taxonomy = df_genomad_taxonomy.reindex(identifiers)
        mask = df_genomad_taxonomy["lineage"].isnull()
        df_genomad_taxonomy.loc[mask,"lineage"] = opts.unclassified_lineage_label
        df_genomad_taxonomy.loc[mask,"taxid"] = opts.unclassified_taxid
        df_genomad_taxonomy["taxid"] = df_genomad_taxonomy["taxid"].astype(int)
        df_genomad_taxonomy.loc[mask,"agreement"] = opts.unclassified_score
        df_genomad_taxonomy.loc[mask,"n_genes_with_taxonomy"] = 0
        df_genomad_taxonomy["n_genes_with_taxonomy"] = df_genomad_taxonomy["n_genes_with_taxonomy"].astype(int)
        df_genomad_taxonomy.index.name = "id_genome"
        update_table = True


    if opts.scaffolds_to_bins:
        print(" * Adding MAG identifiers", file=sys.stdout)
        scaffold_to_bin = pd.read_csv(opts.scaffolds_to_bins, sep="\t", index_col=0, header=None).iloc[:,0]
        df_genomad_taxonomy.insert(loc=0, column="id_contig", value=df_genomad_taxonomy.index)
        df_genomad_taxonomy.index = df_genomad_taxonomy.index.map(lambda x: scaffold_to_bin[x])
        df_genomad_taxonomy.index.name = "id_genome"
        update_table = True

    if update_table:
        df_genomad_taxonomy.to_csv(os.path.join(directories["output"], "viral_taxonomy.tsv"), sep="\t")



    # if (not opts.remove_unclassified) or (opts.scaffolds_to_bins):
    #     validate_file_existence([os.path.join(directories["output"], "viral_taxonomy.tsv")], prologue="Validating the following updated files:")

   


if __name__ == "__main__":
    main(sys.argv[1:])
