#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
from collections import OrderedDict, defaultdict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.09.29"

# .............................................................................
# Notes
# .............................................................................
# * Make batch version that takes in a manifest file
# .............................................................................
# Primordial
# .............................................................................

# Assembly
def preprocess( input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]

    if not os.path.exists(os.path.join(directories["checkpoints"], "0__preprocessing")):

        # Proteins
        df = pd.read_csv(opts.proteins, sep="\t", index_col=0, header=None)
        assert df.shape[1] in {0,1}, "Must either be 1 (for list) or 2 (for table) columns"

        proteins = dict()
        if df.shape[1] == 1:
            proteins = df.iloc[:,0].dropna().to_dict()
        else:
            for path in df.index:
                id_mag = path.split("/")[-1][:-1*(len(opts.proteins_extension) + 1)]
                assert id_mag not in proteins, "--proteins has non-unique MAG identifiers"
                proteins[id_mag] = path 
        proteins = pd.Series(proteins)

        # Write files
        proteins.to_frame().to_csv(os.path.join(directories["project"], "proteins.tsv"), sep="\t", header=None)


        os.makedirs(os.path.join(output_directory, "proteins"), exist_ok=True)

        # Proteins
        for id_mag, path in proteins.items():
                src = os.path.realpath(path)
                dst = os.path.join(os.path.join(output_directory, "proteins", "{}.faa".format(id_mag)))
                os.symlink(src,dst)

                # cmd += [
                # "ln -sf {} {}".format(src,dst),
                # "&&",
                # ]

        # cmd += [
        #     "ls {} > {}".format(
        #         os.path.join(os.path.join(output_directory, "proteins","*.faa" )),
        #         os.path.join(os.path.join(directories["project"], "proteins.list" )),
        #         ),
        # ]
        with open(os.path.join(os.path.join(directories["project"], "proteins.list" )), "w") as f:
            for fp in glob.glob(os.path.join(os.path.join(output_directory, "proteins","*.faa" ))):
                print(fp, file=f)
    return []

# # HMMER
# def get_hmmsearch_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

#     os.environ["TMPDIR"] = directories["tmp"]
#     # Command
#     cmd = [
#         "(",
#         os.environ["hmmsearch"],
#         "--tblout {}".format(os.path.join(output_directory, "output.tsv")),
#         # "-E {}".format(opts.hmmsearch_evalue),
#         "--cpu {}".format(opts.n_jobs),
#         "--seed {}".format(opts.random_state + 1),
#     ]
#     if bool(opts.hmmsearch_threshold):
#         cmd += [ "--{}".format(opts.hmmsearch_threshold)]

#     cmd += [
#         opts.database_hmm,
#         opts.proteins,
#         ">",
#         "/dev/null",
#         ")",
#     ]

#     cmd += [ 
#         "&&",
#         "pigz",
#         "-f",
#         "-p {}".format(opts.n_jobs),
#         os.path.join(output_directory, "output.tsv"),
#     ]    
#     return cmd

# HMMSearch
def get_hmmsearch_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]
    # Command

#     # Command
#     cmd = [
#         "(",
#         os.environ["hmmsearch"],
#         "--tblout {}".format(os.path.join(output_directory, "output.tsv")),
#         # "-E {}".format(opts.hmmsearch_evalue),
#         "--cpu {}".format(opts.n_jobs),
#         "--seed {}".format(opts.random_state + 1),
#     ]
#     if bool(opts.hmmsearch_threshold):
#         cmd += [ "--{}".format(opts.hmmsearch_threshold)]

#     cmd += [
#         opts.database_hmm,
#         opts.proteins,
#         ">",
#         "/dev/null",
#         ")",
#     ]
    cmd = [
"""

for ID in $(cut -f1 %s);
    do echo $ID; 
    FAA=%s/${ID}.faa
    OUT_DIR=%s/${ID}
    mkdir -p ${OUT_DIR}

    # Run HMMSearch
    %s --tblout ${OUT_DIR}/output.tsv --cpu %d --seed %d %s %s ${FAA} > /dev/null

    # Parse HMMSearch
    #%s -i ${OUT_DIR}/output.tsv -a ${FAA} -o ${OUT_DIR} -n ${ID}

    done
"""%(
    input_filepaths[0],
    os.path.join(directories["preprocessing"], "proteins"),
    output_directory,
    os.environ["hmmsearch"],
    opts.n_jobs,
    opts.random_state + 1,
    {True:"--{}".format(opts.hmmsearch_threshold), False:"-E {}".format(opts.hmmsearch_evalue)}[bool(opts.hmmsearch_threshold)],
    opts.database_hmm,
    os.environ["partition_hmmsearch.py"],
    )
    ]
    return cmd



# Symlink
def get_symlink_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    
    # Command
    cmd = ["("]
    for filepath in input_filepaths:
        # cmd.append("ln -f -s {} {}".format(os.path.realpath(filepath), os.path.realpath(output_directory)))
        cmd.append("ln -f -s {} {}".format(os.path.realpath(filepath), output_directory))
        cmd.append("&&")
    cmd[-1] = ")"
    return cmd

# ============
# Run Pipeline
# ============
# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = set([ 
        "partition_hmmsearch.py",
    ])

    
    required_executables={
                # 1
                "hmmsearch",
                "muscle",
                "iqtree",

     } | accessory_scripts

    if opts.path_config == "CONDA_PREFIX":
        executables = dict()
        for name in required_executables:
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
        executables[name] = "python " + os.path.join(opts.script_directory, "scripts", name)
    print(format_header( "Adding executables to path from the following source: {}".format(opts.path_config), "-"), file=sys.stdout)
    for name, executable in executables.items():
        if name in required_executables:
            print(name, executable, sep = " --> ", file=sys.stdout)
            os.environ[name] = executable.strip()
    print("", file=sys.stdout)

# Pipeline
def create_pipeline(opts, directories, f_cmds):

    # .................................................................
    # Primordial
    # .................................................................
    # Commands file
    pipeline = ExecutablePipeline(name=__program__,  f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])

    
    # ==========
    # Preprocessing
    # ==========
    
    program = "preprocessing"
    # Add to directories
    output_directory = directories["preprocessing"] 

    # Info
    step = 0
    description = "Symlink proteins"

    # i/o
    input_filepaths = [opts.proteins]
    output_filepaths = [
        os.path.join(directories["project"], "proteins.list"),
        os.path.join(directories["project"], "proteins.tsv"),
    ]

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = preprocess(**params)

    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
    )

  
    # ==========
    # HMMSearch
    # ==========
    step = 1

    program = "hmmsearch"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "HMMSearch for each MAG"


    # i/o
    input_filepaths = [
        os.path.join(directories["project"], "proteins.tsv"),
        ]
    output_filenames = ["*/*.tsv"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_hmmsearch_cmd(**params)

    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=False,
    )

    #     # =============
    #     # Symlink
    #     # =============
    #     program = "symlink"
    #     # Add to directories
    #     output_directory = directories["output"]

    #     # Info
    #     step = 3
    #     description = "Symlinking relevant output files"

    #     # i/o
    #     input_filepaths = [
    #         os.path.join(directories[("intermediate", "1__fastani")], "clusters.tsv"),
    #         os.path.join(directories[("intermediate", "2__orthofinder")], "proteins_to_orthogroups.tsv"),
    #     ]

    #     output_filenames =  map(lambda fp: fp.split("/")[-1], input_filepaths)
    #     output_filepaths = list(map(lambda fn:os.path.join(directories["output"], fn), output_filenames))

    #     params = {
    #     "input_filepaths":input_filepaths,
    #     "output_filepaths":output_filepaths,
    #     "output_directory":output_directory,
    #     "opts":opts,
    #     "directories":directories,
    #     }

    #     cmd = get_symlink_cmd(**params)
    #     pipeline.add_step(
    #             id=program,
    #             description = description,
    #             step=step,
    #             cmd=cmd,
    #             input_filepaths = input_filepaths,
    #             output_filepaths = output_filepaths,
    #             validate_inputs=True,
    #             validate_outputs=False,
    #     )

    return pipeline

# Configure parameters
def configure_parameters(opts, directories):
    if opts.hmmsearch_threshold.lower() == "e":
        opts.hmmsearch_threshold = None
    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -d <database_hmms> -a <proteins> -o <output_directory>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-d", "--database_hmm", type=str,  help=f"path/to/HMM database of markers")
    parser_io.add_argument("-a","--proteins", type=str, help = "Tab-seperated value table of [id_mag]<tab>[path/to/protein.fasta]")
    parser_io.add_argument("-o","--output_directory", type=str, default="veba_output/phylogeny", help = "path/to/project_directory [Default: veba_output/phylogeny]")
    parser_io.add_argument("--proteins_extension", type=str, default="faa", help = "Fasta file extension for proteins if a list is provided [Default: faa]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # HMMER
    parser_hmmer = parser.add_argument_group('HMMSearch arguments')
    parser_hmmer.add_argument("--hmmsearch_threshold", type=str, default="cut_ga", help="HMMSearch | Threshold {cut_ga, cut_nc, gut_tc, e} [Default:  cut_ga]")
    parser_hmmer.add_argument("--hmmsearch_evalue", type=float, default=10.0, help="Diamond | E-Value [Default: 10.0]")
    parser_hmmer.add_argument("--hmmsearch_options", type=str, default="", help="Diamond | More options (e.g. --arg 1 ) [Default: '']")



    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Directories
    directories = dict()
    directories["project"] = create_directory(opts.output_directory)
    directories["preprocessing"] = create_directory(os.path.join(directories["project"], "preprocessing"))
    directories["output"] = create_directory(os.path.join(directories["project"], "output"))
    directories["log"] = create_directory(os.path.join(directories["project"], "log"))
    directories["tmp"] = create_directory(os.path.join(directories["project"], "tmp"))
    directories["checkpoints"] = create_directory(os.path.join(directories["project"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["project"], "intermediate"))

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

    # Symlink genomes

    # Run pipeline
    with open(os.path.join(directories["project"], "commands.sh"), "w") as f_cmds:
        pipeline = create_pipeline(
                     opts=opts,
                     directories=directories,
                     f_cmds=f_cmds,
        )
        pipeline.compile()
        pipeline.execute(restart_from_checkpoint=opts.restart_from_checkpoint)

if __name__ == "__main__":
    main()
