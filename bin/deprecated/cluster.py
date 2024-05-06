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
__version__ = "2023.1.20"


# Assembly
def preprocess( input_filepaths, output_filepaths, output_directory, directories, opts):
    if not os.path.exists(os.path.join(directories["checkpoints"], "0__preprocessing")):

            # MAGs
        df = pd.read_csv(opts.mags, sep="\t", index_col=0, header=None)
        assert df.shape[1] in {0,1}, "Must either be 1 (for list) or 2 (for table) columns"

        mags = dict()
        if df.shape[1] == 1:
            mags = df.iloc[:,0].dropna().to_dict()
        else:
            for path in df.index:
                id_mag = path.split("/")[-1][:-1*(len(opts.mags_extension) + 1)]
                assert id_mag not in mags, "--mags has non-unique MAG identifiers"
                mags[id_mag] = path 
        mags = pd.Series(mags)
        
        # Proteins
        df = pd.read_csv(opts.proteins, sep="\t", index_col=0, header=None)
        assert df.shape[1] in {0,1}, "Must either be 1 (for list) or 2 (for table) columns"

        proteins = dict()
        if df.shape[1] == 1:
            proteins = df.iloc[:,0].dropna().to_dict()
        else:
            for path in df.index:
                id_mag = path.split("/")[-1][:-1*(len(opts.proteins_extension) + 1)]
                if id_mag in mags:
                    assert id_mag not in proteins, "--proteins has non-unique MAG identifiers"

                    proteins[id_mag] = path 
        proteins = pd.Series(proteins)


        assert  set(mags.index) <= set(proteins.index), "--mags must be a subset of --proteins:\nThe following --mags do not have proteins:\n{}".format("\n * ".join(set(mags.index) - set(proteins.index)))
        proteins = proteins.loc[mags.index]

        # Write files
        mags.to_frame().to_csv(os.path.join(directories["project"], "genomes.tsv"), sep="\t", header=None)
        proteins.to_frame().to_csv(os.path.join(directories["project"], "proteins.tsv"), sep="\t", header=None)

        # # Command
        # cmd = [
        #     "mkdir -p {}".format(os.path.join(output_directory, "genomes")),
        #     "&&",
        #     "mkdir -p {}".format(os.path.join(output_directory, "proteins")),
        #     "&&",
        # ]

        # # MAGs
        # for id_mag, path in mags.items():
        #     src = os.path.realpath(path)
        #     dst = os.path.join(os.path.join(output_directory, "genomes", "{}.fa".format(id_mag)))
        #     cmd += [
        #         "ln -sf {} {}".format(src,dst),
        #         "&&",
        #     ]
        # cmd += [
        #     "ls {} > {}".format(
        #         os.path.join(os.path.join(output_directory, "genomes","*.fa" )),
        #         os.path.join(os.path.join(directories["project"], "genomes.list" )),
        #         ),
        #     "&&",
        # ]

        # # Proteins
        # for id_mag, path in proteins.items():
        #         src = os.path.realpath(path)
        #         dst = os.path.join(os.path.join(output_directory, "proteins", "{}.faa".format(id_mag)))
        #         cmd += [
        #         "ln -sf {} {}".format(src,dst),
        #         "&&",
        #         ]

        # cmd += [
        #     "ls {} > {}".format(
        #         os.path.join(os.path.join(output_directory, "proteins","*.faa" )),
        #         os.path.join(os.path.join(directories["project"], "proteins.list" )),
        #         ),
        # ]
        # Command
        # cmd = [
        #     "mkdir -p {}".format(os.path.join(output_directory, "genomes")),
        #     "&&",
        #     "mkdir -p {}".format(os.path.join(output_directory, "proteins")),
        #     "&&",
        # ]
        os.makedirs(os.path.join(output_directory, "genomes"), exist_ok=True)
        os.makedirs(os.path.join(output_directory, "proteins"), exist_ok=True)

    # MAGs
        for id_mag, path in mags.items():
            src = os.path.realpath(path)
            dst = os.path.join(os.path.join(output_directory, "genomes", "{}.fa".format(id_mag)))
            os.symlink(src,dst)
            # cmd += [
            #     "ln -sf {} {}".format(src,dst),
            #     "&&",
            # ]
        # cmd += [
        #     "ls {} > {}".format(
        #         os.path.join(os.path.join(output_directory, "genomes","*.fa" )),
        #         os.path.join(os.path.join(directories["project"], "genomes.list" )),
        #         ),
        #     "&&",
        # ]
        with open(os.path.join(os.path.join(directories["project"], "genomes.list" )), "w") as f:
            for fp in glob.glob(os.path.join(os.path.join(output_directory, "genomes","*.fa" ))):
                print(fp, file=f)




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


# FastANI
def get_fastani_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # fastANI -t 4 --rl ./mags.list --ql ./mags.list -o fastani_output/output.txt --matrix fastani_output/output.matrix

    
    # Command
    cmd = [
        "(",
        os.environ["fastANI"],
        "-t {}".format(opts.n_jobs),
        "--rl {}".format(input_filepaths[0]),
        "--ql {}".format(input_filepaths[0]),
        "-o {}".format(output_filepaths[0]),
        opts.fastani_options,
        ")",
        "&&",
        "(",
        os.environ["fastani_to_clusters.py"],
        "-i {}".format(output_filepaths[0]),
        "-x fa",
        "-t {}".format(opts.ani_threshold),
        "-g {}".format(input_filepaths[0]),
        {True:"-p {}".format(opts.cluster_prefix), False:""}[bool(opts.cluster_prefix)],
        {True:"-s {}".format(opts.cluster_suffix), False:""}[bool(opts.cluster_suffix)],
        "--export_pickle {}".format(os.path.join(output_directory, "genomes")),
        "--export_edgelist {}".format(os.path.join(output_directory, "genomes")),
        "-o {}".format(output_filepaths[1]),
        "--no_header",
        ")",
        "&&",
        "(",
        os.environ["partition_clusters.py"],
        "-c {}".format(output_filepaths[1]),
        "-a {}".format(input_filepaths[1]),
        "-o {}".format(os.path.join(output_directory, "clusters")),
        {True:"--copy",False:""}[bool(opts.copy_proteins)],
        "--clone_singletons",
        "--clone_label {}".format(opts.clone_label),
        ")",
        "&&",
        os.environ["compile_scaffold_identifiers.py"],
        "-i {}".format(opts.scaffolds_to_bins),
        "-c {}".format(output_filepaths[1]),
        "-o {}".format( os.path.join(output_directory, "scaffolds_to_clusters.tsv")),
               
    ]
    return cmd

# # FastANI
# def get_orthofinder_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
#     # orthofinder -t ${N_JOBS} -a ${N_JOBS} -f ${DIR} -M msa -o 
#     # REMEMBER TO GET DUMMIES FOR SINGLETON CLUSTERS
#     os.environ["TMPDIR"] = directories["tmp"]
#     # Command
#     cmd = [
# """

# for ID in $(cut -f1 %s);
#     do echo $ID; 
#     FAA_DIR=%s/${ID}
#     OUT_DIR=%s/${ID}
#     rm -rf ${OUT_DIR}
#     ORTHOGROUP_FILE=${OUT_DIR}/Results_orthofinder/Orthogroups/Orthogroups.tsv

#     # Run Orthofinder
#     %s -t %d -a %d -f ${FAA_DIR} -o ${OUT_DIR} -n orthofinder -M msa %s

#     # Create manifest file
#     echo "${ID}\t${ORTHOGROUP_FILE}" >> %s

#     done

# # Remove big intermediate files
# rm -rf %s
# rm -rf %s
# rm -rf %s
# rm -rf %s

# """%(
#     os.path.join(input_filepaths[0],"clusters_to_proteins.tsv"),
#     input_filepaths[0],
#     output_directory,
#     os.environ["orthofinder"],
#     opts.n_jobs,
#     opts.n_jobs,
#     opts.orthofinder_options,
#     output_filepaths[0],
#     os.path.join(output_directory, "*", "Results_orthofinder", "WorkingDirectory"),
#     os.path.join(output_directory, "*", "Results_orthofinder", "Orthogroup_Sequences"),
#     os.path.join(output_directory, "*", "Results_orthofinder", "Gene_Trees"),
#     os.path.join(output_directory, "*", "Results_orthofinder", "Single_Copy_Orthologue_Sequences"),

#     ),
#     "(",
#     os.environ["partition_orthogroups.py"],
#     "-i {}".format(output_filepaths[0]),
#     "-c {}".format(input_filepaths[1]),
#     "-a {}".format(input_filepaths[2]),
#     "-x {}".format(opts.proteins_extension),
#     "--clone_label {}".format(opts.clone_label),
#     "--sep _",
#     "-o {}".format(output_filepaths[1]),
#     ")",

    
#     ]
#     return cmd

# MUSCLE
def get_orthofinder_gnuparallel_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command

    cmd = [
"""
# Orthofinder
cut -f1 %s | %s -j %d 'echo "[OrthoFinder] {}" && rm -rf %s/{} && %s -t 1 -a 1 -f %s/{} -o %s/{} -n orthofinder -M msa %s && echo "{}\t%s/{}/Results_orthofinder/Orthogroups/Orthogroups.tsv" >> %s'

# Remove big intermediate files
rm -rf %s
rm -rf %s
rm -rf %s
rm -rf %s


"""%( 
    # Args
    os.path.join(input_filepaths[0],"clusters_to_proteins.tsv"),
    os.environ["parallel"],
    opts.n_jobs,
    output_directory,
    os.environ["orthofinder"],
    input_filepaths[0],
    output_directory,
    opts.orthofinder_options,
    output_directory,
    output_filepaths[0],

    os.path.join(output_directory, "*", "Results_orthofinder", "WorkingDirectory"),
    os.path.join(output_directory, "*", "Results_orthofinder", "Orthogroup_Sequences"),
    os.path.join(output_directory, "*", "Results_orthofinder", "Gene_Trees"),
    os.path.join(output_directory, "*", "Results_orthofinder", "Single_Copy_Orthologue_Sequences"),

    ),
    "(",
    os.environ["partition_orthogroups.py"],
    "-i {}".format(output_filepaths[0]),
    "-c {}".format(input_filepaths[1]),
    "-a {}".format(input_filepaths[2]),
    "-x {}".format(opts.proteins_extension),
    "--clone_label {}".format(opts.clone_label),
    "--sep _",
    "-o {}".format(output_filepaths[1]),
    ")",
    "&&",
    "(",
    "cat",
     output_filepaths[1],
     "|", 
     "cut -f1,5",
     "|", 
     "tail -n +2",
      ">",
      output_filepaths[2],
      ")",
       
    ]

    return cmd

# MUSCLE
def get_orthofinder_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command

    cmd = [
"""
n=1
for ID in $(cut -f1 %s); do
    # Remove directory
    rm -rf %s/${ID}

    START_TIME=${SECONDS}

    # Run Orthofinder
    %s -t %d -a %d -f %s/${ID} -o %s/${ID} -n orthofinder -M msa %s

    # Placehold
    echo "${ID}\t%s/${ID}/Results_orthofinder/Orthogroups/Orthogroups.tsv" >> %s

    # Remove big intermediate files
    rm -rf %s
    rm -rf %s
    rm -rf %s
    rm -rf %s

    END_TIME=${SECONDS}
    RUN_TIME=$((END_TIME-START_TIME))
    echo "*** n=${n} // ${ID} // Duration: ${RUN_TIME} seconds ***"

    n=$(($n+1))
    done

"""%( 
    # Args
    os.path.join(input_filepaths[0],"clusters_to_proteins.tsv"),
    output_directory,

    os.environ["orthofinder"],
    opts.n_jobs,
    opts.n_jobs,
    input_filepaths[0],
    output_directory,
    opts.orthofinder_options,
    output_directory,
    output_filepaths[0],

    os.path.join(output_directory, "${ID}", "Results_orthofinder", "WorkingDirectory"),
    os.path.join(output_directory, "${ID}", "Results_orthofinder", "Orthogroup_Sequences"),
    os.path.join(output_directory, "${ID}", "Results_orthofinder", "Gene_Trees"),
    os.path.join(output_directory, "${ID}", "Results_orthofinder", "Single_Copy_Orthologue_Sequences"),


    ),
    "(",
    os.environ["partition_orthogroups.py"],
    "-i {}".format(output_filepaths[0]),
    "-c {}".format(input_filepaths[1]),
    "-a {}".format(input_filepaths[2]),
    "-x {}".format(opts.proteins_extension),
    "--clone_label {}".format(opts.clone_label),
    "--sep _",
    "-o {}".format(output_filepaths[1]),
    ")",
    "&&",
    "(",
    "cat",
     output_filepaths[1],
     "|", 
     "cut -f1,5",
     "|", 
     "tail -n +2",
      ">",
      output_filepaths[2],
      ")",
       
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

    accessory_scripts = {
        "fastani_to_clusters.py", 
        "partition_clusters.py",
        "partition_orthogroups.py", 
        "compile_scaffold_identifiers.py",
        "edgelist_to_clusters.py",
        }

    required_executables={
                # 1
                "fastANI",
                "orthofinder",
                "mmseqs",
                "parallel",


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
    description = "Symlink genomes and proteins"

    # i/o
    input_filepaths = [opts.mags, opts.proteins]
    output_filepaths = [
        os.path.join(directories["project"], "genomes.list"),
        os.path.join(directories["project"], "proteins.list"),
        os.path.join(directories["project"], "genomes.tsv"),
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
                errors_ok=True,
    )

    # ==========
    # FastANI
    # ==========
    step = 1

    program = "fastani"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "FastANI of genomes"

    # if opts.clusters:
    #     # i/o
    #     input_filepaths = [opts.clusters, os.path.join(directories["project"], "proteins.tsv")]
    #     output_filenames = ["clusters.tsv"]
    #     output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    #     cmd = [
    #         "cp", 
    #         input_filepaths[0], 
    #         output_filepaths[0],
    #     ]

    # else:
    # i/o
    input_filepaths = [
        os.path.join(directories["project"], "genomes.list"),
        os.path.join(directories["project"], "proteins.tsv"),
        ]
    output_filenames = [
        "fastani_output.tsv",
        "clusters.tsv",
        "clusters", 
        "scaffolds_to_clusters.tsv",
        "genomes-ani_{}.edgelist.tsv".format(opts.ani_threshold),
        "genomes-ani_{}.graph.pkl".format(opts.ani_threshold),
        ]

    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_fastani_cmd(**params)

    pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
                errors_ok=True,
    )
    if not opts.cluster_only:
        # ==========
        # OrthoFinder
        # ==========
        step = 2

        program = "orthofinder"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

        # Info
        description = "Orthologue detection for FastANI clusters"


        # i/o
        input_filepaths = [
            os.path.join(directories[("intermediate", "1__fastani")], "clusters"), 
            os.path.join(directories[("intermediate", "1__fastani")], "clusters.tsv"),
            os.path.join(directories["project"], "proteins.tsv"),

            ]
        output_filenames = ["clusters_to_orthogroups.tsv", "identifier_mapping.orthogroups.tsv", "proteins_to_orthogroups.tsv"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }
        if not opts.one_task_per_cpu:
            cmd = get_orthofinder_cmd(**params)
        else:
            cmd = get_orthofinder_gnuparallel_cmd(**params)

        pipeline.add_step(
                    id=program,
                    description = description,
                    step=step,
                    cmd=cmd,
                    input_filepaths = input_filepaths,
                    output_filepaths = output_filepaths,
                    validate_inputs=True,
                    validate_outputs=True,
                    errors_ok=True,
        )

        # =============
        # Symlink
        # =============
        program = "symlink"
        # Add to directories
        output_directory = directories["output"]

        # Info
        step = 3
        description = "Symlinking relevant output files"

        # i/o
        input_filepaths = [
            os.path.join(directories[("intermediate", "1__fastani")], "clusters.tsv"),
            # os.path.join(directories[("intermediate", "1__fastani")], "fastani_output.tsv"),
            os.path.join(directories[("intermediate", "1__fastani")], "genomes-ani_{}.edgelist.tsv".format(opts.ani_threshold)),
            os.path.join(directories[("intermediate", "1__fastani")], "genomes-ani_{}.graph.pkl".format(opts.ani_threshold)),
            os.path.join(directories[("intermediate", "1__fastani")], "scaffolds_to_clusters.tsv"),

            os.path.join(directories[("intermediate", "2__orthofinder")], "identifier_mapping.orthogroups.tsv"),
            os.path.join(directories[("intermediate", "2__orthofinder")], "proteins_to_orthogroups.tsv"),

        ]

        output_filenames =  map(lambda fp: fp.split("/")[-1], input_filepaths)
        output_filepaths = list(map(lambda fn:os.path.join(directories["output"], fn), output_filenames))

        params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        }

        cmd = get_symlink_cmd(**params)
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

    return pipeline

# Configure parameters
def configure_parameters(opts, directories):

    assert "-" not in opts.clone_label, "--clone_label can't have a `-` in it because it gets misinterpreted by the command line arguments"

    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -m <mags> -a <proteins> -o <output_directory> -t 95".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i","--scaffolds_to_bins", type=str, required=True,  help = "path/to/scaffolds_to_bins.tsv, Format: [id_scaffold]<tab>[id_bin], No header")
    parser_io.add_argument("-m","--mags", type=str, help = "Tab-seperated value table of [id_mag]<tab>[path/to/genome.fasta]")
    parser_io.add_argument("-a","--proteins", type=str, help = "Tab-seperated value table of [id_mag]<tab>[path/to/protein.fasta]")
    # parser_io.add_argument("-c","--clusters", type=str, help = "Tab-seperated value table of [id_mag]<tab>[id_cluster].  Use this to skip fastANI step [NOT TESTED WILL PROBABLY FAIL]")
    parser_io.add_argument("-o","--output_directory", type=str, default="veba_output/cluster", help = "path/to/project_directory [Default: veba_output/cluster]")
    parser_io.add_argument("--mags_extension", type=str, default="fa", help = "Fasta file extension for --mags if a list is provided [Default: fa]")
    parser_io.add_argument("--proteins_extension", type=str, default="faa", help = "Fasta file extension for proteins if a list is provided [Default: faa]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--cluster_only", action="store_true", help = "Only run FastANI") #!
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # FastANI
    parser_fastani = parser.add_argument_group('FastANI arguments')
    parser_fastani.add_argument("-t", "--ani_threshold", type=float, default=95.0, help="FastANI | Species-level clustering threshold [Default: 95.0]")
    parser_fastani.add_argument("--fastani_options", type=str, default="", help="FastANI | More options (e.g. --arg 1 ) [Default: '']")
    parser_fastani.add_argument("--cluster_prefix", type=str, default="SLC", help="Cluster prefix [Default: 'SLC")
    parser_fastani.add_argument("--cluster_suffix", type=str, default="", help="Cluster suffix [Default: '")
    parser_fastani.add_argument("--copy_proteins", action="store_true", help="Copy instead of symlink")
    parser_fastani.add_argument("--clone_label", type=str, default="___clone", help="Singleton clone label [Default: '___clone")

    # OrthoFinder
    parser_orthofinder = parser.add_argument_group('OrthoFinder arguments')
    parser_orthofinder.add_argument("--one_task_per_cpu", action="store_true", help="Use GNU parallel to run GNU parallel with 1 task per CPU.  Useful if all clusters are roughly the same size but inefficient if cluster sizes vary.")
    parser_orthofinder.add_argument("--orthofinder_options", type=str, default="", help="OrthoFinder | More options (e.g. --arg 1 ) [Default: '']")

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
