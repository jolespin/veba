#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, site
from collections import OrderedDict, defaultdict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.11.16"


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
        
        # Gene models
        df = pd.read_csv(opts.gene_models, sep="\t", index_col=0, header=None)
        assert df.shape[1] in {0,1}, "Must either be 1 (for list) or 2 (for table) columns"

        gene_models = dict()
        if df.shape[1] == 1:
            gene_models = df.iloc[:,0].dropna().to_dict()
        else:
            for path in df.index:
                id_mag = path.split("/")[-1][:-1*(len(opts.gene_models_extension) + 1)]
                if id_mag in mags:
                    assert id_mag not in gene_models, "--gene_models has non-unique MAG identifiers"

                    gene_models[id_mag] = path 
        gene_models = pd.Series(gene_models)


        assert  set(mags.index) <= set(gene_models.index), "--mags must be a subset of --gene_models:\nThe following --mags do not have gene_models:\n{}".format("\n * ".join(set(mags.index) - set(gene_models.index)))
        gene_models = gene_models.loc[mags.index]

        # Write files
        pd.concat([mags.to_frame(), gene_models.to_frame()], axis=1).to_csv(os.path.join(directories["preprocessing"], "genomes_gene-models.tsv"), sep="\t", header=None)

    return []
    
        # os.makedirs(os.path.join(output_directory, "genomes"), exist_ok=True)
        # os.makedirs(os.path.join(output_directory, "gene_models"), exist_ok=True)
    # # MAGs
    #     for id_mag, path in mags.items():
    #         src = os.path.realpath(path)
    #         dst = os.path.join(os.path.join(output_directory, "genomes", "{}.fa".format(id_mag)))
    #         os.symlink(src,dst)

    #     with open(os.path.join(os.path.join(directories["project"], "genomes.list" )), "w") as f:
    #         for fp in glob.glob(os.path.join(os.path.join(output_directory, "genomes","*.fa" ))):
    #             print(fp, file=f)




    #     # gene_models
    #     for id_mag, path in gene_models.items():
    #             src = os.path.realpath(path)
    #             dst = os.path.join(os.path.join(output_directory, "gene_models", "{}.gff".format(id_mag)))
    #             os.symlink(src,dst)


    #     with open(os.path.join(os.path.join(directories["project"], "gene_models.list" )), "w") as f:
    #         for fp in glob.glob(os.path.join(os.path.join(output_directory, "gene_models","*.gff" ))):
    #             print(fp, file=f)
    # return []




# antiSMASH
def get_antismash_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command




    cmd = [
"""
mkdir -p %s

n=1
while IFS= read -r LINE
do read -r -a ARRAY <<< $LINE
    ID=${ARRAY[0]}
    GENOME=${ARRAY[1]}
    GENE_MODELS=${ARRAY[2]}

    # Remove directory
    rm -rf %s/${ID}

    START_TIME=${SECONDS}

    # Run antiSMASH
    %s --allow-long-headers --verbose --skip-zip-file -c %d --output-dir %s/${ID} --html-title ${ID} --taxon %s --minlength %d --databases %s --hmmdetection-strictness %s --logfile %s/${ID}/log.txt --genefinding-gff3 ${GENE_MODELS} ${GENOME}

    # Genbanks to table
    %s -i %s/${ID} -o %s/${ID}/antismash_features.tsv -s %s/${ID}/synopsis.tsv

    END_TIME=${SECONDS}
    RUN_TIME=$((END_TIME-START_TIME))
    echo "*** n=${n} // ${ID} // Duration: ${RUN_TIME} seconds ***"

    n=$(($n+1))

done < %s

# Concatenate tables
%s -a 0 -e %s/*/antismash_features.tsv > %s/biosynthetic_gene_clusters.features.tsv
%s -a 0 -e %s/*/synopsis.tsv > %s/biosynthetic_gene_clusters.type_counts.tsv

"""%( 
    # Args
    output_directory,
    output_directory,

    # antiSMASH
    os.environ["antismash"],
    opts.n_jobs,
    output_directory,
    opts.taxon,
    opts.minimum_contig_length,
    opts.antismash_database,
    opts.hmmdetection_strictness,
    output_directory,

    # Summary table
    os.environ["antismash_genbanks_to_table.py"],
    output_directory,
    output_directory,
    output_directory,


    input_filepaths[0],

    # Concatenate
    os.environ["concatenate_dataframes.py"],
    output_directory,
    output_directory,

    os.environ["concatenate_dataframes.py"],
    output_directory,
    output_directory,
    ),
    ]

    return cmd



# Symlink
def get_symlink_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    
    # Command
    cmd = ["("]
    for filepath in input_filepaths:
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

    accessory_scripts = {"antismash_genbanks_to_table.py", "concatenate_dataframes.py"}

    required_executables={
                # 1
                "antismash",
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
    description = "Organizing genomes and gene models"

    # i/o
    input_filepaths = [opts.mags, opts.gene_models]
    output_filepaths = [
        # os.path.join(directories["project"], "genomes.list"),
        # os.path.join(directories["project"], "gene_models.list"),
        os.path.join(directories["preprocessing"], "genomes_gene-models.tsv"),
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
    # antiSMASH
    # ==========
    step = 1

    program = "antismash"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Biosynthetic gene cluster detection"


    # i/o
    input_filepaths = [
        os.path.join(directories["preprocessing"], "genomes_gene-models.tsv"),
        ]
    output_filenames = ["biosynthetic_gene_clusters.features.tsv", "biosynthetic_gene_clusters.type_counts.tsv"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_antismash_cmd(**params)


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
    step = 2
    description = "Symlinking relevant output files"

    # i/o
    input_filepaths = [
        os.path.join(directories[("intermediate", "1__antismash")], "biosynthetic_gene_clusters.features.tsv"),
        os.path.join(directories[("intermediate", "1__antismash")], "biosynthetic_gene_clusters.type_counts.tsv"),

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
    assert os.path.exists(opts.antismash_database), "--antismash_database %s does not exist".format(opts.antismash_database)
    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -m <mags> -g <gene_models> -o <output_directory> -t bacteria".format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-m","--mags", type=str, help = "Tab-seperated value table of [id_mag]<tab>[path/to/genome.fasta]")
    parser_io.add_argument("-g","--gene_models", type=str, help = "Tab-seperated value table of [id_mag]<tab>[path/to/gene_models.gff]. Must be gff3.")
    parser_io.add_argument("-o","--output_directory", type=str, default="veba_output/cluster", help = "path/to/project_directory [Default: veba_output/cluster]")
    parser_io.add_argument("--mags_extension", type=str, default="fa", help = "Fasta file extension for --mags if a list is provided [Default: fa]")
    parser_io.add_argument("--gene_models_extension", type=str, default="gff", help = "File extension for gene models if a list is provided [Default: gff]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # antiSMASH
    parser_antismash = parser.add_argument_group('antiSMASH arguments')
    parser_antismash.add_argument("-t", "--taxon", type=str, default="bacteria", help="Taxonomic classification of input sequence {bacteria,fungi} [Default: bacteria]")
    parser_antismash.add_argument("--minimum_contig_length", type=int, default=1500, help="Minimum contig length.  [Default: 1500] ")
    parser_antismash.add_argument("-d", "--antismash_database", type=str, default=os.path.join(site.getsitepackages()[0], "antismash", "databases"), help="antiSMASH | Database directory path [Default: {}]".format(os.path.join(site.getsitepackages()[0], "antismash", "databases")))
    parser_antismash.add_argument("-s", "--hmmdetection_strictness", type=str, default="relaxed", help="antiSMASH | Defines which level of strictness to use for HMM-based cluster detection {strict,relaxed,loose}  [Default: relaxed] ")
    parser_antismash.add_argument("--tta_threshold", type=float, default=0.65, help="antiSMASH | Lowest GC content to annotate TTA codons at [Default: 0.65]")
    parser_antismash.add_argument("--antismash_options", type=str, default="", help="antiSMASH | More options (e.g. --arg 1 ) [Default: '']")

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
