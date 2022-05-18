#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
import pandas as pd
from genopype import * 
from soothsayer_utils import *

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.10.08"

DATABASE="/usr/local/scratch/CORE/jespinoz/db/veba/v2021.10.08"

# Assembly
def preprocess( input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]

    if not os.path.exists(os.path.join(directories["checkpoints"], "0__preprocessing")):
        # Get genomes
        if os.path.isdir(opts.mags):
            filepaths = glob.glob(os.path.join(opts.mags, "*.{}".format(opts.extension)))
            df = pd.Series(filepaths).to_frame("genomes").set_index("genomes")
        else:
            df = pd.read_csv(opts.mags, sep="\t", index_col=0, header=None)
        assert df.shape[1] in {0,1}, "Must either be 1 (for list) or 2 (for table) columns"

        mags = dict()
        if df.shape[1] == 1:
            mags = df.iloc[:,0].dropna().to_dict()
        else:
            for path in df.index:
                id_mag = path.split("/")[-1][:-1*(len(opts.extension) + 1)]
                assert id_mag not in mags, "--mags has non-unique MAG identifiers"
                mags[id_mag] = path 
        mags = pd.Series(mags)
        
        # Write files
        mags.to_frame().to_csv(os.path.join(directories["project"], "genomes.tsv"), sep="\t", header=None)

    
        os.makedirs(os.path.join(output_directory, "genomes"), exist_ok=True)

        # MAGs
        for id_mag, path in mags.items():
            src = os.path.realpath(path)
            dst = os.path.join(os.path.join(output_directory, "genomes", "{}.fa".format(id_mag)))
            os.symlink(src,dst)

        with open(os.path.join(os.path.join(directories["project"], "genomes.list" )), "w") as f:
            for fp in glob.glob(os.path.join(os.path.join(output_directory, "genomes","*.fa" ))):
                print(fp, file=f)


    return []

# MUSCLE
def get_comparesketch_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]
#     # Command

#     cmd = [ # 01:37:07
# """
# for FP in $(cut -f1 %s);
#     do ID=$(basename $FP .fa)
#     OUT_DIR=%s
#     echo $ID
#     %s in=${FP} out=${OUT_DIR}/${ID}.sketch.txt outsketch=${OUT_DIR}/${ID}.sketch overwrite=t format=2 threads=%d %s %s/*.sketch
#     done

# """%( 
#     # comparesketch.sh
#     input_filepaths[0],
#     output_directory,
#     os.environ["comparesketch.sh"],
#     opts.n_jobs,
#     opts.comparesketch_options,
#     os.path.join(opts.database, "Sketches"),
#     )
#     ]
    cmd = [ # 00:35:01
"""

cut -f1 %s | %s -j %d 'ID=$(basename {} .fa) && OUT_DIR=%s && echo "[comparesketch.sh] ${ID}" &&  %s in={} out=${OUT_DIR}/${ID}.txt outsketch=${OUT_DIR}/${ID}.sketch overwrite=t format=2 threads=1 %s %s/*.sketch'

"""%( 
    # comparesketch.sh
    input_filepaths[0],
    os.environ["parallel"],
    opts.n_jobs,
    output_directory,
    os.environ["comparesketch.sh"],
    opts.comparesketch_options,
    os.path.join(opts.database, "Sketches"),
    )
    ]

    return cmd

# Prediction
def get_prediction_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]
    # Command 
    cmd = [ 
        os.environ["predict_domain.py"],
        "-i {}".format(input_filepaths[0]),
        "-o {}".format(output_directory),
        "-d {}".format(os.path.join(opts.database, "NCBITaxonomy", "taxa.sqlite")),
        "-m {}".format(opts.metric),
        "--figsize {}".format(opts.figsize),
        "--colormap {}".format(opts.colormap),
        "--title {}".format(opts.title),
    ]

    return cmd

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
    description = "Symlink genomes"

    # i/o
    input_filepaths = [opts.mags]
    output_filepaths = [
        os.path.join(directories["project"], "genomes.list"),
        os.path.join(directories["project"], "genomes.tsv"),
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
    # comparesketch
    # ==========
    step = 1

    program = "comparesketch"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Compare genome sketches against RefSeq"


    # i/o
    input_filepaths = [
        os.path.join(directories["project"], "genomes.list"),
        ]
    output_filenames = ["*.sketch", "*.txt"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_comparesketch_cmd(**params)

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
    # Predict
    # ==========
    step = 2

    program = "predictions"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Predict domain for each genome"


    # i/o
    input_filepaths = [
        directories[("intermediate",  "1__comparesketch")] ,
        ]
    output_filenames = ["predictions.tsv"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_prediction_cmd(**params)

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

# ============
# Run Pipeline
# ============
# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = set([ 
        "predict_domain.py",
    ])

    
    required_executables={
                "comparesketch.sh",
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
    usage = "{} -i <mags> -o <output> --heatmap_output <pdf> ".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i","--mags", type=str, help = "Directory of MAGs")
    parser_io.add_argument("-x","--extension", type=str, default="fa", help = "Used if directory is specified for --mags [Default: fa]")
    parser_io.add_argument("-o","--output_directory", type=str, default="domain_predictions_output", help = "path/to/output_directory/ [Default: domain_predictions_output]")
    parser_io.add_argument("-d", "--database", type=str, default=DATABASE, help="path/to/database [Default: {}']".format(DATABASE))
    parser_io.add_argument("-m", "--metric", type=str, default="WKID", help="Metric to use for weight.  Choose between {WKID,KID,ANI} [Default: WKID']")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))


    parser_sketch = parser.add_argument_group('Required comparesketch arguments')
    parser_sketch.add_argument("--comparesketch_options", type=str, default="", help="comparesketch | Options [Default:'']")


    parser_plotting = parser.add_argument_group('Required plotting arguments')
    parser_plotting.add_argument("--figsize", type=str, default="infer", help="Matplotlib | figsize [Default: 'infer']")
    parser_plotting.add_argument("--colormap", type=str, default="Greys", help="Matplotlib | colormap [Default: 'Gray']")
    parser_plotting.add_argument("--title", type=str, help="Matplotlib | Heatmap title")


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
  

    # Run pipeline
    with open(os.path.join(directories["project"], "commands.sh"), "w") as f_cmds:
        pipeline = create_pipeline(
                     opts=opts,
                     directories=directories,
                     f_cmds=f_cmds,
        )
        pipeline.compile()
        pipeline.execute(restart_from_checkpoint=opts.restart_from_checkpoint)



    # # Get NCBI Taxonomy server
    # print("Loading NCBI Taxonomy Database:", os.path.join(opts.database, "NCBITaxonomy", "taxa.sqlite"), file=sys.stderr)
    # ncbi = NCBITaxa(
    #     dbfile=os.path.join(opts.database, "NCBITaxonomy", "taxa.sqlite"), 
    # )

    # # Domain from lineage
    # def get_domain_from_taxonid(id_taxon, ncbi=ncbi):
    #     # Force integer
    #     id_taxon = int(id_taxon)

    #     try:
    #         # Get lineage in taxon ids
    #         lineage = ncbi.get_lineage(id_taxon)

    #         # Get ranks for taxon ids
    #         ranks = dict(filter(lambda x: x[1] != "no rank", ncbi.get_rank(lineage).items()))

    #         # Translate ranks
    #         lineage_translated = pd.Series(ncbi.get_taxid_translator(ranks.keys()))
    #         lineage_translated.index = lineage_translated.index.map(lambda x:ranks[x])

    #         # Return superkingdom
    #         return lineage_translated["superkingdom"]
    #     except ValueError:
    #         return np.nan

    # # Get scores
    # mag_to_weights = dict()
    # for fp in glob.glob("./sketches/*.out"):
    #     id_mag = fp.split("/")[-1][:-11]
    #     domain_to_weight = {
    #         "Archaea":0.0,
    #         "Bacteria":0.0,
    #         "Eukaryota":0.0,
    #         "Viruses":0.0,
    #     }
    #     df = pd.read_csv(fp, sep="\t",  skiprows=2)
    #     df["Domain"] = df["TaxID"].map(lambda id_taxon: get_domain_from_taxonid(id_taxon))
    #     df[metric] = df[metric].map(lambda x: float(x[:-1]))
    #     for i, row in df.iterrows():
    #         id_domain = row["Domain"]
    #         weight = row[metric]
    #         domain_to_weight[id_domain] += weight
    #     mag_to_weights[id_mag] = pd.Series(domain_to_weight)
    # df_weights = pd.DataFrame(mag_to_weights).T
    # df_weights.index.name = "Genomes"

    # # Get predictions
    # print("Calculating probabilities for each superkingdom", file=sys.stderr)
    # Y_hat = df_scores/df_scores.sum(axis=1).values.reshape((-1,1))
    # y_hat = Y_hat.idxmax(axis=1)

    # # Add multiindex
    # df_scores.columns = df_scores.columns.map(lambda x: ("Scores", x))
    # Y_hat.columns = Y_hat.columns.map(lambda x: ("Probabilities", x))

    # # Output
    # print("Creating prediction output file:", opts.output, file=sys.stderr)
    # df_output = pd.concat([
    #     y_hat.to_frame(("Prediction", "Domain")),
    #     Y_hat,
    #     df_scores,
    # ], axis=1)
    # df_output.index.name = "Genomes"

    # df_output.to_csv({True:sys.stdout, False:opts.output}[opts.output == "stdout"], sep="\t")

    # # Prokaryotes
    # print("Creating Prokaryote list file:", opts.prokaryotes, file=sys.stderr)
    # with open(os.path.join(opts.output_directory,"prokaryotes.list") "w") as f:
    #     for id_mag in y_hat[y_hat.map(lambda y: y in {"Archaea", "Bacteria"})].index:
    #         print(id_mag, file=f)

    # # Eukaryotes
    # print("Creating Eukaryote list file:", opts.eukaryotes, file=sys.stderr)
    # with open(os.path.join(opts.output_directory,"eukaryotes.list"), "w") as f:
    #     for id_mag in y_hat[y_hat.map(lambda y: y in { "Eukaryota"})].index:
    #         print(id_mag, file=f)

    # # Viruses
    # print("Creating Virus list file:", opts.viruses, file=sys.stderr)
    # with open(os.path.join(opts.output_directory,"viruses.list"), "w") as f:
    #     for id_mag in y_hat[y_hat.map(lambda y: y in { "Viruses"})].index:
    #         print(id_mag, file=f)

    # # Heatmap
    # if opts.heatmap_output:
    #     print("Creating heatmap of prediction probabilities:", opts.heatmap_output, file=sys.stderr)

    #     import matplotlib.pyplot as plt 
    #     import seaborn as sns 

    #     df_heatmap = df_output["Probabilities"]
    #     if opts.figsize == "infer":
    #         n, m = df_heatmap.shape
    #         opts.figsize = (m*1.618, n*0.5)
    #     else:
    #         opts.figsize = opts.figsize.strip().split(",")
    #         opts.figsize = (float(opts.figsize[0]), float(opts.figsize([1])))

    #     with plt.style.context("seaborn-white"):
    #         fig, ax = plt.subplots(figsize=opts.figsize)
    #         sns.heatmap(df_heatmap, vmin=0, vmax=1, cmap=opts.colormap, annot=True, ax=ax, edgecolor="white", linewidth=1, cbar_kws={"label":"Probability"})
    #         ax.set_xticklabels(ax.get_xticklabels(), fontsize=15)
    #         ax.set_yticklabels(ax.get_yticklabels(), fontsize=12)
    #         ax.set_xlabel("Superkingdom", fontsize=15, fontweight="bold")
    #         ax.set_ylabel("Genomes", fontsize=15, fontweight="bold")
    #         if opts.title:
    #             ax.set_title(opts.title, fontsize=15, fontweight="bold")
    #         fig.savefig(opts.heatmap_output, format="pdf", bbox_inches="tight", dpi=300)

if __name__ == "__main__":
    main()
