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
__version__ = "2023.8.30"

# antiSMASH
def get_antismash_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    cmd = [
"""
OUTPUT_DIRECTORY=%s
INTERMEDIATE_DIRECTORY=%s
mkdir -p ${INTERMEDIATE_DIRECTORY}
n=1
while IFS= read -r LINE
do read -r -a ARRAY <<< $LINE
    ID=${ARRAY[0]}
    GENOME=${ARRAY[1]}
    GENE_MODELS=${ARRAY[2]}

    CHECKPOINT="${INTERMEDIATE_DIRECTORY}/${ID}/ANTISMASH_CHECKPOINT"

    CWD=${PWD}

    if [ ! -f ${CHECKPOINT} ]; then

        echo "[Running ${ID}]"

        # Remove directory
        rm -rf ${INTERMEDIATE_DIRECTORY}/${ID}

        START_TIME=${SECONDS}

        # Run antiSMASH
        %s --allow-long-headers --verbose --skip-zip-file -c %d --output-dir ${INTERMEDIATE_DIRECTORY}/${ID} --html-title ${ID} --taxon %s --minlength %d --databases %s --hmmdetection-strictness %s --logfile ${INTERMEDIATE_DIRECTORY}/${ID}/log.txt --genefinding-gff3 ${GENE_MODELS} ${GENOME}

        # Genbanks to table
        %s -i ${INTERMEDIATE_DIRECTORY}/${ID} -o ${INTERMEDIATE_DIRECTORY}/${ID}/antismash_components.tsv.gz -s ${INTERMEDIATE_DIRECTORY}/${ID}/synopsis.tsv.gz -t ${INTERMEDIATE_DIRECTORY}/${ID}/type_counts.tsv.gz --fasta_output ${INTERMEDIATE_DIRECTORY}/${ID}/bgc.components.faa.gz

        # Compile table for Krona graph
        %s -i ${INTERMEDIATE_DIRECTORY}/${ID}/type_counts.tsv.gz -m biosynthetic-local -o ${INTERMEDIATE_DIRECTORY}/${ID}/krona.tsv

        # Create Krona graph
        %s -o ${INTERMEDIATE_DIRECTORY}/${ID}/krona.html -n ${ID} ${INTERMEDIATE_DIRECTORY}/${ID}/krona.tsv

        # Symlink proteins
        FASTA_DIRECTORY=${OUTPUT_DIRECTORY}/fasta/
        SRC_FILES=${INTERMEDIATE_DIRECTORY}/${ID}/bgc.components.faa.gz
        DST=${FASTA_DIRECTORY}/

        for SRC in $(ls ${SRC_FILES})
        do
            if [ -f "${SRC}" ]; then
                SRC=$(realpath --relative-to ${DST} ${SRC})
                ln -sf ${SRC} ${DST}/${ID}.faa.gz
            fi
        done

        # Symlink genbanks
        GENBANK_DIRECTORY=${OUTPUT_DIRECTORY}/genbanks/
        mkdir -p ${GENBANK_DIRECTORY}/${ID}
        SRC_FILES=${INTERMEDIATE_DIRECTORY}/${ID}/*.region*.gbk
        DST=${GENBANK_DIRECTORY}/${ID}/

        for SRC in $(ls ${SRC_FILES})
        do
            SRC=$(realpath --relative-to ${DST} ${SRC})
            ln -sf ${SRC} ${DST}
        done

        # Completed
        echo "Completed: $(date)" > ${INTERMEDIATE_DIRECTORY}/${ID}/ANTISMASH_CHECKPOINT

        # Remove large assembly files
        rm -rf ${INTERMEDIATE_DIRECTORY}/${ID}/assembly.*

        END_TIME=${SECONDS}
        RUN_TIME=$((END_TIME-START_TIME))
        echo "*** n=${n} // ${ID} // Duration: ${RUN_TIME} seconds ***"

        n=$(($n+1))

    else
        echo "[Skipping ${ID}] Found the following file: ${CHECKPOINT}"
    fi  

done < %s

# Concatenate tables
%s -a 0 -e ${INTERMEDIATE_DIRECTORY}/*/antismash_components.tsv.gz | gzip > ${OUTPUT_DIRECTORY}/bgc.components.tsv.gz
%s -a 0 -e ${INTERMEDIATE_DIRECTORY}/*/type_counts.tsv.gz | gzip > ${OUTPUT_DIRECTORY}/bgc.type_counts.tsv.gz
%s -a 0 -e ${INTERMEDIATE_DIRECTORY}/*/synopsis.tsv.gz | gzip > ${OUTPUT_DIRECTORY}/bgc.synopsis.tsv.gz

# Krona
# Compile table for Krona graph
%s -i ${OUTPUT_DIRECTORY}/bgc.type_counts.tsv.gz -m biosynthetic-global -o ${OUTPUT_DIRECTORY}/krona.tsv

# Create Krona graph
%s ${OUTPUT_DIRECTORY}/krona.tsv -o ${OUTPUT_DIRECTORY}/krona.html -n 'antiSMASH'
"""%( 
    # Args
    directories["output"],
    output_directory,

    # antiSMASH
    os.environ["antismash"],
    opts.n_jobs,
    opts.taxon,
    opts.minimum_contig_length,
    opts.antismash_database,
    opts.hmmdetection_strictness,

    # Summary table
    os.environ["antismash_genbanks_to_table.py"],

    # Krona (Local)
    os.environ["compile_krona.py"],

    os.environ["ktImportText"],

    # Symlink (proteins)

    # Symlink (genbanks)

    # Remove large assembly

    # Input
    input_filepaths[0],

    # Concatenate
    os.environ["concatenate_dataframes.py"],

    os.environ["concatenate_dataframes.py"],

    os.environ["concatenate_dataframes.py"],

    # Krona (Global)
    os.environ["compile_krona.py"],


    os.environ["ktImportText"],

    ),
    ]

    return cmd

# Diamond
def get_diamond_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    fields = ["qseqid","sseqid","pident","length","mismatch","qlen","qstart","qend", "qseq", "slen", "sstart","send", "sseq", "evalue","bitscore","qcovhsp","scovhsp", "cigar"]
    
    # Command
    cmd = [
        "mkdir -p {}".format(os.path.join(directories["tmp"], "diamond")),

            "&&",
        
        "cat",
        os.path.join(input_filepaths[0], "*", "bgc.components.faa.gz"),
        "|",
        "gzip -d",
        ">",
        os.path.join(directories["tmp"], "components.concatenated.faa"),

        # MiBIG
            "&&",

        os.environ["diamond"],
        "blastp",
        "--db {}".format(input_filepaths[1]),
        "--query {}".format(os.path.join(directories["tmp"], "components.concatenated.faa")),
        "--threads {}".format(opts.n_jobs),
        "-f 6 {}".format(" ".join(fields)),
        "--evalue {}".format(opts.diamond_evalue),
        "-o {}".format(os.path.join(output_directory, "diamond_output.mibig.no_header.tsv")),
        "--max-target-seqs 1",
        "--tmpdir {}".format(os.path.join(directories["tmp"], "diamond")),
        # "--header 2",
    ]
    if bool(opts.diamond_sensitivity):
        cmd += [ 
            "--{}".format(opts.diamond_sensitivity),
        ]
    cmd += [ 
        opts.diamond_options,

                "&&",
            
        "echo",
        "'{}'".format("\t".join(fields)),
        ">",
        os.path.join(output_directory, "diamond_output.mibig.tsv"),

            "&&",

        "cat",
        os.path.join(output_directory, "diamond_output.mibig.no_header.tsv"),
        ">>",
        os.path.join(output_directory, "diamond_output.mibig.tsv"),

        # VFDB
                "&&",

        "rm -rf {}".format(os.path.join(directories["tmp"], "diamond", "*")),

                "&&",

        os.environ["diamond"],
        "blastp",
        "--db {}".format(input_filepaths[2]),
        "--query {}".format(os.path.join(directories["tmp"], "components.concatenated.faa")),
        "--threads {}".format(opts.n_jobs),
        "-f 6 {}".format(" ".join(fields)),
        "--evalue {}".format(opts.diamond_evalue),
        "-o {}".format(os.path.join(output_directory, "diamond_output.vfdb.no_header.tsv")),
        "--max-target-seqs 1",
        "--tmpdir {}".format(os.path.join(directories["tmp"], "diamond")),
        # "--header 2",
    ]
    if bool(opts.diamond_sensitivity):
        cmd += [ 
            "--{}".format(opts.diamond_sensitivity),
        ]
    cmd += [ 
        opts.diamond_options,

                "&&",
            
        "echo",
        "'{}'".format("\t".join(fields)),
        ">",
        os.path.join(output_directory, "diamond_output.vfdb.tsv"),

            "&&",

        "cat",
        os.path.join(output_directory, "diamond_output.vfdb.no_header.tsv"),
        ">>",
        os.path.join(output_directory, "diamond_output.vfdb.tsv"),

    ]

    cmd += [ 
    
            "&&",

        os.environ["concatenate_dataframes.py"],
        "-a 1",
        "--prepend_column_levels MiBIG,VFDB",
        "--sort_by bitscore",
        "-o {}".format(os.path.join(directories["output"], "homology.tsv.gz")),
        os.path.join(output_directory, "diamond_output.mibig.tsv"),
        os.path.join(output_directory, "diamond_output.vfdb.tsv"),

            "&&",

        "rm",
        os.path.join(directories["tmp"], "components.concatenated.faa"),
        os.path.join(output_directory, "*.no_header.tsv"),

    ]           

    return cmd

# Novelty
def get_novelty_score_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    
    # Command
    cmd = [
        os.environ["bgc_novelty_scorer.py"],
        "-c {}".format(input_filepaths[0]),
        "-s {}".format(input_filepaths[1]),
        "-d {}".format(input_filepaths[2]),
        "-o {}".format(output_filepaths[0]),
        "--pident {}".format(opts.pident),
        "--qcovhsp {}".format(opts.qcovhsp),
        "--scovhsp {}".format(opts.scovhsp),
        "--evalue {}".format(opts.evalue),
    ]
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
        "antismash_genbanks_to_table.py", 
        "concatenate_dataframes.py",
        "bgc_novelty_scorer.py",
        "compile_krona.py",
        }

    required_executables={
                # 1
                "antismash",
                # 2
                "diamond",
                # Krona
                "ktImportText",

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
        executables[name] = "'{}'".format(os.path.join(opts.script_directory, "scripts", name)) # Can handle spaces in path
                    
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
        opts.input,
        ]
    output_filenames = [
        "bgc.components.tsv.gz", 
        "bgc.type_counts.tsv.gz",
        "bgc.synopsis.tsv.gz",
        "krona.html",
        "fasta/",
        "genbanks/",
        ]
    output_filepaths = list(map(lambda filename: os.path.join(directories["output"], filename), output_filenames))

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

    # ==========
    # Diamond
    # ==========
    step = 2

    program = "diamond"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Diamond [MIBiG, VFDB]"

    # i/o
    input_filepaths = [
         directories[("intermediate",  "1__antismash")], 
         os.path.join(opts.veba_database, "Annotate", "MIBiG", "mibig_v3.1.dmnd"),
         os.path.join(opts.veba_database, "Annotate", "VFDB", "VFDB_setA_pro.dmnd"),
        ]
    output_filenames = ["homology.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(directories["output"], filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_diamond_cmd(**params)

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

    # =============
    # Novelty
    # =============
    program = "novelty"
    # Add to directories
    output_directory = directories["output"]

    # Info
    step = 3
    description = "Calculating novelty score for BGCs and updating synopsis"

    # i/o
    input_filepaths = [
        os.path.join(directories["output"], "bgc.components.tsv.gz"),
        os.path.join(directories["output"], "bgc.synopsis.tsv.gz"),
        os.path.join(directories["output"], "homology.tsv.gz"),
    ]

    output_filenames =  ["bgc.synopsis.tsv.gz"] # Overwriting
    output_filepaths = list(map(lambda fn:os.path.join(directories["output"], fn), output_filenames))


    params = {
    "input_filepaths":input_filepaths,
    "output_filepaths":output_filepaths,
    "output_directory":output_directory,
    "opts":opts,
    "directories":directories,
    }

    cmd = get_novelty_score_cmd(**params)
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

    # # =============
    # # Symlink
    # # =============
    # program = "symlink"
    # # Add to directories
    # output_directory = directories[("output",  "features")] = create_directory(os.path.join(directories["output"], "features"))

    # # Info
    # step = 3
    # description = "Symlinking BGC protein features"

    # # i/o
    # input_filepaths = [
    #     os.path.join(directories[("intermediate",  "1__antismash")], "*", "features.faa.gz"),
    # ]

    # output_filenames =  ["*.bgc_features.faa.gz"] # Overwriting
    # output_filepaths = list(map(lambda fn:os.path.join(directories["output"], fn), output_filenames))


    # params = {
    # "input_filepaths":input_filepaths,
    # "output_filepaths":output_filepaths,
    # "output_directory":output_directory,
    # "opts":opts,
    # "directories":directories,
    # }

    # cmd = get_symlink_cmd(**params)
    # pipeline.add_step(
    #         id=program,
    #         description = description,
    #         step=step,
    #         cmd=cmd,
    #         input_filepaths = input_filepaths,
    #         output_filepaths = output_filepaths,
    #         validate_inputs=True,
    #         validate_outputs=True,
    # )

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
    usage = "{} -i <genomes_gene-models.tsv> -o <output_directory> -t bacteria | Suggested input is from `compile_genomes_table.py` script. Use cut -f3,4,7".format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-i", "--input", type=str, required=True, help = "path/to/input.tsv, Format: Must include the follow columns (No header) [id_mag]<tab>[genome]<tab>[gene_models].  Suggested input is from `compile_genomes_table.py` script. Use cut -f3,4,7")
    # parser_io.add_argument("-m","--mags", type=str, help = "Either 1) Tab-seperated value table of [id_mag]<tab>[path/to/genome.fasta]; or 2) a file with a list filepaths for genomes (Must have same extension)")
    # parser_io.add_argument("-g","--gene_models", type=str, help = "Tab-seperated value table of [id_mag]<tab>[path/to/gene_models.gff]. Must be gff3.")
    parser_io.add_argument("-o","--output_directory", type=str, default="veba_output/biosynthetic", help = "path/to/project_directory [Default: veba_output/biosynthetic]")
    # parser_io.add_argument("--mags_extension", type=str, default="fa", help = "Fasta file extension for --mags if a list is provided [Default: fa]")
    # parser_io.add_argument("--gene_models_extension", type=str, default="gff", help = "File extension for gene models if a list is provided [Default: gff]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Databases
    parser_databases = parser.add_argument_group('Database arguments')
    parser_databases.add_argument("--veba_database", type=str,  help=f"VEBA database location.  [Default: $VEBA_DATABASE environment variable]")
    
    # antiSMASH
    parser_antismash = parser.add_argument_group('antiSMASH arguments')
    parser_antismash.add_argument("-t", "--taxon", type=str, default="bacteria", help="Taxonomic classification of input sequence {bacteria,fungi} [Default: bacteria]")
    parser_antismash.add_argument("--minimum_contig_length", type=int, default=1500, help="Minimum contig length.  [Default: 1500] ")
    parser_antismash.add_argument("-d", "--antismash_database", type=str, default=os.path.join(site.getsitepackages()[0], "antismash", "databases"), help="antiSMASH | Database directory path [Default: {}]".format(os.path.join(site.getsitepackages()[0], "antismash", "databases")))
    parser_antismash.add_argument("-s", "--hmmdetection_strictness", type=str, default="relaxed", help="antiSMASH | Defines which level of strictness to use for HMM-based cluster detection {strict,relaxed,loose}  [Default: relaxed] ")
    parser_antismash.add_argument("--tta_threshold", type=float, default=0.65, help="antiSMASH | Lowest GC content to annotate TTA codons at [Default: 0.65]")
    parser_antismash.add_argument("--antismash_options", type=str, default="", help="antiSMASH | More options (e.g. --arg 1 ) [Default: '']")

    # Diamond
    parser_diamond = parser.add_argument_group('Diamond arguments')
    parser_diamond.add_argument("--diamond_sensitivity", type=str, default="", help="Diamond | Sensitivity [Default:  '']")
    parser_diamond.add_argument("--diamond_evalue", type=float, default=0.001, help="Diamond | E-Value [Default: 0.001]")
    parser_diamond.add_argument("--diamond_options", type=str, default="", help="Diamond | More options (e.g. --arg 1 ) [Default: '']")

    # Novelty
    parser_noveltyscore = parser.add_argument_group('Novelty score threshold arguments')
    parser_noveltyscore.add_argument("--pident", type=float, default=0.0, help = "pident lower bound [float:0 ≤ x < 100] [Default: 0]")
    parser_noveltyscore.add_argument("--qcovhsp", type=float, default=0.0, help = "qcovhsp lower bound [float:0 ≤ x < 100] [Default: 0]")
    parser_noveltyscore.add_argument("--scovhsp", type=float, default=0.0, help = "scovhsp lower bound [float:0 ≤ x < 100] [Default: 0]")
    parser_noveltyscore.add_argument("--evalue", type=float, default=1e-3, help = "e-value lower bound [float:0 < x < 1] [Default: 1e-3]")


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
    directories["project"] = create_directory(opts.output_directory)
    directories["preprocessing"] = create_directory(os.path.join(directories["project"], "preprocessing"))
    directories["output"] = create_directory(os.path.join(directories["project"], "output"))
    directories["log"] = create_directory(os.path.join(directories["project"], "log"))
    directories["tmp"] = create_directory(os.path.join(directories["project"], "tmp"))
    directories["checkpoints"] = create_directory(os.path.join(directories["project"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["project"], "intermediate"))
    os.environ["TMPDIR"] = directories["tmp"]
    directories[("output", "fasta")] = create_directory(os.path.join(directories["output"], "fasta"))
    directories[("output", "genbanks")] = create_directory(os.path.join(directories["output"], "genbanks"))


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
