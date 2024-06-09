#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, site
from collections import OrderedDict, defaultdict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from genopype import __version__ as genopype_version
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.3.5"

# antiSMASH
def get_antismash_from_genomes_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    cmd = [
"""
OUTPUT_DIRECTORY=%s
INTERMEDIATE_DIRECTORY=%s
TMP=%s
mkdir -p ${INTERMEDIATE_DIRECTORY}
n=1
while IFS= read -r LINE
do read -r -a ARRAY <<< $LINE
    ID=${ARRAY[0]}
    GENOME=${ARRAY[1]}
    GENE_MODELS=${ARRAY[2]}

    CHECKPOINT="${INTERMEDIATE_DIRECTORY}/${ID}/ANTISMASH_CHECKPOINT"

    if [ ! -f ${CHECKPOINT} ]; then

        echo "[Running ${ID}]"

        # Remove directory
        rm -rf ${INTERMEDIATE_DIRECTORY}/${ID}

        START_TIME=${SECONDS}

        # Extract CDS records (or else antiSMASH will fail)
        GENE_MODELS_CDS_ONLY=${TMP}/gene_models.cds.gff
        grep "CDS" ${GENE_MODELS} > ${GENE_MODELS_CDS_ONLY}

        # Run antiSMASH
        %s --allow-long-headers --verbose --skip-zip-file -c %d --output-dir ${INTERMEDIATE_DIRECTORY}/${ID} --html-title ${ID} --taxon %s --minlength %d --databases %s --hmmdetection-strictness %s --logfile ${INTERMEDIATE_DIRECTORY}/${ID}/log.txt --genefinding-gff3 ${GENE_MODELS_CDS_ONLY} ${GENOME} || (echo "antiSMASH for ${ID} failed" && exit 1)
        rm ${GENE_MODELS_CDS_ONLY}

        # Genbanks to table
        %s -i ${INTERMEDIATE_DIRECTORY}/${ID} -o ${INTERMEDIATE_DIRECTORY}/${ID}/veba_formatted_output

        # Compile table for Krona graph
        %s -i ${INTERMEDIATE_DIRECTORY}/${ID}/veba_formatted_output/bgc_protocluster-types.tsv.gz -m biosynthetic-local -o ${INTERMEDIATE_DIRECTORY}/${ID}/veba_formatted_output/krona.tsv

        # Create Krona graph
        %s -o ${INTERMEDIATE_DIRECTORY}/${ID}/veba_formatted_output/krona.html -n ${ID} ${INTERMEDIATE_DIRECTORY}/${ID}/veba_formatted_output/krona.tsv

        # Symlink sequences
        FASTA_DIRECTORY=${OUTPUT_DIRECTORY}/fasta/
        DST=${FASTA_DIRECTORY}/

        # Proteins
        SRC_FILES=${INTERMEDIATE_DIRECTORY}/${ID}/veba_formatted_output/fasta/components.faa.gz
        for SRC in $(ls ${SRC_FILES})
        do
            if [ -f "${SRC}" ]; then
                SRC=$(realpath --relative-to ${DST} ${SRC})
                ln -sf ${SRC} ${DST}/${ID}.faa.gz
            fi
        done

        # DNA
        SRC_FILES=${INTERMEDIATE_DIRECTORY}/${ID}/veba_formatted_output/fasta/bgcs.fasta.gz
        for SRC in $(ls ${SRC_FILES})
        do
            if [ -f "${SRC}" ]; then
                SRC=$(realpath --relative-to ${DST} ${SRC})
                ln -sf ${SRC} ${DST}/${ID}.fasta.gz
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
%s -a 0 -e ${INTERMEDIATE_DIRECTORY}/*/veba_formatted_output/identifier_mapping.components.tsv.gz | gzip > ${OUTPUT_DIRECTORY}/identifier_mapping.components.tsv.gz
%s -a 0 -e ${INTERMEDIATE_DIRECTORY}/*/veba_formatted_output/bgc_protocluster-types.tsv.gz | gzip > ${OUTPUT_DIRECTORY}/bgc_protocluster-types.tsv.gz
%s -a 0 -e ${INTERMEDIATE_DIRECTORY}/*/veba_formatted_output/identifier_mapping.bgcs.tsv.gz | gzip > ${OUTPUT_DIRECTORY}/identifier_mapping.bgcs.tsv.gz

# Krona
# Compile table for Krona graph
%s -i ${OUTPUT_DIRECTORY}/bgc_protocluster-types.tsv.gz -m biosynthetic-global -o ${OUTPUT_DIRECTORY}/krona.tsv

# Create Krona graph
%s ${OUTPUT_DIRECTORY}/krona.tsv -o ${OUTPUT_DIRECTORY}/krona.html -n 'antiSMASH'
"""%( 
    # Args
    directories["output"],
    output_directory,
    directories["tmp"],

    # antiSMASH
    os.environ["antismash"],
    opts.n_jobs,
    opts.taxon,
    opts.minimum_contig_length,
    opts.antismash_database,
    opts.hmmdetection_strictness,

    # Summary table
    os.environ["biosynthetic_genbanks_to_table.py"],

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

def get_antismash_from_existing_results_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    cmd = [
"""
OUTPUT_DIRECTORY=%s
INTERMEDIATE_DIRECTORY=%s
TMP=%s
mkdir -p ${INTERMEDIATE_DIRECTORY}
n=1
while IFS= read -r LINE
do read -r -a ARRAY <<< $LINE
    ID=${ARRAY[0]}
    ANTISMASH_RESULTS_DIRECTORY=${ARRAY[1]}

    echo "[Running ${ID}]"

    # Create directory
    mkdir -p ${INTERMEDIATE_DIRECTORY}/${ID}

    START_TIME=${SECONDS}

    cp -rf ${ANTISMASH_RESULTS_DIRECTORY}/*.region*.gbk ${INTERMEDIATE_DIRECTORY}/${ID}

    # Genbanks to table
    %s -i ${INTERMEDIATE_DIRECTORY}/${ID} -o ${INTERMEDIATE_DIRECTORY}/${ID}/veba_formatted_output

    # Compile table for Krona graph
    %s -i ${INTERMEDIATE_DIRECTORY}/${ID}/veba_formatted_output/bgc_protocluster-types.tsv.gz -m biosynthetic-local -o ${INTERMEDIATE_DIRECTORY}/${ID}/veba_formatted_output/krona.tsv

    # Create Krona graph
    %s -o ${INTERMEDIATE_DIRECTORY}/${ID}/veba_formatted_output/krona.html -n ${ID} ${INTERMEDIATE_DIRECTORY}/${ID}/veba_formatted_output/krona.tsv

    # Symlink sequences
    FASTA_DIRECTORY=${OUTPUT_DIRECTORY}/fasta/
    DST=${FASTA_DIRECTORY}/

    # Proteins
    SRC_FILES=${INTERMEDIATE_DIRECTORY}/${ID}/veba_formatted_output/fasta/components.faa.gz
    for SRC in $(ls ${SRC_FILES})
    do
        if [ -f "${SRC}" ]; then
            SRC=$(realpath --relative-to ${DST} ${SRC})
            ln -sf ${SRC} ${DST}/${ID}.faa.gz
        fi
    done

    # DNA
    SRC_FILES=${INTERMEDIATE_DIRECTORY}/${ID}/veba_formatted_output/fasta/bgcs.fasta.gz
    for SRC in $(ls ${SRC_FILES})
    do
        if [ -f "${SRC}" ]; then
            SRC=$(realpath --relative-to ${DST} ${SRC})
            ln -sf ${SRC} ${DST}/${ID}.fasta.gz
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


done < %s

# Concatenate tables
%s -a 0 -e ${INTERMEDIATE_DIRECTORY}/*/veba_formatted_output/identifier_mapping.components.tsv.gz | gzip > ${OUTPUT_DIRECTORY}/identifier_mapping.components.tsv.gz
%s -a 0 -e ${INTERMEDIATE_DIRECTORY}/*/veba_formatted_output/bgc_protocluster-types.tsv.gz | gzip > ${OUTPUT_DIRECTORY}/bgc_protocluster-types.tsv.gz
%s -a 0 -e ${INTERMEDIATE_DIRECTORY}/*/veba_formatted_output/identifier_mapping.bgcs.tsv.gz | gzip > ${OUTPUT_DIRECTORY}/identifier_mapping.bgcs.tsv.gz

# Krona
# Compile table for Krona graph
%s -i ${OUTPUT_DIRECTORY}/bgc_protocluster-types.tsv.gz -m biosynthetic-global -o ${OUTPUT_DIRECTORY}/krona.tsv

# Create Krona graph
%s ${OUTPUT_DIRECTORY}/krona.tsv -o ${OUTPUT_DIRECTORY}/krona.html -n 'antiSMASH'
"""%( 
    # Args
    directories["output"],
    output_directory,
    directories["tmp"],


    # Summary table
    os.environ["biosynthetic_genbanks_to_table.py"],

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
        os.path.join(input_filepaths[0], "*","veba_formatted_output", "fasta", "components.faa.gz"),
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
        "-rf",
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
        "-b {}".format(input_filepaths[1]),
        "-d {}".format(input_filepaths[2]),
        "-o {}".format(output_filepaths[0]),
        "--pident {}".format(opts.pident),
        "--qcovhsp {}".format(opts.qcovhsp),
        "--scovhsp {}".format(opts.scovhsp),
        "--evalue {}".format(opts.evalue),
    ]
    return cmd

# MMSEQS2
def get_mmseqs2_protein_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [
        "cat",
        os.path.join(input_filepaths[0],"*","veba_formatted_output", "fasta", "components.faa.gz"),
        "|",
        "gzip",
        "-d",
        ">",
        os.path.join(directories["tmp"], "components.concatenated.faa"),

            "&&",

        "cat",
        os.path.join(directories["tmp"], "components.concatenated.faa"),
        "|",
        'grep "^>"',
        "|",
        "cut -c2-",
        "|",
        'cut -f1 -d " "', # Should the original gene names be used or the new component id?
        ">",
        os.path.join(directories["tmp"], "components.concatenated.faa.list"),

            "&&",

        os.environ["clustering_wrapper.py"],
        "--fasta {}".format(os.path.join(directories["tmp"], "components.concatenated.faa")),
        "--output_directory {}".format(output_directory),
        "--no_singletons" if bool(opts.no_singletons) else "",
        "--algorithm {}".format(opts.algorithm),
        "--n_jobs {}".format(opts.n_jobs),
        "--minimum_identity_threshold {}".format(opts.protein_minimum_identity_threshold),
        "--minimum_coverage_threshold {}".format(opts.protein_minimum_coverage_threshold),
        "--mmseqs2_options='{}'" if bool(opts.mmseqs2_options) else "",
        "--cluster_prefix {}".format(opts.protein_cluster_prefix),
        "--cluster_suffix {}".format(opts.protein_cluster_suffix) if bool(opts.protein_cluster_suffix) else "",
        "--cluster_prefix_zfill {}".format(opts.protein_cluster_prefix_zfill),
        "--basename clusters",
        "--identifiers {}".format(os.path.join(directories["tmp"], "components.concatenated.faa.list")),
        "--no_sequences_and_header",
        "--representative_output_format fasta",

            "&&",

        "DST={}; SRC={}; SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST/component_clusters.tsv".format(
            directories["output"],
            os.path.join(output_directory, "output/clusters.tsv"),
        ),

            "&&",

        "SRC={}; SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST/components.representative_sequences.faa.gz".format(
            os.path.join(output_directory, "output/representative_sequences.fasta.gz"),
        ),

            "&&",

        "cat",
        os.path.join(output_directory, "output/clusters.tsv"),
        "|",
        """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", header=None); df.insert(0, -1, df[0].map(lambda x: "_".join(x.split("_")[:-1]))); df.to_csv(sys.stdout, sep="\t", index=None, header=None)'""",
        "|",
        os.environ["compile_protein_cluster_prevalence_table.py"],
        "--dtype bool",
        "--rows_name id_bgc",
        "-b",
        "|",
        "gzip",
        ">",
        os.path.join(directories[("output", "prevalence_tables")], "components.tsv.gz"),
    
            "&&",

        "rm -rf",
        os.path.join(directories["tmp"], "components.*"),
  
    ]
    return cmd

def get_mmseqs2_nucleotide_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [
        "cat",
        os.path.join(input_filepaths[0],"*","veba_formatted_output", "fasta", "bgcs.fasta.gz"),
        "|"
        "gzip",
        "-d",
        ">",
        os.path.join(directories["tmp"], "bgcs.concatenated.fasta"),

            "&&",

        "cat",
        os.path.join(directories["tmp"], "bgcs.concatenated.fasta"),
        "|",
        'grep "^>"',
        "|",
        "cut -c2-",
        "|",
        'cut -f1 -d " "', # Should the original gene names be used or the new component id?
        ">",
        os.path.join(directories["tmp"], "bgcs.concatenated.fasta.list"),

            "&&",

        os.environ["clustering_wrapper.py"],
        "--fasta {}".format(os.path.join(directories["tmp"], "bgcs.concatenated.fasta")),
        "--output_directory {}".format(output_directory),
        "--no_singletons" if bool(opts.no_singletons) else "",
        "--algorithm {}".format(opts.algorithm),
        "--n_jobs {}".format(opts.n_jobs),
        "--minimum_identity_threshold {}".format(opts.nucleotide_minimum_identity_threshold),
        "--minimum_coverage_threshold {}".format(opts.nucleotide_minimum_coverage_threshold),
        "--mmseqs2_options='{}'" if bool(opts.mmseqs2_options) else "",
        "--cluster_prefix {}".format(opts.nucleotide_cluster_prefix),
        "--cluster_suffix {}".format(opts.nucleotide_cluster_suffix) if bool(opts.nucleotide_cluster_suffix) else "",
        "--cluster_prefix_zfill {}".format(opts.nucleotide_cluster_prefix_zfill),
        "--basename clusters",
        "--identifiers {}".format(os.path.join(directories["tmp"], "bgcs.concatenated.fasta.list")),
        "--no_sequences_and_header",
        "--representative_output_format fasta",

            "&&",

        "DST={}; SRC={}; SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST/bgc_clusters.tsv".format(
            directories["output"],
            os.path.join(output_directory, "output/clusters.tsv"),
        ),

            "&&",

        "SRC={}; SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST/bgcs.representative_sequences.fasta.gz".format(
            os.path.join(output_directory, "output/representative_sequences.fasta.gz"),
        ),

            "&&",

        "cat",
        os.path.join(output_directory, "output/clusters.tsv"),
        "|",
        """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", header=None); df.insert(0, -1, df[0].map(lambda x: x.split("|")[0])); df.to_csv(sys.stdout, sep="\t", index=None, header=None)'""",
        "|",
        os.environ["compile_protein_cluster_prevalence_table.py"],
        "--dtype bool",
        "--columns_name id_bgc-cluster",
        "-b",
        "|",
        "gzip",
        ">",
        os.path.join(directories[("output", "prevalence_tables")], "bgcs.tsv.gz"),

            "&&",

        "rm -rf",
        os.path.join(directories["tmp"], "bgcs.*"),
  
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
        "biosynthetic_genbanks_to_table.py", 
        "concatenate_dataframes.py",
        "bgc_novelty_scorer.py",
        "compile_krona.py",
        "clustering_wrapper.py",
        "compile_protein_cluster_prevalence_table.py",
        }

    required_executables={
                # 1
                "antismash",
                # 2
                "diamond",
                # Krona
                "ktImportText",
                # 
                "mmseqs",

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

    if opts.from_genomes:

        # Info
        description = "Biosynthetic gene cluster detection from genomes"


        # i/o
        input_filepaths = [
            opts.from_genomes,
            ]

        output_filenames = [
            "identifier_mapping.components.tsv.gz", 
            "bgc_protocluster-types.tsv.gz",
            "identifier_mapping.bgcs.tsv.gz",
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

        cmd = get_antismash_from_genomes_cmd(**params)

    if opts.from_antismash:

        # Info
        description = "Compiling biosynthetic gene clusters from existing antiSMASH run"


        # i/o
        input_filepaths = [
            opts.from_antismash,
            ]

        output_filenames = [
            "identifier_mapping.components.tsv.gz", 
            "bgc_protocluster-types.tsv.gz",
            "identifier_mapping.bgcs.tsv.gz",
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

        cmd = get_antismash_from_existing_results_cmd(**params)


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
        os.path.join(directories["output"], "identifier_mapping.components.tsv.gz"),
        os.path.join(directories["output"], "identifier_mapping.bgcs.tsv.gz"),
        os.path.join(directories["output"], "homology.tsv.gz"),
    ]

    output_filenames =  ["identifier_mapping.bgcs.tsv.gz"] # Overwriting
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

    # ==========
    # MMSEQS2 (Proteins)
    # ==========
    if not opts.no_protein_clustering:
        step += 1

        program = "mmseqs2_protein"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

        # Info
        description = "MMSEQS2 clustering of biosynthetic gene clusters (proteins)"


        # i/o
        input_filepaths = [
            directories[("intermediate",  "1__antismash")],
            ]
        output_filenames = [
            "output/clusters.tsv",
        ]
        if opts.representative_output_format == "table":
            output_filenames += ["output/representative_sequences.tsv.gz"]
        if opts.representative_output_format == "fasta":
            output_filenames += ["output/representative_sequences.fasta.gz"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_mmseqs2_protein_cmd(**params)


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
    # MMSEQS2 (Nucleotide)
    # ==========
    if not opts.no_nucleotide_clustering:
        step += 1

        program = "mmseqs2_nucleotide"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

        # Info
        description = "MMSEQS2 clustering of biosynthetic gene clusters (nucleotides)"


        # i/o
        input_filepaths = [
            directories[("intermediate",  "1__antismash")],
            ]
        output_filenames = [
            "output/clusters.tsv",
        ]
        if opts.representative_output_format == "table":
            output_filenames += ["output/representative_sequences.tsv.gz"]
        if opts.representative_output_format == "fasta":
            output_filenames += ["output/representative_sequences.fasta.gz"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_mmseqs2_nucleotide_cmd(**params)


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
    assert bool(opts.from_genomes) != bool(opts.from_antismash), "Must use either -i/--from_genomes or -A/--from_antismash not both"
    assert any([bool(opts.from_genomes),  bool(opts.from_antismash)]), "Must use either -i/--from_genomes or -A/--from_antismash not both"

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
    parser_io.add_argument("-i", "--from_genomes", type=str, help = "path/to/input.tsv. Cannot be used with -A/--from_antismash. Format: Must include the follow columns (No header) [id_genome]<tab>[genome]<tab>[gene_models].  Suggested input is from `compile_genomes_table.py` script. Use cut -f3,4,7")
    parser_io.add_argument("-A", "--from_antismash", type=str,  help = "path/to/input.tsv. Cannot be used with -i/--from_genomes. Format: Must include the follow columns (No header) [id_genome]<tab>[path/to/biosynthetic_results_directory/]")
    parser_io.add_argument("-o","--output_directory", type=str, default="veba_output/biosynthetic", help = "path/to/project_directory [Default: veba_output/biosynthetic]")


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
    parser_antismash.add_argument("--minimum_contig_length", type=int, default=1, help="Minimum contig length.  [Default: 1] ")
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

    # MMSEQS2
    parser_mmseqs2 = parser.add_argument_group('MMSEQS2 arguments')
    parser_mmseqs2.add_argument("-a", "--algorithm", type=str, default="mmseqs-cluster", choices={"mmseqs-cluster", "mmseqs-linclust"}, help="MMSEQS2 | {mmseqs-cluster, mmseqs-linclust} [Default: mmseqs-cluster]")
    parser_mmseqs2.add_argument("-f","--representative_output_format", type=str, default="fasta", help = "Format of output for representative sequences: {table, fasta} [Default: fasta]") # Should fasta be the new default?


    parser_mmseqs2.add_argument("--protein_minimum_identity_threshold", type=float, default=50.0, help="MMSEQS2 | Protein cluster percent identity threshold (Range (0.0, 100.0]) [Default: 50.0]")
    parser_mmseqs2.add_argument("--protein_minimum_coverage_threshold", type=float, default=0.8, help="MMSEQS2 | Protein coverage threshold (Range (0.0, 1.0]) [Default: 0.8]")
    parser_mmseqs2.add_argument("--protein_cluster_prefix", type=str, default="BGCPC-", help="Protein cluster prefix [Default: 'BGCPC-")
    parser_mmseqs2.add_argument("--protein_cluster_suffix", type=str, default="", help="Protein cluster suffix [Default: '")
    parser_mmseqs2.add_argument("--protein_cluster_prefix_zfill", type=int, default=0, help="Protein cluster prefix zfill. Use 7 to match identifiers from OrthoFinder.  Use 0 to add no zfill. [Default: 0]") #7
    parser_mmseqs2.add_argument("--no_protein_clustering", action="store_true", help="No protein clustering") 

    parser_mmseqs2.add_argument("--nucleotide_minimum_identity_threshold", type=float, default=90.0, help="MMSEQS2 | Nucleotide cluster percent identity threshold (Range (0.0, 100.0]) [Default: 90.0]")
    parser_mmseqs2.add_argument("--nucleotide_minimum_coverage_threshold", type=float, default=0.8, help="MMSEQS2 | Nucleotide coverage threshold (Range (0.0, 1.0]) [Default: 0.8]")
    parser_mmseqs2.add_argument("--nucleotide_cluster_prefix", type=str, default="BGCNC-", help="Nucleotide cluster prefix [Default: 'BGCNC-")
    parser_mmseqs2.add_argument("--nucleotide_cluster_suffix", type=str, default="", help="Nucleotide cluster suffix [Default: '")
    parser_mmseqs2.add_argument("--nucleotide_cluster_prefix_zfill", type=int, default=0, help="Nucleotide cluster prefix zfill. Use 0 to add no zfill. [Default: 0]") #7
    parser_mmseqs2.add_argument("--no_nucleotide_clustering", action="store_true", help="No nucleotide clustering") 

    parser_mmseqs2.add_argument("--no_singletons", action="store_true", help="Exclude singletons (not recommended)") 
    parser_mmseqs2.add_argument("--mmseqs2_options", type=str, default="", help="MMSEQS2 | More options (e.g. --arg 1 ) [Default: '']")


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
    directories[("output", "prevalence_tables")] = create_directory(os.path.join(directories["output"], "prevalence_tables"))


    # Info
    print(format_header(__program__, "="), file=sys.stdout)
    print(format_header("Configuration:", "-"), file=sys.stdout)
    print("Python version:", sys.version.replace("\n"," "), file=sys.stdout)
    print("Python path:", sys.executable, file=sys.stdout) #sys.path[2]
    print("GenoPype version:", genopype_version, file=sys.stdout) #sys.path[2]
    print("Script version:", __version__, file=sys.stdout)
    print("Moment:", get_timestamp(), file=sys.stdout)
    print("Directory:", os.getcwd(), file=sys.stdout)
    if "TMPDIR" in os.environ: print(os.environ["TMPDIR"], file=sys.stdout)
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
