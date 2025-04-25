#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, shutil
from collections import OrderedDict, defaultdict

import pandas as pd
import numpy as np

# Soothsayer Ecosystem
from genopype import *
from genopype import __version__ as genopype_version
from soothsayer_utils import *

# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2025.4.25"

# Tiara
def get_tiara_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [
    os.environ["tiara"],
    "-i {}".format(opts.fasta),
    "-o {}".format(os.path.join(output_directory,  "tiara_output.tsv")),
    "--probabilities",
    "-m {}".format(opts.tiara_minimum_length),
    "-t {}".format(opts.n_jobs),
    opts.tiara_options,
    ]
    return cmd

def get_partition_organelle_sequences_single_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [
    os.environ["partition_organelle_sequences.py"],
    "-f {}".format(input_filepaths[0]),
    "-t {}".format(input_filepaths[1]),
    "-n {}".format(opts.name),
    "-o {}".format(directories["output"]),
    "-u {}".format(opts.unknown_organelle_prediction),
    # "--mitochondrion_suffix {}".format(opts.mitochondrion_suffix),
    # "--plastid_suffix {}".format(opts.plastid_suffix),
    # "--unknown_suffix {}".format(opts.unknown_suffix),
    "--verbose",

        "&&",

    "echo",
    opts.name, 
    ">",
    os.path.join(output_directory, "genomes.list"),
"""
OUTPUT_DIRECTORY={}
INTERMEDIATE_DIRECTORY={}
cat $OUTPUT_DIRECTORY/*.fa | grep "^>" | cut -c2- | cut -f1 -d " " > $INTERMEDIATE_DIRECTORY/eukaryotic_contigs.list
cat $OUTPUT_DIRECTORY/mitochondrion/*.fa | grep "^>" | cut -c2- | cut -f1 -d " " > $INTERMEDIATE_DIRECTORY/mitochondria_contigs.list
cat $OUTPUT_DIRECTORY/plastid/*.fa | grep "^>" | cut -c2- | cut -f1 -d " " > $INTERMEDIATE_DIRECTORY/plastid_contigs.list
""".format(
    directories["output"],
    output_directory,
)
    ]
    return cmd

def get_partition_organelle_sequences_multiple_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [
"""
OUTPUT_DIRECTORY={}
INTERMEDIATE_DIRECTORY={}
S2B={}

# Partition sequences
for ID in $(cut -f2 $S2B | sort -u);
do
    SCAFFOLD_LIST={}
    # grep "$ID$" $S2B | cut -f1 > $SCAFFOLD_LIST # The extra $ for a linebreak is needed to distinguish betwee overlapping ids (https://github.com/jolespin/veba/issues/175)
    awk -F'\\t' -v id="$ID" '$2 == id {{print $1}}' $S2B > $SCAFFOLD_LIST
    
    {} grep -f $SCAFFOLD_LIST {} | {} -f stdin -t {} -n $ID -o $OUTPUT_DIRECTORY -u {} --verbose
done

# Store list of genomes
cat $S2B | cut -f2 | sort -u > $INTERMEDIATE_DIRECTORY/genomes.list

# Store list of partitioned contigs
cat $OUTPUT_DIRECTORY/*.fa | grep "^>" | cut -c2- | cut -f1 -d " " > $INTERMEDIATE_DIRECTORY/eukaryotic_contigs.list
cat $OUTPUT_DIRECTORY/mitochondrion/*.fa | grep "^>" | cut -c2- | cut -f1 -d " " > $INTERMEDIATE_DIRECTORY/mitochondrion_contigs.list
cat $OUTPUT_DIRECTORY/plastid/*.fa | grep "^>" | cut -c2- | cut -f1 -d " " > $INTERMEDIATE_DIRECTORY/plastid_contigs.list

rm -rf $SCAFFOLD_LIST
""".format(
        directories["output"],
        output_directory,
        opts.scaffolds_to_bins,
        os.path.join(directories["tmp"], "$ID.list"),
        os.environ["seqkit"],
        input_filepaths[0],

        os.environ["partition_organelle_sequences.py"],
        input_filepaths[1],
        opts.unknown_organelle_prediction,
        # opts.mitochondrion_suffix,
        # opts.plastid_suffix,
        # opts.unknown_suffix,

)
    ]
    return cmd

# MetaEuk
def get_metaeuk_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [
        # Get eukaryotic sequences
        "cat",
        opts.fasta,
        "|",
        os.environ["seqkit"],
        "grep",
        "-f {}".format(input_filepaths[0]),
        ">",
        os.path.join(directories["tmp"], "tmp.fasta"),


            "&&",

        # Placeholder
        "OUTPUT_DIRECTORY={};  for ID in $(cat {}); do >$OUTPUT_DIRECTORY/$ID.faa; >$OUTPUT_DIRECTORY/$ID.ffn; >$OUTPUT_DIRECTORY/$ID.gff; done".format(output_directory, input_filepaths[1]),

            "&&",

        # Run MetaEuk
        os.environ["metaeuk"],
        "easy-predict",
        "--threads {}".format(opts.n_jobs),
        "-s {}".format(opts.metaeuk_sensitivity),
        "-e {}".format(opts.metaeuk_evalue),
        "--split-memory-limit {}".format(opts.metaeuk_split_memory_limit),
        opts.metaeuk_options,
        os.path.join(directories["tmp"], "tmp.fasta"),
        opts.metaeuk_database, # db
        os.path.join(output_directory, "metaeuk"), # output prefix
        os.path.join(directories["tmp"],"metaeuk"),

        # Convert MetaEuk identifiers
            "&&",

        os.environ["compile_metaeuk_identifiers.py"],
        "--cds {}".format(os.path.join(output_directory, "metaeuk.codon.fas")),
        "--protein {}".format(os.path.join(output_directory, "metaeuk.fas")),
        "-o {}".format(output_directory),
        "-b {}".format(opts.basename),
    ]

    # Remove temporary files
    cmd += [

        "&&",
        
    "rm -rf {} {} {}".format(
        os.path.join(output_directory, "*.fas"), # output prefix
        os.path.join(output_directory, "metaeuk.gff"), # output prefix
        os.path.join(directories["tmp"],"metaeuk", "*"),
        os.path.join(directories["tmp"], "tmp.fasta"),
    ),
    ]

    if opts.scaffolds_to_bins:
        cmd += [
        # Partition the gene models and genomes
            "&&",

        os.environ["partition_gene_models.py"],
        "-i {}".format(opts.scaffolds_to_bins),
        # "-f {}".format(opts.fasta),
        "-g {}".format(os.path.join(output_directory, "{}.gff".format(opts.basename))),
        "-d {}".format(os.path.join(output_directory, "{}.ffn".format(opts.basename))),
        "-a {}".format(os.path.join(output_directory, "{}.faa".format(opts.basename))),
        "-o {}".format(os.path.join(output_directory)),
        "--use_mag_as_description",

            "&&",

        "rm -rf",
        os.path.join(output_directory, "{}.*".format(opts.basename)),

        ]
    return cmd

# Pyrodigal
def get_pyrodigal_cmd(input_filepaths, output_filepaths, output_directory, directories, opts, genetic_code):

    cmd = [

        # Placeholder
        "OUTPUT_DIRECTORY={};  for ID in $(cat {}); do >$OUTPUT_DIRECTORY/$ID.faa; >$OUTPUT_DIRECTORY/$ID.ffn; >$OUTPUT_DIRECTORY/$ID.gff; done".format(output_directory, input_filepaths[1]),

            "&&",

        "cat",
        opts.fasta,
        "|",
        os.environ["seqkit"],
        "grep",
        "-f {}".format(input_filepaths[0]),
        "|",

        # Run analysis
        os.environ["pyrodigal"],
        "-p meta",
        "-g {}".format(genetic_code),
        "-f gff",
        "-d {}".format(os.path.join(output_directory, "{}.ffn".format(opts.basename))),
        "-a {}".format(os.path.join(output_directory, "{}.faa".format(opts.basename))),
        "--min-gene {}".format(opts.pyrodigal_minimum_gene_length),
        "--min-edge-gene {}".format(opts.pyrodigal_minimum_edge_gene_length),
        "--max-overlap {}".format(opts.pyrodigal_maximum_gene_overlap_length),
        "-j {}".format(opts.n_jobs),
        "|",
        os.environ["append_geneid_to_prodigal_gff.py"],
        "-a gene_id",
        ">",
        os.path.join(output_directory, "{}.gff".format(opts.basename)),

        ]
    
    if opts.scaffolds_to_bins:
        cmd += [
        # Partition the gene models and genomes
            "&&",

        os.environ["partition_gene_models.py"],
        "-i {}".format(opts.scaffolds_to_bins),
        # "-f {}".format(opts.fasta),
        "-g {}".format(os.path.join(output_directory, "{}.gff".format(opts.basename))),
        "-d {}".format(os.path.join(output_directory, "{}.ffn".format(opts.basename))),
        "-a {}".format(os.path.join(output_directory, "{}.faa".format(opts.basename))),
        "-o {}".format(os.path.join(output_directory)),
        "--use_mag_as_description",

            "&&",

        "rm -rf",
        os.path.join(output_directory, "{}.*".format(opts.basename)),

        ]

    return cmd

# barrnap
def get_barrnap_cmd(input_filepaths, output_filepaths, output_directory, directories, opts, kingdom):
    cmd = [

"""
# Placeholder
OUTPUT_DIRECTORY={}
for ID in $(cat {})
do 
    >$OUTPUT_DIRECTORY/$ID.rRNA
    >$OUTPUT_DIRECTORY/$ID.rRNA.gff
done

# Run analysis  
for GENOME_FASTA in {}
do  
    ID=$(basename $GENOME_FASTA .fa)
    {} --kingdom {} --threads {} --lencutoff {} --reject {} --evalue {} --outseq $OUTPUT_DIRECTORY/$ID.rRNA $GENOME_FASTA | {} > $OUTPUT_DIRECTORY/$ID.rRNA.gff
    rm -rf $GENOME_FASTA.fai
done
""".format(
        output_directory,
        input_filepaths[0],

        input_filepaths[1],
        os.environ["barrnap"],
        kingdom,
        opts.n_jobs,
        opts.barrnap_length_cutoff, 
        opts.barrnap_reject, 
        opts.barrnap_evalue,
        os.environ["append_geneid_to_barrnap_gff.py"],
    ),
    ]
    return cmd

# tRNAscan-SE
def get_trnascan_cmd(input_filepaths, output_filepaths, output_directory, directories, opts, search_mode, trnascan_options):
    cmd = [

"""

# Placeholder
OUTPUT_DIRECTORY={}
for ID in $(cat {})
do 
    >$OUTPUT_DIRECTORY/$ID.tRNA
    >$OUTPUT_DIRECTORY/$ID.tRNA.gff
    >$OUTPUT_DIRECTORY/$ID.tRNA.struct
    >$OUTPUT_DIRECTORY/$ID.tRNA.txt

done

# Run analysis  
for GENOME_FASTA in {}
do  
    ID=$(basename $GENOME_FASTA .fa)
    TRNA_FASTA=$OUTPUT_DIRECTORY/$ID.tRNA
    if [[ -s "$TRNA_FASTA" ]]; 
        then
            echo "[Skipping] [tRNAscan-SE] $GENOME_FASTA because tRNA fasta exists and is not empty"
        else
            echo "[Running] [tRNAscan-SE] $GENOME_FASTA"
            {} {} --forceow --progress --thread {} --fasta $OUTPUT_DIRECTORY/$ID.tRNA --gff $OUTPUT_DIRECTORY/$ID.tRNA.gff --struct $OUTPUT_DIRECTORY/$ID.tRNA.struct {} $GENOME_FASTA > $OUTPUT_DIRECTORY/$ID.tRNA.txt
    fi
done
""".format(
        output_directory,
        input_filepaths[0],

        input_filepaths[1],
        os.environ["tRNAscan-SE"],
        search_mode,
        opts.n_jobs,
        trnascan_options,
    ),
    ]
    return cmd

def get_symlink_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    # Nuclear
    cmd = [
    "DST={}; SRC_DIRECTORY={}; (for SRC in $SRC_DIRECTORY/*.faa $SRC_DIRECTORY/*.ffn $SRC_DIRECTORY/*.ffn $SRC_DIRECTORY/identifier_mapping.metaeuk.tsv $SRC_DIRECTORY/identifier_mapping.tsv; do SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST; done)".format(
        output_directory,
        directories[("intermediate", "2__metaeuk")],
        ),
        
        "&&",

    "DST={}; SRC_DIRECTORY={}; (for SRC in $SRC_DIRECTORY/*.rRNA; do SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST; done)".format(
        output_directory,
        directories[("intermediate", "5__barrnap-nuclear")],
        ),

        "&&",

    "DST={}; SRC_DIRECTORY={}; (for SRC in $SRC_DIRECTORY/*.tRNA; do SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST; done)".format(
        output_directory,
        directories[("intermediate", "8__trnascan-nuclear")],
        ),

        "&&",

    "for ID in $(cat {}); do cat {}/$ID.gff {}/$ID.rRNA.gff {}/$ID.tRNA.gff > {}/$ID.gff; done".format(
        input_filepaths[0],
        directories[("intermediate", "2__metaeuk")],
        directories[("intermediate", "5__barrnap-nuclear")],
        directories[("intermediate", "8__trnascan-nuclear")],
        output_directory,
    )


    ]

    # Mitochondrion
    cmd += [
        "&&",

    "DST={}; SRC_DIRECTORY={}; (for SRC in $SRC_DIRECTORY/*.faa $SRC_DIRECTORY/*.ffn $SRC_DIRECTORY/identifier_mapping.tsv; do SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST; done)".format(
        os.path.join(output_directory,"mitochondrion"),
        directories[("intermediate", "3__pyrodigal-mitochondrion")],
        ),

       "&&",

    "DST={}; SRC_DIRECTORY={}; (for SRC in $SRC_DIRECTORY/*.rRNA; do SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST; done)".format(
        os.path.join(output_directory,"mitochondrion"),
        directories[("intermediate", "6__barrnap-mitochondrion")],
        ),

        "&&",

    "DST={}; SRC_DIRECTORY={}; (for SRC in $SRC_DIRECTORY/*.tRNA; do SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST; done)".format(
        os.path.join(output_directory,"mitochondrion"),
        directories[("intermediate", "9__trnascan-mitochondrion")],
        ),

            "&&",

    "for ID in $(cat {}); do cat {}/$ID.gff {}/$ID.rRNA.gff {}/$ID.tRNA.gff > {}/$ID.gff; done".format(
        input_filepaths[0],
        directories[("intermediate", "3__pyrodigal-mitochondrion")],
        directories[("intermediate", "6__barrnap-mitochondrion")],
        directories[("intermediate", "9__trnascan-mitochondrion")],
        os.path.join(output_directory, "mitochondrion"),
    )
    ]

    # Plastid
    cmd += [
        "&&",

    "DST={}; SRC_DIRECTORY={}; (for SRC in $SRC_DIRECTORY/*.faa $SRC_DIRECTORY/*.ffn $SRC_DIRECTORY/identifier_mapping.tsv; do SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST; done)".format(
        os.path.join(output_directory,"plastid"),
        directories[("intermediate", "4__pyrodigal-plastid")],
        ),


       "&&",

    "DST={}; SRC_DIRECTORY={}; (for SRC in $SRC_DIRECTORY/*.rRNA; do SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST; done)".format(
        os.path.join(output_directory,"plastid"),
        directories[("intermediate", "7__barrnap-plastid")],
        ),

        "&&",

    "DST={}; SRC_DIRECTORY={}; (for SRC in $SRC_DIRECTORY/*.tRNA; do SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST; done)".format(
        os.path.join(output_directory,"plastid"),
        directories[("intermediate", "10__trnascan-plastid")],
        ),


            "&&",

    "for ID in $(cat {}); do cat {}/$ID.gff {}/$ID.rRNA.gff {}/$ID.tRNA.gff > {}/$ID.gff; done".format(
        input_filepaths[0],
        directories[("intermediate", "4__pyrodigal-plastid")],
        directories[("intermediate", "7__barrnap-plastid")],
        directories[("intermediate", "10__trnascan-plastid")],
        os.path.join(output_directory, "plastid"),
    )
    ]
    
    return cmd

def get_stats_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

 # Genomes
    cmd = [ 

        os.environ["seqkit"],
        "stats",
        "-a",
        "-T",
        "-j {}".format(opts.n_jobs),
        os.path.join(output_directory, "*.fa"),
        os.path.join(output_directory, "*", "*.fa"),

        "|",
        # """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: "/".join(x.split("/")[2:])[:-3]); df.to_csv(sys.stdout, sep="\t")'""",
        """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: x.split("output/")[-1][:-3]); df.to_csv(sys.stdout, sep="\t")'""",

        ">",
        os.path.join(output_directory,"genome_statistics.tsv"),

            
        # CDS

            "&&",

        os.environ["seqkit"],
        "stats",
        "-a",
        # "-b",
        "-T",
        "-j {}".format(opts.n_jobs),
        os.path.join(output_directory, "*.ffn"),
        os.path.join(output_directory,"*",  "*.ffn"),

        "|",
        # """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: "/".join(x.split("/")[2:])[:-4]); df.to_csv(sys.stdout, sep="\t")'""",
            """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: x.split("output/")[-1][:-4]); df.to_csv(sys.stdout, sep="\t")'""",

        ">",
        os.path.join(output_directory,"gene_statistics.cds.tsv"),


        # rRNA
            "&&",

        os.environ["seqkit"],
        "stats",
        "-a",
        # "-b",
        "-T",
        "-j {}".format(opts.n_jobs),
        os.path.join(output_directory, "*.rRNA"),
        os.path.join(output_directory,"*",  "*.rRNA"),

        "|",
        # """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: "/".join(x.split("/")[2:])[:-5]); df.to_csv(sys.stdout, sep="\t")'""",
        """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: x.split("output/")[-1][:-5]); df.to_csv(sys.stdout, sep="\t")'""",

        ">",
        os.path.join(output_directory,"gene_statistics.rRNA.tsv"),


        # tRNA
            "&&",

        os.environ["seqkit"],
        "stats",
        "-a",
        # "-b",
        "-T",
        "-j {}".format(opts.n_jobs),
        os.path.join(output_directory, "*.tRNA"),
        os.path.join(output_directory,"*",  "*.tRNA"),

        "|",
        # """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: "/".join(x.split("/")[2:])[:-5]); df.to_csv(sys.stdout, sep="\t")'""",
        """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: x.split("output/")[-1][:-5]); df.to_csv(sys.stdout, sep="\t")'""",

        ">",
        os.path.join(output_directory,"gene_statistics.tRNA.tsv"),
            

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

    # ============
    # Partitioning
    # ============
    if  opts.tiara_results:
        output_filepaths = [opts.tiara_results]
    else:

        step = 0

        program = "tiara"

        program_label = "{}__{}".format(step, program)
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


        # Info
        description = "Predict taxonmic domain of contigs"
        
        # i/o
        input_filepaths = [
            opts.fasta, 
            ]

        output_filepaths = [ 
            os.path.join(output_directory, "tiara_output.tsv"),
        ]


        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_tiara_cmd(**params)


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



    # ==========
    # Partition sequences
    # ==========
    step = 1

    program = "partition"

    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Partitioning nuclear, mitochondrial, and plastid sequences"
    
    # i/o

    input_filepaths = [
        opts.fasta,
        output_filepaths[0],
    ]



    output_filepaths = [
        os.path.join(directories["output"],"*.fa"),
        os.path.join(output_directory,"eukaryotic_contigs.list"),
        os.path.join(output_directory,"genomes.list"),

    ]

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    if opts.scaffolds_to_bins:
        cmd = get_partition_organelle_sequences_multiple_cmd(**params)
    else:
        cmd = get_partition_organelle_sequences_single_cmd(**params)

    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
                errors_ok=False,
                acceptable_returncodes={0},                    
                log_prefix=program_label,
                # whitelist_empty_output_files=[
                #     "mitochondrion/*.fa", 
                #     "plastid/*.fa",
                # ]
    )

    # =============
    # MetaEuk
    # =============
    step = 2

    program = "metaeuk"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))
    # output_directory = directories["output"]

    # Info
    description = "Ab initio eukaryotic gene prediction"

    input_filepaths = [
            os.path.join( directories[("intermediate",  "1__partition")], "eukaryotic_contigs.list"),
            os.path.join( directories[("intermediate",  "1__partition")], "genomes.list"),
            opts.fasta,

        ]
    if opts.scaffolds_to_bins:
        input_filepaths += [
                opts.scaffolds_to_bins,
            ]


    output_filenames = [
        "*.fa",
        "*.faa",
        "*.gff",
        "*.ffn",
        "identifier_mapping.metaeuk.tsv", 
        "metaeuk.headersMap.tsv",
    ]

    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_metaeuk_cmd(**params)

    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=True,
                errors_ok=False,
                log_prefix=program_label,

    )

   # =============
    # Pyrodigal
    # =============
    step = 3

    program = "pyrodigal-mitochondrion"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))
    # output_directory = os.path.join(directories["output"], "mitochondrion")


    # Info
    description = "Ab initio prokaryotic gene prediction [Mitochondrion]"

    input_filepaths = [
            os.path.join( directories[("intermediate",  "1__partition")], "mitochondrion_contigs.list"),
            os.path.join( directories[("intermediate",  "1__partition")], "genomes.list"),
            opts.fasta,
        ]
    if opts.scaffolds_to_bins:
        input_filepaths += [
                opts.scaffolds_to_bins,
            ]


    output_filenames = [
        "*.fa",
        "*.faa",
        "*.gff",
        "*.ffn",
    ]

    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        "genetic_code":opts.pyrodigal_mitochondrial_genetic_code,
    }

    cmd = get_pyrodigal_cmd(**params)

    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=False,
                validate_outputs=False,
                errors_ok=False,
                log_prefix=program_label,

    )

   # =============
    # Pyrodigal
    # =============
    step = 4

    program = "pyrodigal-plastid"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))
    # output_directory = os.path.join(directories["output"], "plastid")

    # Info
    description = "Ab initio prokaryotic gene prediction [Plastid]"

    input_filepaths = [
            os.path.join( directories[("intermediate",  "1__partition")], "plastid_contigs.list"),
            os.path.join( directories[("intermediate",  "1__partition")], "genomes.list"),
            opts.fasta,

        ]
    if opts.scaffolds_to_bins:
        input_filepaths += [
                opts.scaffolds_to_bins,
            ]


    output_filenames = [
        "*.fa",
        "*.faa",
        "*.gff",
        "*.ffn",
    ]

    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        "genetic_code":opts.pyrodigal_plastid_genetic_code,
    }

    cmd = get_pyrodigal_cmd(**params)

    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=False,
                validate_outputs=False,
                errors_ok=False,
                log_prefix=program_label,

    )

   # =============
    # BARRNAP
    # =============
    step = 5

    program = "barrnap-nuclear"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "rRNA gene detection [Nuclear]"

    input_filepaths = [
            os.path.join( directories[("intermediate",  "1__partition")], "genomes.list"),
            os.path.join( directories["output"],  "*.fa"),
        ]

    output_filenames = [
        "*.rRNA",
    ]

    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        "kingdom":"euk",
    }

    cmd = get_barrnap_cmd(**params)

    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=False,
                validate_outputs=False,
                errors_ok=False,
                log_prefix=program_label,

    )

    # =============
    # BARRNAP
    # =============
    step = 6

    program = "barrnap-mitochondrion"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "rRNA gene detection [Mitochondrion]"

    input_filepaths = [
            os.path.join( directories[("intermediate",  "1__partition")], "genomes.list"),
            os.path.join( directories["output"], "mitochondrion", "*.fa"),

        ]

    output_filenames = [
        "*.rRNA",
    ]

    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        "kingdom":"mito",
    }

    cmd = get_barrnap_cmd(**params)

    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=False,
                validate_outputs=False,
                errors_ok=False,
                log_prefix=program_label,

    )

    # =============
    # BARRNAP
    # =============
    step = 7

    program = "barrnap-plastid"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "rRNA gene detection [Plastid]"

    input_filepaths = [
            os.path.join( directories[("intermediate",  "1__partition")], "genomes.list"),
            os.path.join( directories["output"], "plastid", "*.fa"),
        ]

    output_filenames = [
        "*.rRNA",
    ]

    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        "kingdom":"bac",
    }

    cmd = get_barrnap_cmd(**params)

    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=False,
                validate_outputs=False,
                errors_ok=False,
                log_prefix=program_label,

    )

    # =============
    # tRNAscan-SE
    # =============
    step = 8

    program = "trnascan-nuclear"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "tRNA gene detection [Nuclear]"

    input_filepaths = [
            os.path.join( directories[("intermediate",  "1__partition")], "genomes.list"),
            os.path.join( directories["output"], "*.fa"),
        ]

    output_filenames = [
        "*.tRNA",
    ]

    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        "search_mode":"-E",
        "trnascan_options":opts.trnascan_nuclear_options,
    }

    cmd = get_trnascan_cmd(**params)

    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=False,
                validate_outputs=False,
                errors_ok=False,
                log_prefix=program_label,

    )

    # =============
    # tRNAscan-SE
    # =============
    step = 9

    program = "trnascan-mitochondrion"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "tRNA gene detection [Mitochondrion]"

    input_filepaths = [
            os.path.join( directories[("intermediate",  "1__partition")], "genomes.list"),
            os.path.join( directories["output"], "mitochondrion", "*.fa"),
        ]

    output_filenames = [
        "*.tRNA",
    ]

    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        "search_mode":opts.trnascan_mitochondrial_searchmode,
        "trnascan_options":opts.trnascan_mitochondrial_options,

    }

    cmd = get_trnascan_cmd(**params)

    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=False,
                validate_outputs=False,
                errors_ok=False,
                log_prefix=program_label,

    )

    # =============
    # tRNAscan-SE
    # =============
    step = 10

    program = "trnascan-plastid"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "tRNA gene detection [Plastid]"

    input_filepaths = [
            os.path.join( directories[("intermediate",  "1__partition")], "genomes.list"),
            os.path.join( directories["output"], "plastid", "*.fa"),
        ]

    output_filenames = [
        "*.tRNA",
    ]

    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        "search_mode":opts.trnascan_plastid_searchmode,
        "trnascan_options":opts.trnascan_plastid_options,

    }

    cmd = get_trnascan_cmd(**params)

    pipeline.add_step(
                id=program_label,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=False,
                validate_outputs=False,
                errors_ok=False,
                log_prefix=program_label,

    )

    # =============
    # Output
    # =============
    step = 11

    program = "symlink"
    program_label = "{}__{}".format(step, program)
    description = "Merging and symlinking results for output"

    # Add to directories
    output_directory = directories["output"]

    # i/o
    input_filepaths = [ 
        os.path.join( directories[("intermediate",  "1__partition")], "genomes.list"),
    ]

    output_filenames =  [
        # "*.faa",
        # "*.ffn",
        "*.gff",


    ]


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

   # =============
    # Output
    # =============
    step = 12

    program = "stats"
    program_label = "{}__{}".format(step, program)
    description = "Calculating statistics for sequences"

    # Add to directories
    output_directory = directories["output"]

    # i/o
    input_filepaths = [ 
        os.path.join(directories["output"], "*.fa"),
        # os.path.join(directories["output"], "mitochondrion"),
        # os.path.join(directories["output"], "plastid"),

    ]

    output_filenames =  [
        "genome_statistics.tsv",
        "gene_statistics.cds.tsv",
        "gene_statistics.rRNA.tsv",
        "gene_statistics.tRNA.tsv",


    ]


    output_filepaths = list(map(lambda fn:os.path.join(directories["output"], fn), output_filenames))

    
    params = {
    "input_filepaths":input_filepaths,
    "output_filepaths":output_filepaths,
    "output_directory":output_directory,
    "opts":opts,
    "directories":directories,
    }

    cmd = get_stats_cmd(**params)
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
    accessory_scripts = {
        "partition_gene_models.py",
        "compile_metaeuk_identifiers.py",
        "partition_organelle_sequences.py",
        "append_geneid_to_prodigal_gff.py",
        "append_geneid_to_barrnap_gff.py",
    }

    required_executables={
        "seqkit",
        "metaeuk",
        "pyrodigal",
        "barrnap",
        "tRNAscan-SE",
        
     } 
    if not opts.tiara_results:
        required_executables |= {"tiara"}
  

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
        executables[name] = "'{}'".format(os.path.join(opts.script_directory, name)) # Can handle spaces in path



    print(format_header( "Adding executables to path from the following source: {}".format(opts.path_config), "-"), file=sys.stdout)
    for name, executable in executables.items():
        if name in required_executables:
            print(name, executable, sep = " --> ", file=sys.stdout)
            os.environ[name] = executable.strip()
    print("", file=sys.stdout)


# Configure parameters
def configure_parameters(opts, directories):
    assert bool(opts.name) != bool(opts.scaffolds_to_bins), "--name and --scaffolds_to_bins are mutually exclusive.  Use --name if you are modeling genes on a single assembly and --scaffolds_to_bins in batch (faster for multiple assemblies)"
    if opts.name:
        opts.basename = opts.name
    else:
        opts.basename = "gene_models"

    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -f <scaffolds.fasta> -d <metaeuk_database> -i <scaffolds_to_bins.tsv>  -o <output_directory>".format(__program__)

    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-f","--fasta", type=str, required=True, help = "path/to/scaffolds.fasta")
    parser_io.add_argument("-t","--tiara_results", type=str, required=False, help = "path/to/scaffolds.fasta [Optional]")
    parser_io.add_argument("-n","--name", type=str, required=False, help = "path/to/scaffolds.fasta [Cannot be used with --scaffolds_to_bins]")
    parser_io.add_argument("-i","--scaffolds_to_bins", type=str, required=False,  help = "path/to/scaffolds_to_bins.tsv, [Optional] Format: [id_scaffold]<tab>[id_bin], No header. [Cannot be used with --name]")
    parser_io.add_argument("-o","--output_directory", type=str, default="eukaryotic_gene_modeling_output", help = "path/to/project_directory [Default: eukaryotic_gene_modeling_output]")
    parser_io.add_argument("-d", "--metaeuk_database", type=str,  required=True, help=f"MetaEuk/MMSEQS2 database (E.g., $VEBA_DATABASE/Classify/MicroEuk/MicroEuk50)")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    # parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Tiara
    parser_organelle = parser.add_argument_group('Organelle arguments')
    parser_organelle.add_argument("-u","--unknown_organelle_prediction", type=str, default="nuclear", help = "{unknown, nuclear} [Default: nuclear]")
    # parser_organelle.add_argument("--mitochondrion_suffix", type=str, default=".mtDNA", help = "Mitochondrion suffix [Default: .mtDNA]")
    # parser_organelle.add_argument("--plastid_suffix", type=str, default=".plastid", help = "Plastid suffix [Default: .plastid]")
    # parser_organelle.add_argument("--unknown_suffix", type=str, default=".unknown", help = "Unknown suffix [Default: .unknown]")
    parser_organelle.add_argument("--tiara_minimum_length", type=int, default=3000, help="Tiara | Minimum contig length. Anything lower than 3000 is not recommended. [Default: 3000]")
    parser_organelle.add_argument("--tiara_options", type=str, default="", help="Tiara | More options (e.g. --arg 1 ) [Default: '']")

    # MetaEuk
    parser_metaeuk = parser.add_argument_group('MetaEuk arguments')
    parser_metaeuk.add_argument("--metaeuk_sensitivity", type=float, default=4.0, help="MetaEuk | Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive  [Default: 4.0]")
    parser_metaeuk.add_argument("--metaeuk_evalue", type=float, default=0.01, help="MetaEuk | List matches below this E-value (range 0.0-inf) [Default: 0.01]")
    parser_metaeuk.add_argument("--metaeuk_split_memory_limit", type=str, default="36G", help="MetaEuk | Set max memory per split. E.g. 800B, 5K, 10M, 1G. Use 0 to use all available system memory. (Default value is experimental) [Default: 36G]")
    parser_metaeuk.add_argument("--metaeuk_options", type=str, default="", help="MetaEuk | More options (e.g. --arg 1 ) [Default: ''] https://github.com/soedinglab/metaeuk")

    # Pyrodigal
    parser_pyrodigal = parser.add_argument_group('Pyrodigal arguments (Mitochondria)')
    parser_pyrodigal.add_argument("--pyrodigal_minimum_gene_length", type=int, default=90, help="Pyrodigal | Minimum gene length [Default: 90]")
    parser_pyrodigal.add_argument("--pyrodigal_minimum_edge_gene_length", type=int, default=60, help="Pyrodigal | Minimum edge gene length [Default: 60]")
    parser_pyrodigal.add_argument("--pyrodigal_maximum_gene_overlap_length", type=int, default=60, help="Pyrodigal | Maximum gene overlap length [Default: 60]")
    parser_pyrodigal.add_argument("--pyrodigal_mitochondrial_genetic_code", type=int, default=4, help="Pyrodigal -g translation table (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi/) [Default: 4] (The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code))")
    parser_pyrodigal.add_argument("--pyrodigal_plastid_genetic_code", type=int, default=11, help="Pyrodigal -g translation table (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi/) [Default: 11] (The Bacterial, Archaeal and Plant Plastid Code))")

    # rRNA
    parser_barrnap = parser.add_argument_group('barrnap arguments')
    parser_barrnap.add_argument("--barrnap_length_cutoff", type=float, default=0.8,  help="barrnap | Proportional length threshold to label as partial [Default: 0.8]")
    parser_barrnap.add_argument("--barrnap_reject", type=float, default=0.25,  help="barrnap | Proportional length threshold to reject prediction [Default: 0.25]")
    parser_barrnap.add_argument("--barrnap_evalue", type=float, default=1e-6,  help="barrnap | Similarity e-value cut-off [Default: 1e-6]")

    # tRNA
    parser_trnascan = parser.add_argument_group('tRNAscan-SE arguments')
    parser_trnascan.add_argument("--trnascan_nuclear_options", type=str, default="", help="tRNAscan-SE | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/UCSC-LoweLab/tRNAscan-SE")
    parser_trnascan.add_argument("--trnascan_mitochondrial_searchmode", type=str, default="-O", help="tRNAscan-SE | Search mode [Default: '-O'] | Current best option according to developer: https://github.com/UCSC-LoweLab/tRNAscan-SE/issues/24")
    parser_trnascan.add_argument("--trnascan_mitochondrial_options", type=str, default="", help="tRNAscan-SE | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/UCSC-LoweLab/tRNAscan-SE")
    parser_trnascan.add_argument("--trnascan_plastid_searchmode", type=str, default="-O", help="tRNAscan-SE | Search mode [Default: '-O'] | https://github.com/UCSC-LoweLab/tRNAscan-SE")
    parser_trnascan.add_argument("--trnascan_plastid_options", type=str, default="", help="tRNAscan-SE | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/UCSC-LoweLab/tRNAscan-SE")


    # Options
    opts = parser.parse_args()

    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Threads
    if opts.n_jobs == -1:
        from multiprocessing import cpu_count 
        opts.n_jobs = cpu_count()
    assert opts.n_jobs >= 1, "--n_jobs must be ≥ 1.  To select all available threads, use -1."

    # Directories
    directories = dict()
    directories["project"] = create_directory(opts.output_directory)
    directories["output"] = create_directory(os.path.join(directories["project"], "output"))
    directories["log"] = create_directory(os.path.join(directories["project"], "log"))
    directories["tmp"] = create_directory(os.path.join(directories["project"], "tmp"))
    directories["checkpoints"] = create_directory(os.path.join(directories["project"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["project"], "intermediate"))
    os.environ["TMPDIR"] = directories["tmp"]

    # Temporary
    opts.mitochondrion_suffix = ""
    opts.plastid_suffix = ""
    opts.unknown_suffix = ""

    # Info
    print(format_header(__program__, "="), file=sys.stdout)
    print(format_header("Configuration:", "-"), file=sys.stdout)
    print("Python version:", sys.version.replace("\n"," "), file=sys.stdout)
    print("Python path:", sys.executable, file=sys.stdout) #sys.path[2]
    print("GenoPype version:", genopype_version, file=sys.stdout) #sys.path[2]
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

    # shutil.rmtree(directories["intermediate"])
    
   


if __name__ == "__main__":
    main(sys.argv[1:])
