#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, random
from collections import OrderedDict, defaultdict

import pandas as pd
import numpy as np

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.7.8"

# DATABASE_METAEUK="/usr/local/scratch/CORE/jespinoz/db/veba/v1.0/Classify/Eukaryotic/eukaryotic"


def get_preprocess_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # checkv end_to_end ${FASTA} ${OUT_DIR} -t ${N_JOBS} --restart

    cmd = [
        "(",
        "cat",
        opts.fasta,
        "|",
        os.environ["seqkit"],
        "seq",
        "-m {}".format(opts.minimum_contig_length),
        "|",
        os.environ["seqkit"],
        "grep",
        "--pattern-file {}".format(opts.contig_identifiers),
        ">",
        output_filepaths[0],
        ")",
    ]
    return cmd


def get_binning_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, prefix, seed):
    cmd = [
        # Copy fasta as unbinned.fasta in output in case no euks are found
        "cat",
        opts.fasta,
        "|",
        os.environ["seqkit"],
        "seq",
        "-m {}".format(opts.minimum_contig_length),
        ">",
        os.path.join(directories["output"], "unbinned.fasta"),

        "&&",
        # Binning wrapper
        os.environ["binning_wrapper.py"],
        "-f {}".format(opts.fasta), # scaffolds.fasta
        "-b {}".format(" ".join(opts.bam)), # mapped.sorted.bam
        "-o {}".format(output_directory),
        "-a {}".format(opts.algorithm),
        "-p {}".format(opts.n_jobs),
        "--minimum_contig_length {}".format(opts.minimum_contig_length), # mininimum contig length must be >= 1500 nts for MetaBat2
        "--minimum_genome_length {}".format(opts.minimum_genome_length),
        "--bin_prefix {}".format(prefix),
        "--random_state {}".format(seed),
        ]
    if opts.algorithm == "metabat2":
        if bool(opts.metabat2_options):
            cmd += [ "--metabat2_options {}".format(opts.metabat2_options)]
    if opts.algorithm == "concoct":
        cmd += [ 
            "--concoct_fragment_length {}".format(opts.concoct_fragment_length),
            "--concoct_overlap_length {}".format(opts.concoct_overlap_length),
            ]
        if bool(opts.concoct_options):
            cmd += ["--concoct_options {}".format(opts.concoct_options)]

    cmd += [ 
        # # Delete extra files
        # "&&",
        # "rm {} {} {} {}".format(
        #     os.path.join(output_directory, "bins.list"),
        #     os.path.join(output_directory, "binned.list"),
        #     os.path.join(output_directory, "unbinned.list"),
        #     os.path.join(output_directory, "genome_statistics.tsv"),
        # ),
        # Move scaffolds_to_bins.tsv to tmp
        "&&",

        "mv {} {}".format(
            os.path.join(output_directory, "scaffolds_to_bins.tsv"),
            directories["tmp"],
        ),

        # Non-Eukaryotic
        "&&",
        "cat",
        os.path.join(os.path.join(output_directory, "bins"), "*.fa"),
        "|",
        os.environ["seqkit"],
        "seq",
        "-m {}".format(opts.tiara_minimum_length),
        ">" ,
        os.path.join(directories["tmp"], "scaffolds.binned.gte{}.fasta".format(opts.tiara_minimum_length)),


        # Tiara
        "&&",
        "mkdir -p {}".format(os.path.join(output_directory, "consensus_domain_classification")),
        "&&",
        "mkdir -p {}".format(os.path.join(output_directory,  "bins", "non-eukaryota")),

        "&&",
        os.environ["tiara"],
        "-i {}".format(os.path.join(directories["tmp"], "scaffolds.binned.gte{}.fasta".format(opts.tiara_minimum_length))),
        "-o {}".format(os.path.join(output_directory, "consensus_domain_classification", "tiara_output.tsv")),
        "--probabilities",
        "-m {}".format(opts.tiara_minimum_length),
        "-t {}".format(opts.n_jobs),
        opts.tiara_options,

        # Predict domain
        "&&",
        os.environ["consensus_domain_classification.py"],
        "-i {}".format(os.path.join(directories["tmp"], "scaffolds_to_bins.tsv")),
        "-t {}".format(os.path.join(output_directory, "consensus_domain_classification", "tiara_output.tsv")),
        "-o {}".format(os.path.join(output_directory, "consensus_domain_classification")),
        "--logit_transform {}".format(opts.logit_transform),

        # Move all (this is a hack)
        "&&",


        "mv {} {}".format(os.path.join(output_directory, "bins", "*.fa"), os.path.join(output_directory,  "bins", "non-eukaryota")),
        # "for FP in %s; do BN=$(basename ${FP} .fa); mv $FP %s/%s${BN}.fa; done"%(
        #     os.path.join(output_directory, "bins", "*.fa"),
        #     os.path.join(output_directory,  "bins", "non-eukaryota"),
        #     prefix,
        #     ),
        # Move just Eukaryota back
        "&&",
        "for ID_GENOME in $(cat %s); do mv %s/${ID_GENOME}.* %s; done"%(
            os.path.join(output_directory, "consensus_domain_classification", "eukaryota.list"),
            os.path.join(output_directory,  "bins", "non-eukaryota"),
            os.path.join(output_directory,  "bins"),
            ),
        # Non-Eukaryotic scaffolds
        # "&&",
        # "> {}".format(os.path.join(directories["tmp"], "non-eukaryota.scaffolds.fasta")), # Create empty file
        "&&",
        "(for FP in {}; do cat $FP >> {}; done) 2> /dev/null || > {}".format( # Handle edge cases where there aren't any non-eukaryota bins
            os.path.join(output_directory, "bins", "non-eukaryota", "*.fa"),
            os.path.join(directories["tmp"], "non-eukaryota.scaffolds.fasta"),
            os.path.join(directories["tmp"], "non-eukaryota.scaffolds.fasta"),
            ), 
        "&&",
        "cat",
        os.path.join(directories["tmp"], "non-eukaryota.scaffolds.fasta"),
        "|",
        "grep",
        '"^>"',
        "|",
        "cut -c2-",
        ">",
        os.path.join(output_directory, "non-eukaryota.scaffolds.list"),

        # Remove non-eukaryotic scaffolds
        "&&",
        "cat",
        os.path.join(os.path.join(directories["tmp"], "scaffolds_to_bins.tsv")),
        "|",
        "grep",
        "-v",
        "-f {}".format(os.path.join(output_directory,  "non-eukaryota.scaffolds.list")),
        ">",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),

        # EXPERIMENTAL
        "&&",
        "[ -s {} ] || (echo 'No eukaryotic bins' && exit 1)".format(os.path.join(output_directory, "scaffolds_to_bins.tsv")),

        # Save non-eukaryotic scaffolds separately
        "&&",
        "cat",
        os.path.join(os.path.join(directories["tmp"], "scaffolds_to_bins.tsv")),
        "|",
        "grep",
        "-f {}".format(os.path.join(output_directory, "non-eukaryota.scaffolds.list")),
        ">",
        os.path.join(output_directory,  "non-eukaryota.scaffolds_to_bins.tsv"),
        "||",
        # "true",
        ">",
        os.path.join(output_directory,  "non-eukaryota.scaffolds_to_bins.tsv"), # Make empty file
        "&&",
        "cut",
        "-f1",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),
        ">",
        os.path.join(output_directory, "binned.list"),

        "&&",
        # Unique bins
        "cut -f2",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),
        "|",
        "sort",
        "-u",
        ">",
        os.path.join(output_directory, "bins.list"),
        "&&",

        "(",
        "rm -rf {} {} {} {} {}".format(
            os.path.join(output_directory,"bins"), # Remove because the info is in the scaffolds to bins file
            # os.path.join(output_directory,  "bins", "non-eukaryota", "*.fa"), # Remove because the info is in the scaffolds to bins file
            os.path.join(directories["tmp"], "non-eukaryota.scaffolds.fasta"),
            os.path.join(os.path.join(directories["tmp"], "scaffolds_to_bins.tsv")),
            os.path.join(directories["tmp"], "scaffolds.binned.gte{}.fasta".format(opts.tiara_minimum_length)),
            os.path.join(output_directory, "itermediate"), 
        ),
        ")",
    ]

    return cmd



def get_metaeuk_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [
        # Get the eukaryotic contigs
        "cat",
        opts.fasta,
        "|",
        os.environ["seqkit"],
        "seq",
        "-m {}".format(opts.minimum_contig_length),
        "|",
        os.environ["seqkit"],
        "grep",
        "--pattern-file {}".format(os.path.join(directories[("intermediate",  "1__binning_{}".format(opts.algorithm))], "binned.list")),
        ">",
        os.path.join(directories["tmp"], "scaffolds.binned.eukaryotic.fasta"), # contigs

        # Run MetaEuk
        "&&",
        os.environ["metaeuk"],
        "easy-predict",
        "--threads {}".format(opts.n_jobs),
        "-s {}".format(opts.metaeuk_sensitivity),
        "-e {}".format(opts.metaeuk_evalue),
        opts.metaeuk_options,
        os.path.join(directories["tmp"], "scaffolds.binned.eukaryotic.fasta"), # contigs
        opts.metaeuk_database, # db
        os.path.join(output_directory, "metaeuk"), # output prefix
        os.path.join(directories["tmp"],"metaeuk"),

        # Convert MetaEuk identifiers
        "&&",
        os.environ["compile_metaeuk_identifiers.py"],
        "--cds {}".format(os.path.join(output_directory, "metaeuk.codon.fas")),
        "--protein {}".format(os.path.join(output_directory, "metaeuk.fas")),
        "-o {}".format(output_directory),
        "-b gene_models",

        # Partition the gene models and genomes
        "&&",
        os.environ["partition_gene_models.py"],
        "-i {}".format(input_filepaths[1]),
        "-f {}".format(os.path.join(directories["tmp"], "scaffolds.binned.eukaryotic.fasta")),
        "-g {}".format(os.path.join(output_directory, "gene_models.gff")),
        "-d {}".format(os.path.join(output_directory, "gene_models.ffn")),
        "-a {}".format(os.path.join(output_directory, "gene_models.faa")),
        "-o {}".format(os.path.join(output_directory, "genomes")),

        # Remove temporary files
        "&&",
        "rm -rf {} {} {} {} {}".format(
            os.path.join(output_directory, "*.fas"), # output prefix
            os.path.join(output_directory, "metaeuk.gff"), # output prefix
            os.path.join(output_directory, "gene_models.*"), # output prefix
            os.path.join(directories["tmp"], "scaffolds.binned.eukaryotic.fasta"),
            os.path.join(directories["tmp"],"metaeuk", "*"),
        )
    ]
    return cmd

# def get_busco_offline_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

def get_busco_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [
"""
# Create busco output directory
rm -rf %s
mkdir -p %s

# Iterate through protein fasta from MetaEuk
for FP in %s;

    # Get MAG identifier from protein fasta filename
    do ID_GENOME=$(basename $FP .faa);

    # Create MAG-specific subdirectory within busco output
    OUT_DIR=%s/${ID_GENOME}
    mkdir -p ${OUT_DIR}

    echo $FP
    echo $ID_GENOME
    echo $OUT_DIR

    # BUSCO Command
    %s --force -i $FP -o $OUT_DIR -m protein --auto-lineage-euk -c %d --evalue %f --download_path %s

    # Remove big intermediate files

    echo "Removing run_*. .. ... ..... ........"
    rm -rf ${OUT_DIR}/run_*

    echo "Removing auto_lineage. .. ... ..... ........"
    rm -rf ${OUT_DIR}/auto_lineage

# End for-loop
done

# Removing temporary busco files
rm -rf %s
"""%( 
    # Args
    os.path.join(output_directory, "busco_output"),
    os.path.join(output_directory, "busco_output"),
    os.path.join(directories[("intermediate",  "2__metaeuk")], "genomes", "*.faa"),
    os.path.join(output_directory, "busco_output"),
    os.environ["busco"],
    opts.n_jobs,
    opts.busco_evalue,
    os.path.join(directories["tmp"],"busco"),
    os.path.join(directories["tmp"],"busco", "*"),
    ),
    os.environ["merge_busco_json.py"],
    "-i {}".format(os.path.join(output_directory, "busco_output")),
    "-j {}".format(os.path.join(output_directory, "busco_results.json")),
    "-o {}".format(os.path.join(output_directory, "busco_results.tsv")),
    "&&",
    os.environ["filter_busco_results.py"],
    "-i {}".format(os.path.join(output_directory, "busco_results.tsv")),
    "-g {}".format(directories[("intermediate",  "2__metaeuk")]),
    "-o {}".format(os.path.join(output_directory, "filtered")),
    "-m {}".format(opts.minimum_contig_length),
    "-f {}".format(opts.fasta),
    "--completeness {}".format(opts.busco_completeness),
    "--contamination {}".format(opts.busco_contamination),
    "--unbinned",
    ]


    
    return cmd

def get_featurecounts_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # ORF-Level Counts
    cmd = [

        "cat", 
        os.path.join(directories[("intermediate",  "2__metaeuk")], "genomes", "*.gff"),
        ">",
        os.path.join(directories["tmp"], "gene_models.eukaryotic.gff"),
        "&&",
        "mkdir -p {}".format(os.path.join(directories["tmp"], "featurecounts")),
        "&&",

    "(",
        os.environ["featureCounts"],
        "-G {}".format(opts.fasta),
        "-a {}".format(os.path.join(directories["tmp"], "gene_models.eukaryotic.gff")),
        "-o {}".format(os.path.join(output_directory, "featurecounts.orfs.tsv")),
        "-F GTF",
        "--tmpDir {}".format(os.path.join(directories["tmp"], "featurecounts")),
        "-T {}".format(opts.n_jobs),
        "-g gene_id",
        "-t CDS",
        opts.featurecounts_options,
        " ".join(opts.bam),
    ")",
    "&&",
    "gzip -f {}".format(os.path.join(output_directory, "featurecounts.orfs.tsv")),
    "&&",
    "rm -rf {}".format(os.path.join(directories["tmp"], "gene_models.eukaryotic.gff")),
    ]
    return cmd

def get_output_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [
        "rm -rf {}".format(os.path.join(output_directory, "*")),
    ]


    for fp in input_filepaths:
        fn = fp.split("/")[-1]
        cmd += [ 
            "&&",
            "ln -sf",
            os.path.realpath(os.path.join(fp)),
            os.path.join(output_directory,fn),
        ]
        
        
    cmd += [ 
        "&&",
        # Statistics
        "(",
        os.environ["seqkit"],
        "stats",
        "-a",
        "-b",
        "-T",
        "-j {}".format(opts.n_jobs),
        os.path.join(output_directory, "genomes", "*.fa"),
        ">",
        os.path.join(output_directory,"genome_statistics.tsv"),
        ")",
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
    pipeline = ExecutablePipeline(name=__program__, description=opts.name, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])
    iteration = 1
    steps = dict()

    seed = opts.random_state
    if opts.random_state == -1:
        seed = random.choice(range(1,10000))
    if opts.random_state == 0:
        seed = iteration

    if bool(opts.contig_identifiers):
        # ==========
        # Subset
        # ==========
        step = 0

        program = "preprocessing"

        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories["preprocessing"] = create_directory(os.path.join(directories["project"], "preprocessing"))


        # Info
        description = "Splitting contigs from fasta"
        
        # i/o
        input_filepaths = [
            opts.fasta,
            opts.contig_identifiers,
        ]

        output_filenames = [
            "input.fasta",
        ]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_preprocess_cmd(**params)
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
        )

        opts.fasta = os.path.join(directories["preprocessing"], "input.fasta")



    # ==========
    # Binning
    # ==========
    step = 1

    program = "binning_{}".format(opts.algorithm)

    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Binning via {}".format(opts.algorithm)
    
    # i/o
    output_filepaths = [
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),
    ]


    input_filepaths = [opts.fasta, *opts.bam]

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        "prefix":"{}__{}__{}.{}__".format(opts.name, opts.algorithm.upper(), "E", iteration),
        "seed":seed,
    }

    cmd = get_binning_cmd(**params)
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
    )


    steps[program] = step

    # =============
    # MetaEuk
    # =============
    step = 2

    program = "metaeuk"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Ab initio eukaryotic gene prediction"


    input_filepaths = [
            # os.path.join(directories[("intermediate",  "1__binning_{}".format(opts.algorithm))], "bins"),
            opts.fasta,
            os.path.join(directories[("intermediate",  "1__binning_{}".format(opts.algorithm))], "scaffolds_to_bins.tsv"),
        ]

    output_filenames = [
        "genomes/*.fa",
        "genomes/*.faa",
        "genomes/*.gff",
        "genomes/*.ffn",
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

    steps[program] = step

    # =============
    # BUSCO
    # =============
    step = 3

    program = "busco"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Lineage-specific quality assessment"


    input_filepaths = [
            os.path.join(directories[("intermediate",  "2__metaeuk")], "genomes", "*.faa"),
        ]

    output_filenames = [
        # "busco_output/*/*.json",
        # "busco_results.json",
        # "busco_results.tsv",
        "filtered/busco_results.filtered.tsv",
        "filtered/genomes/*.fa",


    ]

    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_busco_cmd(**params)

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

    steps[program] = step

    # ==========
    # featureCounts
    # ==========
    step = 4

    # Info
    program = "featurecounts"
    program_label = "{}__{}".format(step, program)
    description = "Counting reads"

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # i/o
    input_filepaths = [ 
        opts.fasta,
        os.path.join(directories[("intermediate",  "2__metaeuk")], "genomes", "*.gff"),
        *opts.bam,
    ]

    output_filenames = ["featurecounts.orfs.tsv.gz"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_featurecounts_cmd(**params)
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
                # acceptable_returncodes={0,1,2},                    
                log_prefix=program_label,
    )
    steps[program] = step

    # =============
    # Output
    # =============
    step = 5

    program = "output"
    program_label = "{}__{}".format(step, program)
    description = "Merging results for output"

    # Add to directories
    output_directory = directories["output"]

    # i/o
    input_filepaths = [ 

        # BUSCO
        os.path.join(directories[("intermediate",  "3__busco")],  "filtered", "identifier_mapping.metaeuk.tsv"),
        os.path.join(directories[("intermediate",  "3__busco")],  "filtered", "busco_results.filtered.tsv"),
        os.path.join(directories[("intermediate",  "3__busco")], "filtered", "scaffolds_to_bins.tsv"),
        os.path.join(directories[("intermediate",  "3__busco")],  "filtered", "bins.list"),
        os.path.join(directories[("intermediate",  "3__busco")],  "filtered", "binned.list"),
        os.path.join(directories[("intermediate",  "3__busco")],  "filtered", "unbinned.list"),
        os.path.join(directories[("intermediate",  "3__busco")],  "filtered", "unbinned.fasta"),
        os.path.join(directories[("intermediate",  "3__busco")],  "filtered", "genomes"),

        # featureCounts
        os.path.join(directories[("intermediate",  "4__featurecounts")], "featurecounts.orfs.tsv.gz"),
    ]

    output_filenames =  [

        "scaffolds_to_bins.tsv", 
        "binned.list",
        "unbinned.fasta",
        "unbinned.list",
        "genomes",
        "bins.list",
        "genome_statistics.tsv",
        "featurecounts.orfs.tsv.gz",
        "identifier_mapping.metaeuk.tsv",
        "busco_results.filtered.tsv",
    ]


    output_filepaths = list(map(lambda fn:os.path.join(directories["output"], fn), output_filenames))

    
    params = {
    "input_filepaths":input_filepaths,
    "output_filepaths":output_filepaths,
    "output_directory":output_directory,
    "opts":opts,
    "directories":directories,
    }

    cmd = get_output_cmd(**params)
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
    steps[program] = step

    return pipeline



# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """
    accessory_scripts = {
        "binning_wrapper.py",
        "scaffolds_to_bins.py",
        "partition_gene_models.py",
        "compile_metaeuk_identifiers.py",
        "append_geneid_to_prodigal_gff.py",
        "consensus_domain_classification.py",
        "merge_busco_json.py",
        "filter_busco_results.py",
    }

    required_executables={
                "seqkit",
                "metaeuk",
                "busco",
                "featureCounts",
                "tiara",


     } | accessory_scripts

    if opts.algorithm == "metabat2":
        required_executables |= {
                "metabat2",
                "coverm",
        }

    if opts.algorithm == "concoct":
        required_executables |= {
                "concoct",
                "cut_up_fasta.py", 
                "concoct_coverage_table.py",
                "merge_cutup_clustering.py", 
                "extract_fasta_bins.py", 
        }

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
            executables[name] = "python " + os.path.join(opts.script_directory, "scripts", name)
        else: 
            executables[name] = os.path.join(opts.script_directory, "scripts", name)


    print(format_header( "Adding executables to path from the following source: {}".format(opts.path_config), "-"), file=sys.stdout)
    for name, executable in executables.items():
        if name in required_executables:
            print(name, executable, sep = " --> ", file=sys.stdout)
            os.environ[name] = executable.strip()
    print("", file=sys.stdout)


# Configure parameters
def configure_parameters(opts, directories):
    # assert opts.reference_assembly is not None, "Must include --reference_assembly"
    assert_acceptable_arguments(opts.algorithm, {"concoct", "metabat2"})

    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -f <scaffolds.fasta> -b <mapped.sorted.bam> -n <name> -o <output_directory>".format(__program__)

    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-f","--fasta", type=str, required=True, help = "path/to/scaffolds.fasta")
    parser_io.add_argument("-b","--bam", type=str, nargs="+", required=True, help = "path/to/mapped.sorted.bam files separated by spaces. ")
    parser_io.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    parser_io.add_argument("-l","--contig_identifiers", type=str,  help = "path/to/contigs.list [Optional]")
    parser_io.add_argument("-o","--project_directory", type=str, default="veba_output/binning/eukaryotic", help = "path/to/project_directory [Default: veba_output/binning/eukaryotic]")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')

    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
    # parser_utility.add_argument("--remove_temporary_fasta", action="store_true", help="If contig identifiers were provided and a fasta is generated, remove this file")

    # parser_utility.add_argument("-c", "--CONDA_PREFIX", type=str, default=None, help = "Set a conda environment")

    # Databases
    parser_databases = parser.add_argument_group('Database arguments')
    parser_databases.add_argument("--veba_database", type=str,  help=f"VEBA database location.  [Default: $VEBA_DATABASE environment variable]")

    # Binning
    parser_binning = parser.add_argument_group('Binning arguments')
    parser_binning.add_argument("-a", "--algorithm", type=str, default="metabat2", help="Binning algorithm: {concoct, metabat2}  [Default: metabat2] ")
    parser_binning.add_argument("-m", "--minimum_contig_length", type=int, default=1500, help="Minimum contig length.  [Default: 1500] ")
    parser_binning.add_argument("-s", "--minimum_genome_length", type=int, default=2000000, help="Minimum genome length.  [Default: 2000000] ")

    parser_binning.add_argument("--concoct_fragment_length", type=int, default=10000, help="CONCOCT | Fragment length [Default: 10000] ")
    parser_binning.add_argument("--concoct_overlap_length", type=int, default=0, help="CONCOCT | Fragment overlap length [Default: 0] ")
    parser_binning.add_argument("--concoct_options", type=str, default="", help="CONCOCT | More options (e.g. --arg 1 ) [Default: '']")
    parser_binning.add_argument("--metabat2_options", type=str, default="", help="MetaBat2 | More options (e.g. --arg 1 ) [Default: ''] | https://bitbucket.org/berkeleylab/metabat/src/master/")


    # Tiara
    parser_domain = parser.add_argument_group('Domain classification arguments')
    parser_domain.add_argument("--logit_transform", type=str, default="softmax", help = " Transformation for consensus_domain_classification: {softmax, tss} [Default: softmax")
    parser_domain.add_argument("--tiara_minimum_length", type=int, default=3000, help="Tiara | Minimum contig length. Anything lower than 3000 is not recommended. [Default: 3000]")
    parser_domain.add_argument("--tiara_options", type=str, default="", help="Tiara | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/ibe-uw/tiara")

    # MetaEuk
    parser_metaeuk = parser.add_argument_group('MetaEuk arguments')
    parser_metaeuk.add_argument("--metaeuk_sensitivity", type=float, default=4.0, help="MetaEuk | Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive  [Default: 4.0]")
    parser_metaeuk.add_argument("--metaeuk_evalue", type=float, default=0.01, help="MetaEuk | List matches below this E-value (range 0.0-inf) [Default: 0.01]")
    # parser_metaeuk.add_argument("--metaeuk_database", type=str, default=DATABASE_METAEUK, help="MetaEuk | More options (e.g. --arg 1 ) [Default: {}]".format(DATABASE_METAEUK))
    parser_metaeuk.add_argument("--metaeuk_options", type=str, default="", help="MetaEuk | More options (e.g. --arg 1 ) [Default: ''] https://github.com/soedinglab/metaeuk")

    # BUSCO
    parser_busco = parser.add_argument_group('BUSCO arguments')
    # parser_busco.add_argument("--busco_offline", type=str, help="BUSCO | Offline database path")
    parser_busco.add_argument("--busco_completeness", type=float, default=50.0, help = "BUSCO completeness [Default: 50.0]")
    parser_busco.add_argument("--busco_contamination", type=float, default=10.0, help = "BUSCO contamination [Default: 10.0]")

    parser_busco.add_argument("--busco_evalue", type=float, default=0.001, help="BUSCO | E-value cutoff for BLAST searches. Allowed formats, 0.001 or 1e-03 [Default: 1e-03]")

    # featureCounts
    parser_featurecounts = parser.add_argument_group('featureCounts arguments')
    parser_featurecounts.add_argument("--featurecounts_options", type=str, default="", help="featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/")


    # Options
    opts = parser.parse_args()

    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Database
    if opts.veba_database is None:
        assert "VEBA_DATABASE" in os.environ, "Please set the following environment variable 'export VEBA_DATABASE=/path/to/veba_database' or provide path to --veba_database"
        opts.veba_database = os.environ["VEBA_DATABASE"]
    opts.metaeuk_database = os.path.join(opts.veba_database, "Classify", "Microeukaryotic", "microeukaryotic")


    # Directories
    directories = dict()
    directories["project"] = create_directory(opts.project_directory)
    directories["sample"] = create_directory(os.path.join(directories["project"], opts.name))
    directories["preprocessing"] = create_directory(os.path.join(directories["sample"], "preprocessing"))
    directories["output"] = create_directory(os.path.join(directories["sample"], "output"))
    directories["log"] = create_directory(os.path.join(directories["sample"], "log"))
    directories["tmp"] = create_directory(os.path.join(directories["sample"], "tmp"))
    directories["checkpoints"] = create_directory(os.path.join(directories["sample"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["sample"], "intermediate"))
    os.environ["TMPDIR"] = directories["tmp"]

    # Info
    print(format_header(__program__, "="), file=sys.stdout)
    print(format_header("Configuration:", "-"), file=sys.stdout)
    print(format_header("Name: {}".format(opts.name), "."), file=sys.stdout)
    print("Python version:", sys.version.replace("\n"," "), file=sys.stdout)
    print("Python path:", sys.executable, file=sys.stdout) #sys.path[2]
    print("Script version:", __version__, file=sys.stdout)
    print("VEBA Database:", opts.veba_database, file=sys.stdout)    
    print("Moment:", get_timestamp(), file=sys.stdout)
    print("Directory:", os.getcwd(), file=sys.stdout)
    print("Commands:", list(filter(bool,sys.argv)),  sep="\n", file=sys.stdout)
    configure_parameters(opts, directories)
    sys.stdout.flush()


    # Run pipeline
    with open(os.path.join(directories["sample"], "commands.sh"), "w") as f_cmds:
        pipeline = create_pipeline(
                    opts=opts,
                    directories=directories,
                    f_cmds=f_cmds,
        )
        pipeline.compile()
        pipeline.execute(restart_from_checkpoint=opts.restart_from_checkpoint)
   


if __name__ == "__main__":
    main()
