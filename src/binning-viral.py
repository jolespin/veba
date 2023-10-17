#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
from collections import OrderedDict, defaultdict

import pandas as pd
import numpy as np

# Soothsayer Ecosystem
from genopype import *
from genopype import __version__ as genopype_version
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.10.16"

# geNomad
def get_genomad_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [
        "cat",
        opts.fasta,
        "|",
        os.environ["seqkit"],
        "seq",
        "-m {}".format(opts.minimum_contig_length),
        "|",
        # seqkit replace -p "\s.+"
        os.environ["seqkit"],
        "replace",
        '-p "\s.+"',
        ">",
        os.path.join(directories["tmp"], "input.fasta"),

            "&&",

        "cp",
        os.path.join(directories["tmp"], "input.fasta"),
        os.path.join(directories["output"], "unbinned.fasta"),

            "&&",

        os.environ["genomad"],
        "end-to-end",
        "--cleanup",
        "--threads {}".format(opts.n_jobs),
        "--verbose",
        "--enable-score-calibration",
        "--disable-find-proviruses" if not opts.include_provirus_detection else "",
        "--sensitivity {}".format(opts.sensitivity),
        "--splits {}".format(opts.splits),
        "--composition {}".format(opts.composition),
        "--min-score {}".format(opts.minimum_score),
        "--max-fdr {}".format(opts.genomad_qvalue),
        "--min-plasmid-marker-enrichment {}".format(opts.minimum_plasmid_marker_enrichment),
        "--min-virus-marker-enrichment {}".format(opts.minimum_virus_marker_enrichment),
        "--min-plasmid-hallmarks {}".format(opts.minimum_plasmid_hallmarks),
        "--min-virus-hallmarks {}".format(opts.minimum_virus_hallmarks),
        "--max-uscg {}".format(opts.maximum_universal_single_copy_genes),
        opts.genomad_options,
        os.path.join(directories["tmp"], "input.fasta"), # Input
        os.path.join(output_directory, "intermediate"), # Output
        os.path.join(opts.veba_database, "Classify", "geNomad"),

           "&&",

        "cp -rf", 
        os.path.join(output_directory, "intermediate", "input_summary", "input_virus_summary.tsv"), 
        os.path.join(output_directory, "virus_summary.tsv"), 

            "&&",

        "cp -rf", 
        os.path.join(output_directory, "intermediate", "input_summary", "input_plasmid_summary.tsv"), 
        os.path.join(output_directory, "plasmid_summary.tsv"), 

            "&&",

        "cp -rf", 
        os.path.join(output_directory, "intermediate", "input_annotate", "input_taxonomy.tsv"), 
        os.path.join(output_directory, "viral_taxonomy.tsv"), 

            "&&",

        "cat",
        os.path.join(output_directory, "virus_summary.tsv"), 
        "|",
        "tail -n +2",
        "|",
        "cut -f1",
        ">",
        os.path.join(output_directory, "virus_scaffolds.list"), 

            "&&",
            
        "cat",
        os.path.join(output_directory, "plasmid_summary.tsv"), 
        "|",
        "tail -n +2",
        "|",
        "cut -f1",
        ">",
        os.path.join(output_directory, "plasmid_scaffolds.list"), 

    ]

    return cmd

# VirFinder
def get_virfinder_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [

        "cat",
        opts.fasta,
        "|",
        os.environ["seqkit"],
        "seq",
        "-m {}".format(opts.minimum_contig_length),
        "|",
        os.environ["seqkit"],
        "replace",
        '-p "\s.+"',
        ">",
        os.path.join(directories["tmp"], "input.fasta"),

            "&&",
            
        "cp",
        os.path.join(directories["tmp"], "input.fasta"),
        os.path.join(directories["output"], "unbinned.fasta"),

            "&&",

        os.environ["virfinder_wrapper.r"],
        "-f {}".format(os.path.join(directories["tmp"], "input.fasta")),
        "-o {}".format(output_filepaths[0]),
        "-q",
        opts.virfinder_options,

            "&&",

        "cat",
        output_filepaths[0],
        "|",
        "python",
        "-c",
        r""" "import sys, pandas as pd; print(*pd.read_csv(sys.stdin, sep='\t', index_col=0).query('{} < {}').index, sep='\n')" """.format("qvalue" if opts.use_qvalue else "pvalue", opts.virfinder_pvalue),
        ">",
        output_filepaths[1],

            "&&",

        "cat",
        opts.fasta,
        "|",
        os.environ["seqkit"],
        "grep",
        "--pattern-file {}".format(output_filepaths[1]),
        ">",
        os.path.join(directories["tmp"], "input.fasta"),

            "&&",

    #     os.environ["genomad"],
    #     "annotate",
    #     "--cleanup",
    #     "--threads {}".format(opts.n_jobs),
    #     "--splits {}".format(opts.splits),
    #     "--evalue {}".format(opts.mmseqs2_evalue),
    #     "--sensitivity {}".format(opts.sensitivity),
    # ]

    # if opts.use_minimal_database_for_taxonomy:
    #     cmd += [
    #         "--use-minimal-db",
    #     ]

    # cmd += [ 
    #     os.path.join(directories["tmp"], "input.fasta"), # Input
    #     os.path.join(output_directory, "intermediate"), # Output
    #     os.path.join(opts.veba_database, "Classify", "geNomad"),

    #         "&&",

    #     "cp -rf", 
    #     os.path.join(output_directory, "intermediate", "input_annotate", "input_taxonomy.tsv"), 
    #     os.path.join(output_directory, "viral_taxonomy.tsv"), 

    os.environ["genomad_taxonomy_wrapper.py"],
    "-f {}".format( os.path.join(directories["tmp"], "input.fasta")),
    "-o {}".format(output_directory),
    "--veba_database {}".format(opts.veba_database),
    "--sensitivity {}".format(opts.sensitivity),
    "--splits {}".format(opts.splits),
    "--mmseqs2_evalue {}".format(opts.mmseqs2_evalue),
    "--n_jobs {}".format(opts.n_jobs),
    "--unclassified_lineage_label Unclassified",
    "--unclassified_taxid -1",
    "--unclassified_score 1.0",
    "--use_minimal_database_for_taxonomy" if opts.use_minimal_database_for_taxonomy else "",
    ]


    return cmd

def get_checkv_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # checkv end_to_end ${FASTA} ${OUT_DIR} -t ${N_JOBS} --restart

    cmd = [
        "cat",
        opts.fasta,
        "|",
        # os.environ["seqkit"],
        # "seq",
        # "-m {}".format(opts.minimum_contig_length),
        # "|",
        os.environ["seqkit"],
        "grep",
        "--pattern-file {}".format(input_filepaths[1]),
        ">",
        os.path.join(directories["tmp"], "input.fasta"),

            "&&",

        "rm -rf {}".format(output_directory),

            "&&",

        os.environ["checkv"],
        "end_to_end",
        os.path.join(directories["tmp"], "input.fasta"),
        output_directory,
        "-t {}".format(opts.n_jobs),
        "-d {}".format(opts.checkv_database),
        "--restart",
        opts.checkv_options, 

            "&&",

        os.environ["filter_checkv_results.py"],
        "-i {}".format(os.path.join(output_directory, "quality_summary.tsv")),
        "-c {}".format(os.path.join(output_directory, "completeness.tsv")),
        "-f {}".format(opts.fasta),
        "-o {}".format(os.path.join(output_directory, "filtered")),
        "-m {}".format(opts.minimum_contig_length),
        "--unbinned",
        "-p {}__{}__Virus.".format(opts.name, opts.algorithm.upper()),
        "--multiplier_viral_to_host_genes {}".format(opts.multiplier_viral_to_host_genes),
        "--completeness {}".format(opts.checkv_completeness),
        "--checkv_quality {}".format(opts.checkv_quality),
        "--miuvig_quality {}".format(opts.miuvig_quality),
        "--genomad_virus_taxonomy {}".format(input_filepaths[2]),

    ] 

    if opts.include_provirus_detection:
        cmd += ["--include_provirus_detection"]

    if opts.algorithm == "genomad":
        cmd += [ 
        "--genomad_virus_summary {}".format(input_filepaths[3]),
        ]
    
    return cmd

# Prodigal
def get_prodigal_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [
        "cat",
        os.path.join(input_filepaths[0], "*.fa"),
        "|",
        os.environ["prodigal-gv"],
        "-p meta",
        "-g {}".format(opts.prodigal_genetic_code),
        "-f gff",
        "-d {}".format(os.path.join(output_directory, "gene_models.ffn")),
        "-a {}".format(os.path.join(output_directory, "gene_models.faa")),
        "|",
        os.environ["append_geneid_to_prodigal_gff.py"],
        "-a gene_id",
        ">",
        os.path.join(output_directory, "gene_models.gff"),

            "&&",
            
        os.environ["partition_gene_models.py"],
        "-i {}".format(input_filepaths[1]),
        "-g {}".format(os.path.join(output_directory, "gene_models.gff")),
        "-d {}".format(os.path.join(output_directory, "gene_models.ffn")),
        "-a {}".format(os.path.join(output_directory, "gene_models.faa")),
        "-o {}".format(os.path.join(directories[("intermediate",  "2__checkv")],"filtered", "genomes")),

"""

OUTPUT_DIRECTORY={}

for GENOME_FASTA in {};
do
    ID=$(basename $GENOME_FASTA .fa)
    DIR_GENOME=$(dirname $GENOME_FASTA)
    GFF_CDS=$DIR_GENOME/$ID.gff
    GFF_OUTPUT=$OUTPUT_DIRECTORY/$ID.gff
    >$GFF_OUTPUT.tmp
    {} -f $GENOME_FASTA -o $GFF_OUTPUT.tmp -n $ID -c $GFF_CDS -d Virus
    mv $GFF_OUTPUT.tmp $GFF_OUTPUT
done

""".format(
        os.path.join(directories[("intermediate",  "2__checkv")],"filtered", "genomes"),
        os.path.join(directories[("intermediate",  "2__checkv")],"filtered", "genomes", "*.fa"),
        os.environ["compile_gff.py"],
    ),


        "rm -rf",
        os.path.join(output_directory, "gene_models.gff"),
        os.path.join(output_directory, "gene_models.ffn"),
        os.path.join(output_directory, "gene_models.faa"),

    ]
    return cmd

def get_featurecounts_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # ORF-Level Counts
    cmd = [
    "mkdir -p {}".format(os.path.join(directories["tmp"], "featurecounts")),
    "&&",
    "cat",
    os.path.join(directories[("intermediate",  "2__checkv")],"filtered", "genomes","*.gff"),
    ">",
    os.path.join(directories["tmp"], "featurecounts", "gene_models.gff"),
    "&&",
        os.environ["featureCounts"],
        "-G {}".format(input_filepaths[0]),
        "-a {}".format(os.path.join(directories["tmp"], "featurecounts", "gene_models.gff")),
        "-o {}".format(os.path.join(output_directory, "featurecounts.orfs.tsv")),
        "-F GTF",
        "--tmpDir {}".format(os.path.join(directories["tmp"], "featurecounts")),
        "-T {}".format(opts.n_jobs),
        "-g gene_id",
        "-t CDS",
        "-L" if opts.long_reads else "-p --countReadPairs",
        opts.featurecounts_options,
        " ".join(opts.bam),
    "&&",
    "gzip -f {}".format(os.path.join(output_directory, "featurecounts.orfs.tsv")),
    "&&",
    "rm -rf {}".format(os.path.join(directories["tmp"], "featurecounts","*" )),

    ]
    return cmd

def get_output_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # checkm lineage_wf --tab_table -f checkm_output/${ID}/output.tab --pplacer_threads ${N_JOBS} -t ${N_JOBS} -x fa -r ${BINS} ${OUT_DIR}
    cmd = [
        "rm -rf {}".format(os.path.join(output_directory, "*")),
            "&&",
    ]

    cmd += [
    "DST={}; (for SRC in {}; do SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST; done)".format(
        output_directory,
        " ".join(input_filepaths), 
        )
    ]


    cmd += [ 
            "&&",
            
        # Statistics
        os.environ["seqkit"],
        "stats",
        "-a",
        "-b",
        "-T",
        "-j {}".format(opts.n_jobs),
        os.path.join(output_directory, "genomes", "*.fa"),
                "|",

        """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: x[:-3]); df.to_csv(sys.stdout, sep="\t")'"""
        ">",
        os.path.join(output_directory,"genome_statistics.tsv"),

            "&&", 

        # CDS
        os.environ["seqkit"],
        "stats",
        "-a",
        "-b",
        "-T",
        "-j {}".format(opts.n_jobs),
        os.path.join(output_directory, "genomes", "*.ffn"),

        "|",

        """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: x[:-4]); df.to_csv(sys.stdout, sep="\t")'"""
        ">",
        os.path.join(output_directory,"gene_statistics.cds.tsv"),

        "&&",
        "rm -rf {}".format(os.path.join(directories["tmp"], "*")),

    ]
    
    if opts.algorithm == "genomad":
        cmd += [ 
            "&&",

        "cp -rf",
        os.path.join(directories[("intermediate",  "1__genomad")], "plasmid_*"),
        output_directory,
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

    if opts.algorithm == "genomad":
        # ==========
        # geNomad
        # ==========
        step = 1

        program = "genomad"

        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


        # Info
        description = "Viral and plasmid identification with geNomad"
        
        # i/o
        input_filepaths = [opts.fasta]

        output_filenames = [
            "virus_summary.tsv",
            "virus_scaffolds.list",
            "viral_taxonomy.tsv",
        ]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_genomad_cmd(**params)
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
    if opts.algorithm == "virfinder":
        # ==========
        # VirFinder
        # ==========
        step = 1

        program = "virfinder"

        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


        # Info
        description = "Viral identification with VirFinder and taxonomy classification with geNomad"
        
        # i/o
        input_filepaths = [opts.fasta]



        output_filenames = [
            "virfinder_output.tsv",
            "virus_scaffolds.list",
            "viral_taxonomy.tsv",
        ]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_virfinder_cmd(**params)
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


    # ==========
    # CheckV
    # ==========
    step = 2

    program = "checkv"

    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Viral verification with CheckV"
    # i/o
    input_filepaths = [
        opts.fasta,
        os.path.join(directories[("intermediate",  "1__{}".format(opts.algorithm))], "virus_scaffolds.list"),
        os.path.join(directories[("intermediate",  "1__{}".format(opts.algorithm))], "viral_taxonomy.tsv"),
    ]
    if opts.algorithm == "genomad":
        input_filepaths += [ 
        os.path.join(directories[("intermediate",  "1__genomad")], "virus_summary.tsv"),
        ]

    output_filenames = [
        # "quality_summary.tsv",
        "filtered/binned.list",
        "filtered/checkv_results.filtered.tsv",
        # "filtered/unbinned.list",
        # "filtered/unbinned.fasta",
        # "filtered/viral_summary.filtered.tsv",
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_checkv_cmd(**params)
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

    # ==========
    # Prodigal
    # ==========
    step = 3

    program = "prodigal-gv"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Viral gene calls via Prodigal-GV"
    # i/o
    input_filepaths = [
        os.path.join(directories[("intermediate",  "2__checkv")],"filtered", "genomes"),
        os.path.join(directories[("intermediate",  "2__checkv")],"filtered", "scaffolds_to_bins.tsv"),

    ]

    output_filenames = [
        os.path.join(output_directory, "*.gff"),
        os.path.join(output_directory, "*.faa"),
        os.path.join(output_directory, "*.ffn"),
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_prodigal_cmd(**params)
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
    # featureCounts
    # ==========
    if opts.bam is not None:
        step += 1

        # Info
        program = "featurecounts"
        program_label = "{}__{}".format(step, program)
        description = "Counting reads"

        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

        # i/o
        input_filepaths = [ 
            opts.fasta,
            os.path.join(directories[("intermediate",  "2__checkv")], "filtered","genomes","*.gff"),
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
                    acceptable_returncodes={0},                    
                    log_prefix=program_label,
        )
    
    # =============
    # Output
    # =============
    step += 1

    program = "output"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories["output"]

    # Info

    description = "Merging results for output"

    # i/o
    input_filenames = [ 
        "scaffolds_to_bins.tsv", 
        "binned.list",
        "unbinned.fasta",
        "unbinned.list",
        "genomes",
        "bins.list",
        "checkv_results.filtered.tsv",
        # "quality_summary.filtered.tsv",
        # "viral_summary.filtered.tsv",

    ]
    input_filepaths = list(map(lambda fn:os.path.join(directories[("intermediate",  "2__checkv")],"filtered", fn), input_filenames))

    if opts.bam is not None:
        input_filepaths += [ 
            os.path.join(directories[("intermediate",  "4__featurecounts")], "featurecounts.orfs.tsv.gz")
        ]

        # "-g {}".format(input_filepaths[1]),
        # "-d {}".format(input_filepaths[2]),
        # "-a {}".format(input_filepaths[3]),

    output_filenames =  [
        "scaffolds_to_bins.tsv", 
        "binned.list",
        # "unbinned.fasta",
        # "unbinned.list",
        "genomes",
        "bins.list",
        "genome_statistics.tsv",
        "checkv_results.filtered.tsv",

    ]
    if opts.bam is not None:
        output_filenames += [ 
            "featurecounts.orfs.tsv.gz",
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
            validate_inputs=False,
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
        "append_geneid_to_prodigal_gff.py",
        "filter_checkv_results.py",
        "virfinder_wrapper.r",
        "genomad_taxonomy_wrapper.py",
        "compile_gff.py",
        # "subset_table.py",
        # "concatenate_dataframes.py",
    }

    required_executables={
                "prodigal-gv",
                "genomad",
                "checkv",
                "seqkit",
                "featureCounts",
     } | accessory_scripts

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
        executables[name] = "'{}'".format(os.path.join(opts.script_directory, "scripts", name)) # Can handle spaces in path
    

    print(format_header( "Adding executables to path from the following source: {}".format(opts.path_config), "-"), file=sys.stdout)
    for name, executable in executables.items():
        if name in required_executables:
            print(name, executable, sep = " --> ", file=sys.stdout)
            os.environ[name] = executable.strip()
    print("", file=sys.stdout)


# Configure parameters
def configure_parameters(opts, directories):
    assert_acceptable_arguments(opts.algorithm, {"genomad", "virfinder"})
    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -f <scaffolds.fasta> -l <contig_identifiers> -n <name> -o <output_directory> [Requires at least 20GB]".format(__program__)

    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-f","--fasta", type=str, required=True, help = "path/to/scaffolds.fasta")
    parser_io.add_argument("-n","--name", type=str, help="Name of sample", required=True)
    parser_io.add_argument("-o","--project_directory", type=str, default="veba_output/binning/viral", help = "path/to/project_directory [Default: veba_output/binning/viral]")
    parser_io.add_argument("-b","--bam", type=str, nargs="+", required=False, help = "path/to/mapped.sorted.bam files separated by spaces.")

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

    # Binning
    parser_binning = parser.add_argument_group('Binning arguments')
    parser_binning.add_argument("-a", "--algorithm", type=str, default="genomad", help="Binning algorithm to use: {genomad, virfinder}  [Default: genomad]")
    parser_binning.add_argument("-m", "--minimum_contig_length", type=int, default=1500, help="Minimum contig length.  [Default: 1500] ")
    parser_binning.add_argument("--include_provirus_detection", action="store_true", help="Include provirus viral detection")

    # Gene models
    parser_genemodels = parser.add_argument_group('Gene model arguments')
    parser_genemodels.add_argument("--prodigal_genetic_code", type=str, default=11, help="Prodigal-GV -g translation table (https://github.com/apcamargo/prodigal-gv) [Default: 11]")

    # geNomad
    parser_genomad = parser.add_argument_group('geNomad arguments\nUsing --relaxed mode by default.  Adjust settings according to the following table: https://portal.nersc.gov/genomad/post_classification_filtering.html#default-parameters-and-presets')
    parser_genomad.add_argument("--genomad_qvalue", type=float, default=1.0, help = "Maximum accepted false discovery rate. [Default: 1.0; 0.0 < x ≤ 1.0]")
    parser_genomad.add_argument("--sensitivity", type=float, default=4.0, help = "MMseqs2 marker search sensitivity. Higher values will annotate more proteins, but the search will be slower and consume more memory. [Default: 4.0; x ≥ 0.0]")
    parser_genomad.add_argument("--splits", type=int, default=0, help = "Split the data for the MMseqs2 search. Higher values will reduce memory usage, but will make the search slower. If the MMseqs2 search is failing, try to increase the number of splits. Also used for VirFinder. [Default: 0; x ≥ 0]")
    parser_genomad.add_argument("--composition", type=str, default="auto", help = "Method for estimating sample composition. (auto|metagenome|virome) [Default: auto]")
    parser_genomad.add_argument("--minimum_score", type=float, default=0.0, help = "Minimum score to flag a sequence as virus or plasmid. By default, the sequence is classified as virus/plasmid if its virus/plasmid score is higher than its chromosome score, regardless of the value. [Default: 0; 0.0 ≤ x ≤ 1.0]")
    parser_genomad.add_argument("--minimum_plasmid_marker_enrichment", type=float, default=-100, help = "Minimum allowed value for the plasmid marker enrichment score, which represents the total enrichment of plasmid markers in the sequence. Sequences with multiple plasmid markers will have higher values than the ones that encode few or no markers.[Default: -100]")
    parser_genomad.add_argument("--minimum_virus_marker_enrichment", type=float, default=-100, help = "Minimum allowed value for the virus marker enrichment score, which represents the total enrichment of plasmid markers in the sequence. Sequences with multiple plasmid markers will have higher values than the ones that encode few or no markers. [Default: -100]")
    parser_genomad.add_argument("--minimum_plasmid_hallmarks", type=int, default=0, help = "Minimum number of plasmid hallmarks in the identified plasmids.  [Default: 0; x ≥ 0]")
    parser_genomad.add_argument("--minimum_virus_hallmarks", type=int, default=0, help = "Minimum number of virus hallmarks in the identified viruses.  [Default: 0; x ≥ 0]")
    parser_genomad.add_argument("--maximum_universal_single_copy_genes", type=int, default=100, help = "Maximum allowed number of universal single copy genes (USCGs) in a virus or a plasmid. Sequences with more than this number of USCGs will not be classified as viruses or plasmids, regardless of their score.  [Default: 100]")
    parser_genomad.add_argument("--genomad_options", type=str, default="", help="geNomad | More options (e.g. --arg 1 ) [Default: '']")

    # VirFinder
    parser_virfinder = parser.add_argument_group('VirFinder arguments')
    parser_virfinder.add_argument("--virfinder_pvalue", type=float, default=0.05, help="VirFinder statistical test threshold [Default: 0.05]")
    parser_virfinder.add_argument("--mmseqs2_evalue", type=float, default=1e-3, help = "Maximum accepted E-value in the MMseqs2 search. Used by genomad annotate when VirFinder is used as binning algorithm [Default: 1e-3]")
    parser_virfinder.add_argument("--use_qvalue", action="store_true", help="Use qvalue (FDR) instead of pvalue")
    parser_virfinder.add_argument("--use_minimal_database_for_taxonomy", action="store_true", help="Use a smaller marker database to annotate proteins. This will make execution faster but sensitivity will be reduced.")
    parser_virfinder.add_argument("--virfinder_options", type=str, default="", help="VirFinder | More options (e.g. --arg 1 ) [Default: '']")

    # CheckV
    parser_checkv = parser.add_argument_group('CheckV arguments')
    parser_checkv.add_argument("--checkv_options", type=str, default="", help="CheckV | More options (e.g. --arg 1 ) [Default: '']")
    parser_checkv.add_argument("--multiplier_viral_to_host_genes", type=int, default=5, help = "Minimum number of viral genes [Default: 5]")
    parser_checkv.add_argument("--checkv_completeness", type=float, default=50.0, help = "Minimum completeness [Default: 50.0]")
    parser_checkv.add_argument("--checkv_quality", type=str, default="High-quality,Medium-quality,Complete", help = "Comma-separated string of acceptable arguments between {High-quality,Medium-quality,Complete} [Default: High-quality,Medium-quality,Complete]")
    parser_checkv.add_argument("--miuvig_quality", type=str, default="High-quality,Medium-quality,Complete", help = "Comma-separated string of acceptable arguments between {High-quality,Medium-quality,Complete} [Default: High-quality,Medium-quality,Complete]")

    # featureCounts
    parser_featurecounts = parser.add_argument_group('featureCounts arguments')
    parser_featurecounts.add_argument("--long_reads", action="store_true", help="featureCounts | Use this if long reads are being used")
    parser_featurecounts.add_argument("--featurecounts_options", type=str, default="", help="featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/")

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
    opts.checkv_database = os.path.join(opts.veba_database, "Classify", "CheckV")

    # Directories
    directories = dict()
    directories["project"] = create_directory(opts.project_directory)
    directories["sample"] = create_directory(os.path.join(directories["project"], opts.name))
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
    print("GenoPype version:", genopype_version, file=sys.stdout) #sys.path[2]
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
