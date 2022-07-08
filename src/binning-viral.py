#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
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


# # VirSorter2
# def get_virsorter2_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
#     # cat [] | seqkit seq -M [] | seqkit grep --pattern-file .list > unbinned.fasta
#     cmd = [
#         "(",
#         "cat",
#         opts.fasta,
#         "|",
#         os.environ["seqkit"],
#         "seq",
#         "-m {}".format(opts.minimum_contig_length),
#         "|",
#         os.environ["seqkit"],
#         "grep",
#         "--pattern-file {}".format(input_filepaths[0]),
#         ">",
#         os.path.join(output_directory, "unbinned.fasta"),
#         ")",
#         "&&"
#         "(",
#         os.environ["virsorter"],
#         "run",
#         "--tmpdir {}".format(directories["tmp"]),
#         "--rm-tmpdir",
#         # "-d {}".format(opts.virsorter2_database),
#         {True:"", False:"--provirus-off"}[bool(opts.include_provirus)],
#         "--use-conda-off",
#         "-j {}".format(opts.n_jobs),
#         "-w {}".format(output_directory),
#         "-i {}".format(os.path.join(output_directory, "unbinned.fasta")),
#         "--min-length {}".format(opts.minimum_contig_length),
#         ")",
#         "&&",
#         "rm {}".format(os.path.join(output_directory, "unbinned.fasta")),
#     ]
#     return cmd

# def get_viralverify_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
#     cmd = [
#         "(",
#         "cat",
#         opts.fasta,
#         "|",
#         os.environ["seqkit"],
#         "seq",
#         "-m {}".format(opts.minimum_contig_length),
#         "|",
#         os.environ["seqkit"],
#         "grep",
#         "--pattern-file {}".format(input_filepaths[0]),
#         ">",
#         os.path.join(directories["tmp"], "unbinned.fasta"),
#         ")",
#         "&&",
#         "rm -rf {}".format(output_directory),
#         "&&",
#         "(",
#         os.environ["viral_verify"],
#         "-i {}".format(os.path.join(directories["tmp"], "unbinned.fasta")),
#         "-o {}".format(output_directory),
#         "-H {}".format(opts.pfam_database),
#         "-t {}".format(opts.n_jobs),
#         opts.viralverify_options, 
#         ")",
#         "&&",
#         "rm {}".format(os.path.join(directories["tmp"], "unbinned.fasta")),
#     ]
#     return cmd

def get_preprocess_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # checkv end_to_end ${FASTA} ${OUT_DIR} -t ${N_JOBS} --restart
    if bool(opts.contig_identifiers):
        cmd = [
            "(",
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
            "|",
            os.environ["seqkit"],
            "grep",
            "--pattern-file {}".format(opts.contig_identifiers),
            ">",
            output_filepaths[0],
            ")",
        ]
    else:
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
            output_filepaths[0],
        ]

    return cmd

def get_virfinder_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [
        "(",
        os.environ["VirFinder_wrapper.R"],
        "-f {}".format(input_filepaths[0]),
        "-o {}".format(output_filepaths[0]),
        opts.virfinder_options,
        ")",
        "&&",
        "(",
        "cat",
        output_filepaths[0],
        "|",
        "python",
        "-c",
        r""" "import sys, pandas as pd; print(*pd.read_csv(sys.stdin, sep='\t', index_col=0).query('pvalue < {}').index, sep='\n')" """.format(opts.virfinder_pvalue),
        ">",
        output_filepaths[1],
        ")",
    ]
    if bool(opts.contig_identifiers) & bool(opts.remove_temporary_fasta):
        cmd += [ 
            "&&",
            "rm -rf {}".format(os.path.join(directories["preprocessing"], "input.fasta")),
        ]
    return cmd

def get_checkv_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
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
        "--pattern-file {}".format(input_filepaths[1]),
        ">",
        os.path.join(directories["tmp"], "virfinder_contigs_for_checkv.fasta"),
        ")",
        "&&",
        "rm -rf {}".format(output_directory),
        "&&",
        "(",
        os.environ["checkv"],
        "end_to_end",
        os.path.join(directories["tmp"], "virfinder_contigs_for_checkv.fasta"),
        output_directory,
        "-t {}".format(opts.n_jobs),
        "-d {}".format(opts.checkv_database),
        "--restart",
        opts.checkv_options, 
        ")",
        "&&",
        "(",
        os.environ["filter_checkv_results.py"],
        "-i {}".format(os.path.join(output_directory, "quality_summary.tsv")),
        "-c {}".format(os.path.join(output_directory, "completeness.tsv")),
        "-f {}".format(os.path.join(directories["tmp"], "virfinder_contigs_for_checkv.fasta")),
        "-o {}".format(os.path.join(output_directory, "filtered")),
        "-m {}".format(opts.minimum_contig_length),
        "--unbinned",
        "-p {}__VIRFINDER__Virus.".format(opts.name),
        "--multiplier_viral_to_host_genes {}".format(opts.multiplier_viral_to_host_genes),
        "--completeness {}".format(opts.checkv_completeness),
        "--checkv_quality {}".format(opts.checkv_quality),
        "--miuvig_quality {}".format(opts.miuvig_quality),
    ] 

    if opts.include_provirus:
        cmd += ["--include_provirus"]

    cmd += [")"]

    
    if opts.remove_temporary_fasta:
        cmd += [ 
        "&&",
        "rm -rf {}".format(os.path.join(directories["tmp"], "virfinder_contigs_for_checkv.fasta")),
        ]

    return cmd

# Prodigal
def get_prodigal_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [
        "(",
        "cat",
        os.path.join(input_filepaths[0], "*.fa"),
        "|",
        os.environ["prodigal"],
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
        ")",
        "&&",
        "(",
        os.environ["partition_gene_models.py"],
        "-i {}".format(input_filepaths[1]),
        "-g {}".format(os.path.join(output_directory, "gene_models.gff")),
        "-d {}".format(os.path.join(output_directory, "gene_models.ffn")),
        "-a {}".format(os.path.join(output_directory, "gene_models.faa")),
        "-o {}".format(os.path.join(directories[("intermediate",  "2__checkv")],"filtered", "genomes")),
        ")",
        "&&",
        "(",
        "rm -f {}".format(os.path.join(output_directory, "gene_models.gff")),
        "&&",
        "rm -f {}".format(os.path.join(output_directory, "gene_models.ffn")),
        "&&",
        "rm -f {}".format(os.path.join(output_directory, "gene_models.faa")),
        ")",

    ]
    return cmd

def get_featurecounts_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # ORF-Level Counts
    cmd = [
    "mkdir -p {}".format(os.path.join(directories["tmp"], "featurecounts")),
    "&&",
        os.environ["featureCounts"],
        "-G {}".format(input_filepaths[0]),
        "-a {}".format(input_filepaths[1]),
        "-o {}".format(os.path.join(output_directory, "featurecounts.orfs.tsv")),
        "-F GTF",
        "--tmpDir {}".format(os.path.join(directories["tmp"], "featurecounts")),
        "-T {}".format(opts.n_jobs),
        "-g gene_id",
        "-t CDS",
        opts.featurecounts_options,
        " ".join(opts.bam),
    "&&",
    "gzip -f {}".format(os.path.join(output_directory, "featurecounts.orfs.tsv")),
    ]
    return cmd

def get_output_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # checkm lineage_wf --tab_table -f checkm_output/${ID}/output.tab --pplacer_threads ${N_JOBS} -t ${N_JOBS} -x fa -r ${BINS} ${OUT_DIR}
    cmd = [
        "rm -rf {}".format(os.path.join(output_directory, "*")),

        # # For pipeline purposes, most likely this file will be overwritten later in this step
        # "cat",
        # opts.fasta,
        # "|",
        # os.environ["seqkit"],
        # "seq",
        # "-m {}".format(opts.minimum_contig_length),
        # ">" ,
        # output_filepaths[2],
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
        "&&",
        "rm -rf {}".format(os.path.join(directories["tmp"], "*")),
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

    # ==========
    # Preprocess
    # ==========
    step = 0

    program = "preprocessing"

    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories["tmp"] 


    # Info
    description = "Preprocessing fasta"
    
    # i/o
    if bool(opts.contig_identifiers):

        input_filepaths = [
            opts.fasta,
            opts.contig_identifiers,
        ]
    else:
        input_filepaths = [
            opts.fasta,
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
                log_prefix=program_label,

    )
        


    # ==========
    # VirFinder
    # ==========
    step = 1

    program = "virfinder"

    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Viral identification with VirFinder"
    
    # i/o
    input_filepaths = [os.path.join(directories["tmp"], "input.fasta")]



    output_filenames = [
        "virfinder_output.tsv",
        "binned.list",
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
        os.path.join(directories[("intermediate",  "1__virfinder")], "binned.list")
    ]

    output_filenames = [
        "quality_summary.tsv",
        "filtered/binned.list",
        "filtered/unbinned.list",
        "filtered/quality_summary.filtered.tsv"
        # "filtered/unbinned.fasta",
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

    program = "prodigal"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Viral gene calls via Prodigal"
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
            os.path.join(directories[("intermediate",  "3__prodigal")], "gene_models.gff"),
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
        "quality_summary.filtered.tsv",

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
        "unbinned.fasta",
        "unbinned.list",
        "genomes",
        "bins.list",
        "genome_statistics.tsv",
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
        # "scaffolds_to_bins.py",
        # "check_scaffolds_to_bins.py",
        "partition_gene_models.py",
        "append_geneid_to_prodigal_gff.py",
        "filter_checkv_results.py",
        "VirFinder_wrapper.R",
    }

    required_executables={
                "prodigal",
                "checkv",
                "seqkit",
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

    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -f <scaffolds.fasta> -l <contig_identifiers> -n <name> -o <output_directory>".format(__program__)

    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-f","--fasta", type=str, required=True, help = "path/to/scaffolds.fasta")
    parser_io.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    parser_io.add_argument("-l","--contig_identifiers", type=str,  help = "path/to/contigs.list [Optional]")
    parser_io.add_argument("-o","--project_directory", type=str, default="veba_output/binning/viral", help = "path/to/project_directory [Default: veba_output/binning/viral]")
    parser_io.add_argument("-b","--bam", type=str, nargs="+", required=False, help = "path/to/mapped.sorted.bam files separated by spaces.")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')

    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
    parser_utility.add_argument("--remove_temporary_fasta", action="store_true", help="If contig identifiers were provided and a fasta is generated, remove this file")

    # parser_utility.add_argument("-c", "--CONDA_PREFIX", type=str, default=None, help = "Set a conda environment")
    # Databases
    parser_databases = parser.add_argument_group('Database arguments')
    parser_databases.add_argument("--veba_database", type=str, default=None, help=f"VEBA database location.  [Default: $VEBA_DATABASE environment variable]")

    # Binning
    parser_binning = parser.add_argument_group('Binning arguments')
    parser_binning.add_argument("-m", "--minimum_contig_length", type=int, default=1500, help="Minimum contig length.  [Default: 1500] ")

    parser_genemodels = parser.add_argument_group('Gene model arguments')
    parser_genemodels.add_argument("--prodigal_genetic_code", type=str, default=11, help="Prodigal -g translation table [Default: 11]")

    

  
    parser_virus = parser.add_argument_group('Virus arguments')
    parser_virus.add_argument("--include_provirus", action="store_true", help="Include provirus viral detection")
    # parser_virus.add_argument("--virsorter2_database", type=str, default="/usr/local/scratch/CORE/jespinoz/db/virsorter/v3/", help="VirSorter2 | More options (e.g. --arg 1 ) [Default: '/usr/local/scratch/CORE/jespinoz/db/virsorter/v3/']")
    # parser_virus.add_argument("--virsorter2_groups", type=str, default="dsDNAphage,NCLDV,ssDNA,lavidaviridae", help="VirSorter2 | More options (e.g. --arg 1 ) [Default: 'dsDNAphage,NCLDV,ssDNA,lavidaviridae']")
    # parser_virus.add_argument("--virsorter2_options", type=str, default="", help="VirSorter2 | More options (e.g. --arg 1 ) [Default: '']")
    # parser_virus.add_argument("--pfam_database", type=str, default="/usr/local/scratch/CORE/jespinoz/db/pfam/v33.1/Pfam-A.hmm", help="PFAM | More options (e.g. --arg 1 ) [Default: '/usr/local/scratch/CORE/jespinoz/db/pfam/v33.1/Pfam-A.hmm']")
    # parser_virus.add_argument("--viralverify_options", type=str, default="", help="ViralVerify | More options (e.g. --arg 1 ) [Default: '']")
    parser_virus.add_argument("--virfinder_pvalue", type=float, default=0.05, help="VirFinder p-value threshold [Default: 0.05]")
    parser_virus.add_argument("--virfinder_options", type=str, default="", help="VirFinder | More options (e.g. --arg 1 ) [Default: '']")
    # parser_virus.add_argument("--checkv_database", type=str, default="/usr/local/scratch/CORE/jespinoz/db/checkv/checkv-db-v1.0", help="CheckV | More options (e.g. --arg 1 ) [Default: '/usr/local/scratch/CORE/jespinoz/db/checkv/checkv-db-v1.0']")
    parser_virus.add_argument("--checkv_options", type=str, default="", help="CheckV | More options (e.g. --arg 1 ) [Default: '']")
    parser_virus.add_argument("--multiplier_viral_to_host_genes", type=int, default=5, help = "Minimum number of viral genes [Default: 5]")
    parser_virus.add_argument("--checkv_completeness", type=float, default=50.0, help = "Minimum completeness [Default: 50.0]")
    parser_virus.add_argument("--checkv_quality", type=str, default="High-quality,Medium-quality,Complete", help = "Comma-separated string of acceptable arguments between {High-quality,Medium-quality,Complete} [Default: High-quality,Medium-quality,Complete]")
    parser_virus.add_argument("--miuvig_quality", type=str, default="High-quality,Medium-quality,Complete", help = "Comma-separated string of acceptable arguments between {High-quality,Medium-quality,Complete} [Default: High-quality,Medium-quality,Complete]")

    # Options
    opts = parser.parse_args()

    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    # Database
    if opts.veba_database is None:
        assert "VEBA_DATABASE" in os.environ, "Please set the following environment variable 'export VEBA_DATABASE=/path/to/veba_database' or provide path to --veba_database"
    else:
        opts.veba_database = os.environ["VEBA_DATABASE"]
    opts.checkv_database = os.path.join(opts.veba_database, "Classify", "CheckV")

    # Directories
    directories = dict()
    directories["project"] = create_directory(opts.project_directory)
    directories["sample"] = create_directory(os.path.join(directories["project"], opts.name))
    # directories["preprocessing"] = create_directory(os.path.join(directories["sample"], "preprocessing"))
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
