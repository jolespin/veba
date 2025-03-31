#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob, random, warnings
from collections import OrderedDict, defaultdict

import pandas as pd

# Soothsayer Ecosystem
from genopype import *
from genopype import __version__ as genopype_version
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2025.3.31"

# Assembly
def get_coverage_cmd( input_filepaths, output_filepaths, output_directory, directories, opts):

    # Command
    cmd = [


    # For pipeline purposes
    "cat",
    opts.fasta,
    "|",
    os.environ["seqkit"],
    "seq",
    "-m {}".format(opts.minimum_contig_length),
    ">" ,
    os.path.join(directories["output"], "unbinned.fasta"),

    "&&",

    # Coverage for MetaBat2 
    os.environ["coverm"],
    "contig",
    "--threads {}".format(opts.n_jobs),
    "--methods metabat",
    "--bam-files",
    " ".join(opts.bam),
    ">",
    output_filepaths[0],

        "&&",

    # # Coverage for MetaCoaAG
    # "cut -f1,4",
    # output_filepaths[0],
    # "|",
    # "tail -n +2",
    # ">",
    # output_filepaths[1],
    """
    
    python -c "import pandas as pd; df = pd.read_csv('{}', sep='\t', index_col=0); df.loc[:,df.columns.map(lambda x: x.endswith(('.bam','.sam')))].to_csv('{}', sep='\t', header=None)"
    
    """.format( 
            output_filepaths[0],
            output_filepaths[1],
    )
    ]
    
    if "metadecoder" in opts.algorithms:
        cmd += [
            
            
        os.environ["metadecoder"],
        "coverage",
        "-b",
        " ".join(opts.bam),
        "-o {}".format(output_filepaths[2]),
        "--threads {}".format(opts.n_jobs),
        opts.metadecoder_coverage_options,
        ]

    return cmd


# Pyrodigal
def get_pyrodigal_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):

    cmd = [
        "cat",
        input_filepaths[0],
        
        "|",
        
        os.environ["seqkit"],
        "seq",
        "-m {}".format(opts.minimum_contig_length),
        
        "|",
                
        os.environ["pyrodigal"],
        "-p meta",
        "-g {}".format(opts.pyrodigal_genetic_code),
        "-f gff",
        "-d {}".format(os.path.join(output_directory, "gene_models.ffn")),
        "-a {}".format(os.path.join(output_directory, "gene_models.faa")),
        "--min-gene {}".format(opts.pyrodigal_minimum_gene_length),
        "--min-edge-gene {}".format(opts.pyrodigal_minimum_edge_gene_length),
        "--max-overlap {}".format(opts.pyrodigal_maximum_gene_overlap_length),
        "-j {}".format(opts.n_jobs),
 
        "|",
        
        os.environ["append_geneid_to_prodigal_gff.py"],
        "-a gene_id",
        ">",
        os.path.join(output_directory, "gene_models.gff"),

    ]
    return cmd

# Metabat2
def get_metabat2_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, prefix, seed):

    cmd = [
        f"rm -rv {output_directory}",
        "&&",
        f"mkdir -p {output_directory}",
        "&&",
        os.environ["binning_wrapper.py"],
        "-a metabat2",
        "-f {}".format(input_filepaths[0]), # scaffolds.fasta
        "-c {}".format(input_filepaths[1]),
        "-o {}".format(output_directory),
        "-m {}".format(opts.minimum_contig_length), # mininimum contig length must be >= 1500 nts for MetaBat2 which is handled by the wrapper
        "-s {}".format(opts.minimum_genome_length), 
        "--n_jobs {}".format(opts.n_jobs),
        "--random_state {}".format(seed),
        "--bin_prefix {}".format(prefix),
        "--remove_bins",
        "--remove_intermediate_files",
        "--metabat2_options {}".format(opts.metabat2_options) if bool(opts.metabat2_options) else "",
    ]
    return cmd

# SemiBin2
def get_semibin2_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, prefix, seed, biome):
    cmd = [
        f"rm -rv {output_directory}",
        "&&",
        f"mkdir -p {output_directory}",
        "&&",
        os.environ["binning_wrapper.py"],
        "-a semibin2",
        "-f {}".format(input_filepaths[0]), # scaffolds.fasta
        "-c {}".format(input_filepaths[1]),
        "-o {}".format(output_directory),
        "-m {}".format(opts.minimum_contig_length), # mininimum contig length must be >= 1000 nts for SemiBin2 which is handled by the wrapper
        "-s {}".format(opts.minimum_genome_length), 
        "--n_jobs {}".format(opts.n_jobs),
        "--random_state {}".format(seed),
        "--bin_prefix {}".format(prefix),
        "--remove_bins",
        "--remove_intermediate_files",
        "--semibin2_orf_finder {}".format(opts.semibin2_orf_finder),
        # "--proteins {}".format(input_filepaths[2]), # Doesn't actually do anything in SemiBin2: https://github.com/BigDataBiology/SemiBin/issues/185
        "--semibin2_biome {}".format(biome),
        "--semibin2_engine {}".format(opts.semibin2_engine),
        "--semibin2_options {}".format(opts.semibin2_options) if bool(opts.semibin2_options) else "",
    ]
    if opts.long_reads:
        cmd += [
            "--semibin2_sequencing_type long_read",
        ]
    else:
        cmd += [
            "--semibin2_sequencing_type short_read",
        ]
    return cmd

# MetaDecoder
def get_metadecoder_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, prefix, seed):

    cmd = [
        f"rm -rv {output_directory}",
        "&&",
        f"mkdir -p {output_directory}",
        "&&",
        os.environ["binning_wrapper.py"],
        "-a metadecoder",
        "-f {}".format(input_filepaths[0]), # scaffolds.fasta
        "-c {}".format(input_filepaths[1]), # Needs special coverage format
        "-o {}".format(output_directory),
        "-m {}".format(opts.minimum_contig_length), # mininimum contig length must be >= 1000 nts for SemiBin2 which is handled by the wrapper
        "-s {}".format(opts.minimum_genome_length), 
        "--n_jobs {}".format(opts.n_jobs),
        "--random_state {}".format(seed),
        "--bin_prefix {}".format(prefix),
        "--remove_bins",
        "--remove_intermediate_files",
        "--proteins {}".format(input_filepaths[2]),
        "--metadecoder_coverage_options {}".format(opts.metadecoder_coverage_options) if bool(opts.metadecoder_coverage_options) else "",
        "--metadecoder_cluster_options {}".format(opts.metadecoder_cluster_options) if bool(opts.metadecoder_cluster_options) else "",
    ]
    return cmd

# MetaCoAG
def get_metacoag_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, prefix, seed):

    cmd = [
        f"rm -rv {output_directory}",
        "&&",
        f"mkdir -p {output_directory}",
        "&&",
        os.environ["binning_wrapper.py"],
        "-a metacoag",
        "-f {}".format(input_filepaths[0]), # scaffolds.fasta
        "-c {}".format(input_filepaths[1]), # coverage
        "-o {}".format(output_directory),
        "-m {}".format(opts.minimum_contig_length), # mininimum contig length must be >= 1000 nts for SemiBin2 which is handled by the wrapper
        "-s {}".format(opts.minimum_genome_length), 
        "--n_jobs {}".format(opts.n_jobs),
        "--random_state {}".format(seed),
        "--bin_prefix {}".format(prefix),
        "--remove_bins",
        "--remove_intermediate_files",
        "--proteins {}".format(input_filepaths[2]),
        "--metacoag_assembler {}".format(opts.metacoag_assembler),
        "--metacoag_graph {}".format(opts.metacoag_graph),
        "--metacoag_paths {}".format(opts.metacoag_paths),
        "--metacoag_options {}".format(opts.metacoag_options) if bool(opts.metacoag_options) else "",
    ]
    return cmd


# # MaxBin2
# def get_maxbin2_107_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, prefix, seed):

#     cmd = [
#         os.environ["binning_wrapper.py"],
#         "-a maxbin2",
#         "--maxbin2_markerset 107",
#         "-f {}".format(input_filepaths[0]), # scaffolds.fasta
#         "-c {}".format(input_filepaths[1]),
#         "-o {}".format(output_directory),
#         "-m {}".format(opts.minimum_contig_length),
#         "-s {}".format(opts.minimum_genome_length), 
#         "--n_jobs {}".format(opts.n_jobs),
#         "--random_state {}".format(seed),
#         "--bin_prefix {}".format(prefix),
#         "--remove_bins",
#         "--remove_intermediate_files",
#         "--maxbin2_options {}".format(opts.maxbin2_options) if bool(opts.maxbin2_options) else "",

#     ]
#     return cmd


# def get_maxbin2_40_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, prefix, seed):

#     cmd = [
#         os.environ["binning_wrapper.py"],
#         "-a maxbin2",
#         "--maxbin2_markerset 40",
#         "-f {}".format(input_filepaths[0]), # scaffolds.fasta
#         "-c {}".format(input_filepaths[1]),
#         "-o {}".format(output_directory),
#         "-m {}".format(opts.minimum_contig_length),
#         "-s {}".format(opts.minimum_genome_length), 
#         "--n_jobs {}".format(opts.n_jobs),
#         "--random_state {}".format(seed),
#         "--bin_prefix {}".format(prefix),
#         "--remove_bins",
#         "--remove_intermediate_files",
#         "--maxbin2_options {}".format(opts.maxbin2_options) if bool(opts.maxbin2_options) else "",
#     ]
#     return cmd


# # CONCOCT
# def get_concoct_cmd( input_filepaths, output_filepaths, output_directory, directories, opts, prefix, seed):

#     cmd = [
#         os.environ["binning_wrapper.py"],
#         "--concoct_fragment_length {}".format(opts.concoct_fragment_length),
#         "--concoct_overlap_length {}".format(opts.concoct_overlap_length),
#         "-a concoct",
#         "-f {}".format(input_filepaths[0]), # scaffolds.fasta
#         "-b {}".format(" ".join(opts.bam)),
#         "-o {}".format(output_directory),
#         "-m {}".format(opts.minimum_contig_length),
#         "-s {}".format(opts.minimum_genome_length), 
#         "--n_jobs {}".format(opts.n_jobs),
#         "--random_state {}".format(seed),
#         "--bin_prefix {}".format(prefix),
#         "--remove_bins",
#         "--remove_intermediate_files",
#         "--concoct_options {}".format(opts.concoct_options) if bool(opts.concoct_options) else "",

#     ]
#     return cmd


# Binette
def get_binette_cmd(input_filepaths, output_filepaths, output_directory, directories, opts, prefix):
    
    # Get non-empty scaffolds to bins
    cmd = [
        f"rm -rv {output_directory}",
        "&&",
        f"mkdir -p {output_directory}",
    ]
    
    cmd += [
    """
    
S2B_FILES="{}"

NON_EMPTY_S2B_FILES=""
for file in $S2B_FILES; do
    if [ -s "$file" ] && [ -n "$(grep -v '^[[:space:]]*$' "$file")" ]; then
        NON_EMPTY_S2B_FILES="$NON_EMPTY_S2B_FILES $file"
    fi
done
    
    """.format(" ".join(input_filepaths[4:]))
    ]

    cmd += [
        # Run Binnette
        os.environ["binette"],
        "-b",
        "${NON_EMPTY_S2B_FILES}",
        "-c",
        input_filepaths[0],
        "-p",
        input_filepaths[3],
        "-o",
        output_directory,
        "-t",
        opts.n_jobs,
        "-w",
        opts.binnette_contamination_weight,
        "--checkm2_db",
        os.path.join(opts.veba_database, "Classify", "CheckM2", "uniref100.KO.1.dmnd"),
        "--low_mem" if opts.checkm2_low_memory else "",
        opts.binnette_options,
        
            "&&",
            
        # Scaffolds to bins intermediate
        os.environ["scaffolds_to_bins.py"],
        "-i",
        os.path.join(output_directory, "final_bins"),
        ">",
        os.path.join(output_directory, "final_bins", "scaffolds_to_bins.tsv"),
        
            "&&",
            
        # -----------------------dev----------------------------            
        # Tiara
        "mkdir -p {}".format(os.path.join(output_directory, "consensus_domain_classification")),
        
        "&&",

        # Subset binned contigs by length filter
        "cat",
        os.path.join(output_directory, "final_bins", "*.fa"),
        "|",
        os.environ["seqkit"],
        "seq",
        "-m {}".format(opts.tiara_minimum_length),
        
        "|",

        # Pipe sequences into Tiara
        os.environ["tiara"],
        "-o {}".format(os.path.join(output_directory, "consensus_domain_classification", "tiara_output.tsv")),
        "--probabilities",
        "-m {}".format(opts.tiara_minimum_length),
        "-t {}".format(opts.n_jobs),
        opts.tiara_options,

        # Predict domain
        "&&",
        
        os.environ["consensus_domain_classification.py"],
        "-i {}".format(os.path.join(directories["output"], "scaffolds_to_bins.tsv")),
        "-t {}".format(os.path.join(output_directory, "consensus_domain_classification", "tiara_output.tsv")),
        "-o {}".format(os.path.join(output_directory, "consensus_domain_classification")),
        "--logit_transform {}".format(opts.logit_transform),
        # -----------------------dev----------------------------
        
        "&&",
        
        
        # Add filtered directory
        os.environ["filter_binette_results.py"],
        "-i",
        output_directory,
        "-o",
        os.path.join(output_directory, "filtered"),
        "-f",
        input_filepaths[0],
        # "--unbinned",
        "-m",
        opts.minimum_contig_length,
        "--completeness",
        opts.checkm2_completeness,
        "--contamination",
        opts.checkm2_contamination,
        "--bin_prefix",
        prefix,
        
 
            "&&",
            
        # Scaffolds to Bins
        "ls",
        os.path.join(output_directory, "final_bins", "*.fa"),
        "|",
        os.environ["scaffolds_to_bins.py"],
        "-x",
        "fa",
        ">",
        os.path.join(output_directory, "scaffolds_to_bins.tsv"),
        
            "&&",
            
        # Partition gene models
        os.environ["partition_gene_models.py"],
        "-i",
        os.path.join(output_directory, "filtered", "scaffolds_to_bins.tsv"),
        "-g",
        input_filepaths[1], #os.path.join(directories[("intermediate",  "2__pyrodigal")], "gene_models.gff"),
        "-d",
        input_filepaths[2], #os.path.join(directories[("intermediate",  "2__pyrodigal")], "gene_models.ffn"),
        "-a",
        input_filepaths[3], #os.path.join(directories[("intermediate",  "2__pyrodigal")], "gene_models.faa"),
        "-o",
        os.path.join(output_directory, "filtered", "genomes"),
        "--use_mag_as_description",

            "&&",
            
        # Cleanup
        "mv",
        os.path.join(output_directory, "final_bins_quality_reports.tsv"),
        os.path.join(output_directory, "quality_reports.tsv"),
        
        "&&",
        
        "rm",
        "-rfv",
        os.path.join(output_directory, "final_bins"),
        os.path.join(output_directory, "temporary_files", "assembly_proteins.faa.gz"),
        os.path.join(output_directory, "temporary_files", "*.fxi"),
        
        "&&",
        
        "cat",
        input_filepaths[0],
        "|",
        os.environ["seqkit"],
        "seq",
        "-m {}".format(opts.minimum_contig_length),
        "|",
        os.environ["seqkit"],
        "grep",
        "--pattern-file {}".format(os.path.join(output_directory, "filtered", "unbinned.list")),
        ">",
        output_filepaths[-1],

    ]
    return cmd
    
# Tiara
# def get_tiara_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
#     cmd = [
#         "mkdir -p {}".format(os.path.join(output_directory, "consensus_domain_classification")),
        
#         # Get scaffolds to bins
#         "&&",
        
#         "cat",
#         " ".join(input_filepaths[1:]),
#         ">",
#         os.path.join(directories["output"], "scaffolds_to_bins.tsv"),
        
#         # Get binned contigs
        
#         "&&",
        
#         "cut",
#         "-f",
#         "1",
#         os.path.join(directories["output"], "scaffolds_to_bins.tsv"),
#         ">",
#         os.path.join(directories["output"], "binned.list"),
        
#         # Subset binned contigs by length filter
#         "&&",
#         os.environ["seqkit"],
#         "seq",
#         "-m {}".format(opts.tiara_minimum_length),
#         input_filepaths[0],
        
#         "|",
        
#         os.environ["seqkit"],
#         "grep",
#         "-f",
#         os.path.join(directories["output"], "binned.list"),
        
#         "|",

#         # Pipe sequences into Tiara
#         os.environ["tiara"],
#         "-o {}".format(os.path.join(output_directory, "consensus_domain_classification", "tiara_output.tsv")),
#         "--probabilities",
#         "-m {}".format(opts.tiara_minimum_length),
#         "-t {}".format(opts.n_jobs),
#         opts.tiara_options,

#         # Predict domain
#         "&&",
        
#         os.environ["consensus_domain_classification.py"],
#         "-i {}".format(os.path.join(directories["output"], "scaffolds_to_bins.tsv")),
#         "-t {}".format(os.path.join(output_directory, "consensus_domain_classification", "tiara_output.tsv")),
#         "-o {}".format(os.path.join(output_directory, "consensus_domain_classification")),
#         "--logit_transform {}".format(opts.logit_transform),
#     ]
#     return cmd

# def get_consolidate_prokaryotic_and_remove_eukaryotic_genomes_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
#     cmd = [
        
#     ]
#     return cmd

# barrnap
def get_barrnap_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    cmd = [
        # "cat",
        # input_filepaths[0],
        # ">",
        # os.path.join(directories["tmp"], "genomes_to_domain.tsv"),

"""
OUTPUT_DIRECTORY={}
FP={} 
for DOMAIN in $(cut -f2 $FP | sort -u); 
do 
    DOMAIN_ABBREVIATION=$(echo $DOMAIN | python -c 'import sys; print(sys.stdin.read().lower()[:3])')
    
    # Get MAGs for each domain (not all will have passed QC)
    for ID in $(cat $FP | grep $DOMAIN | cut -f1)
    do 
        GENOME_FASTA=$(ls {}) || GENOME_FASTA=""
        if [ -e "$GENOME_FASTA" ]; then
            >$OUTPUT_DIRECTORY/$ID.rRNA
            >$OUTPUT_DIRECTORY/$ID.rRNA.gff
            {} --kingdom $DOMAIN_ABBREVIATION --threads {} --lencutoff {} --reject {} --evalue {} --outseq $OUTPUT_DIRECTORY/$ID.rRNA $GENOME_FASTA  | {} > $OUTPUT_DIRECTORY/$ID.rRNA.gff
            #rm $GENOME_FASTA.fai
        fi
    done
done
""".format(
        output_directory,
        input_filepaths[0],
        os.path.join(os.path.split(input_filepaths[1])[0],"$ID.fa"), # Is this inefficient?
        os.environ["barrnap"],
        opts.n_jobs,
        opts.barrnap_length_cutoff, 
        opts.barrnap_reject, 
        opts.barrnap_evalue,
        os.environ["append_geneid_to_barrnap_gff.py"],
    ),
    ]
    return cmd

# tRNAscan-SE
def get_trnascan_cmd(input_filepaths, output_filepaths, output_directory, directories, opts): # Should use GNU parallel


            # os.environ["cmsearch"],
            # "-o {}".format(os.path.join(output_directory, "{}.out".format(id_model))),
            # "-A {}".format(os.path.join(output_directory, "{}.aln".format(id_model))),
            # "--tblout {}".format(os.path.join(output_directory, "{}.tblout".format(id_model))),
            # "--cpu {}".format(opts.n_jobs),
            # "-g",
            # "--notrunc",
            # "--mid",
            # opts.cmsearch_options,
            # filepath,
            # input_filepaths[0],

    cmd = [
        # "cat",
        # input_filepaths[0],
        # ">",
        # os.path.join(directories["tmp"], "genomes_to_domain.tsv"),


"""
OUTPUT_DIRECTORY={}
FP={} 
for DOMAIN in $(cut -f2 $FP | sort -u); 
do 
    DOMAIN_ABBREVIATION=$(echo $DOMAIN | python -c 'import sys; print(sys.stdin.read().upper()[:1])')
    
    # Get MAGs for each domain (not all will have passed QC)
    for ID in $(cat $FP | grep $DOMAIN | cut -f1)
    do 
        GENOME_FASTA=$(ls {}) || GENOME_FASTA=""
        if [ -e "$GENOME_FASTA" ]; then
            >$OUTPUT_DIRECTORY/$ID.tRNA
            >$OUTPUT_DIRECTORY/$ID.tRNA.gff
            >$OUTPUT_DIRECTORY/$ID.tRNA.struct
            >$OUTPUT_DIRECTORY/$ID.tRNA.txt
            
            TRNA_FASTA=$OUTPUT_DIRECTORY/$ID.tRNA

            if [[ -s "$TRNA_FASTA" ]]; 
                then
                    echo "[Skipping] [tRNAscan-SE] $GENOME_FASTA because tRNA fasta exists and is not empty"
                else
                    echo "[Running] [tRNAscan-SE] $GENOME_FASTA"
                    {} -$DOMAIN_ABBREVIATION --forceow --progress --threads {} --fasta $OUTPUT_DIRECTORY/$ID.tRNA --gff $OUTPUT_DIRECTORY/$ID.tRNA.gff --struct $OUTPUT_DIRECTORY/$ID.tRNA.struct {} $GENOME_FASTA > $OUTPUT_DIRECTORY/$ID.tRNA.txt
            fi    
        fi  
    done
done
""".format(
        output_directory,
        input_filepaths[0],
        os.path.join(os.path.split(input_filepaths[1])[0],"$ID.fa"),
        os.environ["tRNAscan-SE"],
        opts.n_jobs,
        opts.trnascan_options,
    )
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
        "-T {}".format(min(64, opts.n_jobs)), # The maximum number of threads featureCounts can use is 64 so any more will throw this error: "Value for argumant -T is out of range: 1 to 64"
        "-g gene_id",
        "-t CDS",
        "-L" if opts.long_reads else "-p --countReadPairs",
        opts.featurecounts_options,
        " ".join(opts.bam),
    "&&",
    "gzip -f {}".format(os.path.join(output_directory, "featurecounts.orfs.tsv")),
    ]
    return cmd


def get_consolidate_cmd(input_filepaths, output_filepaths, output_directory, directories, opts, step):


    cmd = [
        """
mkdir -p {}
S2B=$(ls {}) || (echo 'No genomes have been detected' && exit 1)

""".format(
    os.path.join(output_directory, "genomes"),
    os.path.join(directories["intermediate"], "*__binette",  "filtered", "scaffolds_to_bins.tsv"),
    ),
    ]

    cmd += [ 

        # scaffolds_to_bins.tsv
        # "cat",
        # os.path.join(directories["intermediate"], "*__binette",  "filtered", "scaffolds_to_bins.tsv"), 
        # ">",
        # os.path.join(output_directory, "scaffolds_to_bins.tsv"),

            # "&&",

        # bins.list
        "cat",
        os.path.join(directories["intermediate"], "*__binette",  "filtered", "bins.list"), 
        ">",
        os.path.join(output_directory, "bins.list"),

            # "&&",

        # # binned.list
        # "cat",
        # os.path.join(directories["intermediate"], "*__binette",  "filtered", "binned.list"), 
        # ">",
        # os.path.join(output_directory, "binned.list"),

            "&&",

        # checkm2_results.filtered.tsv 
        os.environ["concatenate_dataframes.py"],
        "-a 0",
        # "-e",
        os.path.join(directories["intermediate"], "*__binette",  "filtered", "checkm2_results.filtered.tsv"), 
        ">",
        os.path.join(output_directory, "checkm2_results.filtered.tsv"),
    ]

    # Genomes (.fa, .ffn, .faa, .gff)

    cmd += [ 
            "&&",

        "DST={}; for SRC in {}; do SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST; done".format(
        os.path.join(output_directory,"genomes"),
        os.path.join(directories["intermediate"], "*__binette", "filtered", "genomes", "*"),
    ),
    ]

    # GFF
    cmd += [ 
"""

DIR_RRNA={}
DIR_TRNA={}
OUTPUT_DIRECTORY={}
mkdir -p $OUTPUT_DIRECTORY

for GENOME_FASTA in $(ls {});
do
    ID=$(basename $GENOME_FASTA .fa)
    DIR_GENOME=$(dirname $GENOME_FASTA)
    GFF_CDS=$DIR_GENOME/$ID.gff
    GFF_RRNA=$DIR_RRNA/$ID.rRNA.gff
    GFF_TRNA=$DIR_TRNA/$ID.tRNA.gff
    GFF_OUTPUT=$OUTPUT_DIRECTORY/$ID.gff
    >$GFF_OUTPUT.tmp
    {} -f $GENOME_FASTA -o $GFF_OUTPUT.tmp -n $ID -c $GFF_CDS -r $GFF_RRNA -t $GFF_TRNA -d Prokaryotic
    mv $GFF_OUTPUT.tmp $GFF_OUTPUT
done

""".format(
    directories[("intermediate", "{}__barrnap".format(step-3))],
    directories[("intermediate", "{}__trnascan-se".format(step-2))],
    os.path.join(output_directory,"genomes"),
    os.path.join(directories["intermediate"], "*__binette", "filtered", "genomes", "*.fa"),
    os.environ["compile_gff.py"],
    )
    ]

    # rRNA

    cmd += [ 

        "DST={}; for SRC in $(ls {}); do SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST; done".format(
        os.path.join(output_directory,"genomes"),
        os.path.join(directories[("intermediate", "{}__barrnap".format(step-3))], "*.rRNA"),
    ),
    ]

    # tRNA

    cmd += [ 
            "&&",

        "DST={}; for SRC in $(ls {}); do SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST; done".format(
        os.path.join(output_directory,"genomes"),
        os.path.join(directories[("intermediate", "{}__trnascan-se".format(step-2))], "*.tRNA"),
    ),
    ]


    
    # featureCounts
    cmd += [
            "&&",
        "SRC={}; DST={}; SRC=$(realpath --relative-to $DST $SRC); ln -sf $SRC $DST".format(
            os.path.join(directories[("intermediate", "{}__featurecounts".format(step-1))], "featurecounts.orfs.tsv.gz"),
            output_directory,
        ),
        ]
        
    # SeqKit
    cmd += [ 
            "&&",
        # Statistics
        # Assembly
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
            
        # rRNA
        os.environ["seqkit"],
        "stats",
        "-a",
        "-b",
        "-T",
        "-j {}".format(opts.n_jobs),
        os.path.join(output_directory, "genomes", "*.rRNA"),

        "|",

        """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: x[:-5]); df.to_csv(sys.stdout, sep="\t")'"""
        ">",
        os.path.join(output_directory,"gene_statistics.rRNA.tsv"),

                "&&",
            
        # tRNA
        os.environ["seqkit"],
        "stats",
        "-a",
        "-b",
        "-T",
        "-j {}".format(opts.n_jobs),
        os.path.join(output_directory, "genomes", "*.tRNA"),

        "|",

        """python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: x[:-5]); df.to_csv(sys.stdout, sep="\t")'"""
        ">",
        os.path.join(output_directory,"gene_statistics.tRNA.tsv"),

        # Binned/Unbinned
            "&&",

        "cat",
        opts.fasta,
        "|",
        os.environ["seqkit"],
        "grep",
        "-v",
        "-f {}".format(os.path.join(output_directory, "binned.list")),
        "-j {}".format(opts.n_jobs),
        "|",
        os.environ["seqkit"],
        "seq",
        "-j {}".format(opts.n_jobs),
        "-m {}".format(opts.minimum_contig_length),
        ">",
        os.path.join(directories["tmp"], "unbinned.fasta"),

            "&&",

        "mv",
        os.path.join(directories["tmp"], "unbinned.fasta"),
        os.path.join(output_directory,"unbinned.fasta"), 
    
            "&&",

        "rm -rf {}".format(os.path.join(directories["tmp"],"*")),
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
    # Coverage
    # ==========
    step = 1

    program = "coverage"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))

    # Info
    description = "Calculating coverage for assembly"

    input_filepaths = [
            *opts.bam,
        ]

    output_filenames = ["coverage_metabat.tsv", "coverage_noheader.tsv"]
    if "metadecoder" in opts.algorithms:
        output_filenames.append("coverage_metadecoder.tsv")

    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_coverage_cmd(**params)

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
    # Pyrodigal
    # ==========
    step = 2

    program = "pyrodigal"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Gene calls via Pyrodigal"
    # i/o
    input_filepaths = [
        opts.fasta,
    ]

    output_filenames = [
        "gene_models.gff",
        "gene_models.faa",
        "gene_models.ffn",
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_pyrodigal_cmd(**params)
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

    # Iterative mode
    input_fasta = opts.fasta
    for iteration in range(1, opts.n_iter+1):
        steps = dict()

        seed = opts.random_state
        if opts.random_state == -1:
            seed = random.choice(range(1,10000))
        if opts.random_state == 0:
            seed = iteration

        for algorithm in opts.algorithms:
            # ==========
            # MetaBat2
            # ==========
            if algorithm == "metabat2":

                step  += 1

                program = "binning_metabat2"
                program_label = "{}__{}".format(step, program)
                
                # Add to directories
                output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


                # Info
                description = "Binning via MetaBat2 [Iteration={}]".format(iteration)

                # i/o
                input_filepaths = [
                    input_fasta,
                    os.path.join(directories[("intermediate",  "1__coverage")], "coverage_metabat.tsv"),
                ]

                output_filenames = [
                    # "bin*.fa", 
                    "scaffolds_to_bins.tsv",
                    ]
                output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

                params = {
                    "input_filepaths":input_filepaths,
                    "output_filepaths":output_filepaths,
                    "output_directory":output_directory,
                    "opts":opts,
                    "directories":directories,
                    "prefix":"{}__METABAT2__{}.{}__".format(opts.name, "P", iteration),
                    "seed":seed,
                }



                cmd = get_metabat2_cmd(**params)
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
                            acceptable_returncodes={0,1,2},
                            log_prefix=program_label,

                )

                steps[program] = step
                
            # ==========
            # SemiBin2
            # ==========
            if isinstance(algorithm, tuple):
                algorithm, biome = algorithm

                step  += 1

                program = f"binning_semibin2-{biome}"
                program_label = "{}__{}".format(step, program)
                
                # Add to directories
                output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


                # Info
                description = "Binning via SemiBin2 ({}) [Iteration={}]".format(biome, iteration)

                # i/o
                input_filepaths = [
                    input_fasta,
                    os.path.join(directories[("intermediate",  "1__coverage")], "coverage_metabat.tsv"),
                    os.path.join(directories[("intermediate",  "2__pyrodigal")], "gene_models.faa"),
                ]

                output_filenames = [
                    # "bin*.fa", 
                    "scaffolds_to_bins.tsv",
                    ]
                output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

                params = {
                    "input_filepaths":input_filepaths,
                    "output_filepaths":output_filepaths,
                    "output_directory":output_directory,
                    "opts":opts,
                    "directories":directories,
                    "prefix":"{}__SEMIBIN2-{}__{}.{}__".format(opts.name, biome.upper(),"P", iteration),
                    "seed":seed,
                    "biome":biome,
                }



                cmd = get_semibin2_cmd(**params)
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
                            acceptable_returncodes={0,1,2},
                            log_prefix=program_label,

                )

                steps[program] = step
                
            # ==========
            # MetaDecoder
            # ==========
            if algorithm == "metadecoder":

                step  += 1

                program = "binning_metadecoder"
                program_label = "{}__{}".format(step, program)
                
                # Add to directories
                output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


                # Info
                description = "Binning via MetaDecoder [Iteration={}]".format(iteration)

                # i/o
                input_filepaths = [
                    input_fasta,
                    os.path.join(directories[("intermediate",  "1__coverage")], "coverage_metadecoder.tsv"),
                    os.path.join(directories[("intermediate",  "2__pyrodigal")], "gene_models.faa"),
                ]

                output_filenames = [
                    # "bin*.fa", 
                    "scaffolds_to_bins.tsv",
                    ]
                output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

                params = {
                    "input_filepaths":input_filepaths,
                    "output_filepaths":output_filepaths,
                    "output_directory":output_directory,
                    "opts":opts,
                    "directories":directories,
                    "prefix":"{}__METADECODER__{}.{}__".format(opts.name, "P", iteration),
                    "seed":seed,
                }

                cmd = get_metadecoder_cmd(**params)
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
                            acceptable_returncodes={0,1,2},
                            log_prefix=program_label,

                )

                steps[program] = step
                
            # ==========
            # MetaCoAG
            # ==========
            if algorithm == "metacoag":

                step  += 1

                program = "binning_metacoag"
                program_label = "{}__{}".format(step, program)
                
                # Add to directories
                output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


                # Info
                description = "Binning via MetaCoAG [Iteration={}]".format(iteration)

                # i/o
                input_filepaths = [
                    input_fasta,
                    os.path.join(directories[("intermediate",  "1__coverage")], "coverage_noheader.tsv"),
                    os.path.join(directories[("intermediate",  "2__pyrodigal")], "gene_models.faa"),
                ]
                if opts.metacoag_graph != "auto":
                    input_filepaths.append(opts.metacoag_graph)
                if opts.metacoag_paths != "auto":
                    input_filepaths.append(opts.metacoag_paths)

                output_filenames = [
                    # "bin*.fa", 
                    "scaffolds_to_bins.tsv",
                    ]
                output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

                params = {
                    "input_filepaths":input_filepaths,
                    "output_filepaths":output_filepaths,
                    "output_directory":output_directory,
                    "opts":opts,
                    "directories":directories,
                    "prefix":"{}__METACOAG__{}.{}__".format(opts.name, "P", iteration),
                    "seed":seed,
                }

                cmd = get_metacoag_cmd(**params)
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
                            acceptable_returncodes={0,1,2},
                            log_prefix=program_label,

                )

                steps[program] = step



        # ==========
        # Binette
        # ==========
        step  += 1

        program = "binette"
        program_label = "{}__{}".format(step, program)
        # Add to directories
        output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


        # Info
        description = "Evaluation via Binette [Iteration={}]".format(iteration)
        # i/o
        input_filepaths = [
            input_fasta,
            os.path.join(directories[("intermediate",  "2__pyrodigal")], "gene_models.gff"),
            os.path.join(directories[("intermediate",  "2__pyrodigal")], "gene_models.ffn"),
            os.path.join(directories[("intermediate",  "2__pyrodigal")], "gene_models.faa"),
        ]
        for algorithm in opts.algorithms:
            program = "binning_{}".format(algorithm)
            s2b = os.path.join(directories[("intermediate","{}__{}".format(steps[program], program))], "scaffolds_to_bins.tsv")
            input_filepaths.append(s2b)
            

        output_filenames = [
            "filtered/checkm2_results.filtered.tsv",
            "filtered/*.list",
        ]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))
        output_filepaths += [ 
            os.path.join(directories["tmp"], "unbinned_{}.fasta".format(iteration)),
        ]


        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
            "prefix":"{}__BINETTE__{}.{}__".format(opts.name, "P", iteration),

        } 

        cmd = get_binette_cmd(**params)
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
                    acceptable_returncodes={0,1,2},                    
                    log_prefix=program_label,
                    # acceptable_returncodes= {0,1},

        )


        steps[program] = step

    
        # Reset input_fasta
        input_fasta = os.path.join(directories["tmp"], "unbinned_{}.fasta".format(iteration))
        
    # # ==========
    # # Tiara
    # # ==========
    # step  += 1

    # program = "tiara"
    # program_label = "{}__{}".format(step, program)
    
    # # Add to directories
    # output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # # Info
    # description = "Consensus domain classification"
    # # i/o
    # input_filepaths = [
    #     opts.fasta,
    #     os.path.join(directories["intermediate"], "*__binette",  "filtered", "scaffolds_to_bins.tsv"),
    # ]


    # output_filepaths = [
    #     os.path.join(directories["output"], "scaffolds_to_bins.tsv"),
    #     os.path.join(directories["output"], "binned.list"),
    #     os.path.join(output_directory, "tiara_results.tsv"),
    #     os.path.join(output_directory,  "consensus_domain_classification", "predictions.tsv"),


    # ]

    # params = {
    #     "input_filepaths":input_filepaths,
    #     "output_filepaths":output_filepaths,
    #     "output_directory":output_directory,
    #     "opts":opts,
    #     "directories":directories,
    # } 

    # cmd = get_tiara_cmd(**params)
    # pipeline.add_step(
    #             id=program_label,
    #             description = description,
    #             step=step,
    #             cmd=cmd,
    #             input_filepaths = input_filepaths,
    #             output_filepaths = output_filepaths,
    #             validate_inputs=False,
    #             validate_outputs=False,
    #             errors_ok=False,
    #             acceptable_returncodes={0},                    
    #             log_prefix=program_label,
    #             # acceptable_returncodes= {0,1},

    # )


    steps[program] = step
    
    # ==========
    # barrnap
    # ==========
    step  += 1

    program = "barrnap"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Detecting rRNA genes"
    # i/o
    input_filepaths = [
        os.path.join(directories["intermediate"], f"{step - 1}__tiara",  "consensus_domain_classification", "predictions.tsv"),
        os.path.join(directories["intermediate"], "*__binette", "filtered", "genomes", "*"),
    ]

    output_filenames = [
        "*.rRNA",
        "*.rRNA.gff",
    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
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
                acceptable_returncodes={0},                    
                log_prefix=program_label,
                # acceptable_returncodes= {0,1},

    )


    steps[program] = step

    # ==========
    # tRNASCAN-se
    # ==========
    step  += 1

    program = "trnascan-se"
    program_label = "{}__{}".format(step, program)
    # Add to directories
    output_directory = directories[("intermediate",  program_label)] = create_directory(os.path.join(directories["intermediate"], program_label))


    # Info
    description = "Detecting tRNA genes"
    # i/o
    input_filepaths = [
        os.path.join(directories["intermediate"], f"{step - 2}__tiara",  "consensus_domain_classification", "predictions.tsv"),
        os.path.join(directories["intermediate"], "*__binette", "filtered", "genomes", "*"),
    ]

    output_filenames = [
        "*.tRNA",
        "*.tRNA.gff",
        "*.tRNA.struct",

    ]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
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
                acceptable_returncodes={0},                    
                log_prefix=program_label,
                # acceptable_returncodes= {0,1},

    )


    steps[program] = step


    # ==========
    # featureCounts
    # ==========
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
        os.path.join(directories[("intermediate",  "2__pyrodigal")], "gene_models.gff"),
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
                acceptable_returncodes={0,1,2},                    
                log_prefix=program_label,
    )
    
    
    # =============
    # Consolidate
    # =============
    step += 1

    program = "consolidate"
    program_label = "{}__{}".format(step, program)

    # Add to directories
    output_directory = directories["output"]

    # Info

    description = "Consolidate output files"

    # i/o
    input_filepaths = [
        # os.path.join(directories["intermediate"], "*__binette",  "filtered", "scaffolds_to_bins.tsv"),
        # os.path.join(directories["intermediate"], "*__binette",  "filtered", "bins.list"),
        # os.path.join(directories["intermediate"], "*__binette",  "filtered", "binned.list"),
        os.path.join(directories["intermediate"], "*__binette",  "filtered", "checkm2_results.filtered.tsv"),
        os.path.join(directories["intermediate"], "*__binette", "filtered", "genomes", "*"),
        os.path.join(directories[("intermediate", "{}__featurecounts".format(step-1))], "featurecounts.orfs.tsv.gz"),

        # Can't assume these are not empty
        # os.path.join(directories[("intermediate", "{}__trnascan-se".format(step-2))], "*.tRNA"), 
        # os.path.join(directories[("intermediate", "{}__barrnap".format(step-3))], "*.rRNA"),
    ]

    output_filenames =  [
        "scaffolds_to_bins.tsv", 
        "bins.list",
        "binned.list",
        "unbinned.fasta",
        "genomes",
        "checkm2_results.filtered.tsv",
        "featurecounts.orfs.tsv.gz",
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
    "step":step,
    }

    cmd = get_consolidate_cmd(**params)
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
                "binning_wrapper.py",
                "scaffolds_to_bins.py",
                # "check_scaffolds_to_bins.py",
                "filter_binette_results.py",
                "partition_gene_models.py",
                "append_geneid_to_prodigal_gff.py",
                "append_geneid_to_barrnap_gff.py",
                "consensus_domain_classification.py",
                "concatenate_dataframes.py",
                "subset_table.py",
                "compile_gff.py",
                }

    required_executables={
                "coverm",
                "binette",
                "pyrodigal",
                "checkm2",
                "tiara",
                "seqkit",
                "featureCounts",
                "barrnap",
                "tRNAscan-SE",
                # "parallel",
 
     } | accessory_scripts

    
    for algorithm in opts.algorithms:
        if isinstance(algorithm, tuple):
            if algorithm[0] != ("semibin2"):
                raise f"If algorithm is tuple, then the first element must be semibin2: {algorithm}"
            required_executables.add("SemiBin2")
        elif algorithm == "metadecoder":
            required_executables.add("metadecoder")
        elif algorithm == "metacoag":
            required_executables.add("metacoag")
  
    if opts.path_config == "CONDA_PREFIX":
        executables = dict()
        for name in sorted(required_executables):
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
    # Set environment variables
    semibin2_biomes = {'ocean', 'wastewater', 'global', 'pig_gut', 'human_oral', 'cat_gut', 'soil', 'chicken_caecum', 'human_gut', 'built_environment', 'dog_gut', 'mouse_gut', 'NONE'}
    opts.algorithms = set(map(str.strip, opts.algorithms.split(",")))
    choices = {"metabat2", "semibin2", "metadecoder", "metacoag"} | set(map(lambda biome: f"semibin2-{x}", semibin2_biomes))
    assert opts.algorithms <= choices, "Unrecognized algorithm(s): {}".format(opts.algorithms - choices)
    algorithms_formatted = list()
    for algorithm in opts.algorithms:
        if algorithm.startswith("semibin2"):
            if "-" in algorithm:
                biome = algorithm.split("-")[1]
            else:
                warnings.warn("No biome specified for semibin2.  Defaulting to global.")
                biome = "global"
            algorithms_formatted.append(("semibin2", biome))
        else:
            algorithms_formatted.append(algorithm)
    opts.algorithms = sorted(algorithms_formatted)
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
    parser_io = parser.add_argument_group('Required I/O arguments')
    parser_io.add_argument("-f","--fasta", type=str, required=True, help = "path/to/scaffolds.fasta")
    parser_io.add_argument("-b","--bam", type=str, nargs="+", required=True, help = "path/to/mapped.sorted.bam files separated by spaces.")
    parser_io.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    parser_io.add_argument("-o","--project_directory", type=str, default="veba_output/binning/prokaryotic/", help = "path/to/project_directory [Default: veba_output/binning/prokaryotic/]")
    parser_io.add_argument("-L", "--long_reads", action="store_true", help="Use this if long reads are being used")


    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')

    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("-I", "--n_iter", type=int, default=3, help = "Number of iterations to run binning. Use -1 for brute force (-1 isn't ready yet) [Default: 3]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Use -1 for completely random. Use 0 for consecutive random states.  Use any other positive integer for the same random state for all iterations [Default: 0]")
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
    parser_utility.add_argument("--tmpdir", type=str, help="Set temporary directory") 

    # Databases
    parser_databases = parser.add_argument_group('Database arguments')
    parser_databases.add_argument("--veba_database", type=str,  help=f"VEBA database location.  [Default: $VEBA_DATABASE environment variable]")

    # Binning
    parser_binning = parser.add_argument_group('Binning arguments')
    parser_binning.add_argument("-a", "--algorithms", type=str, default="metabat2,semibin2-global,metadecoder,metacoag", help='Comma separated list of binning algorithms.  Choose from {"metabat2", "metadecoder", "metacoag", "semibin2-[biome]"} where [biome] is one of [ocean, wastewater, global, pig_gut, human_oral, cat_gut, soil, chicken_caecum, human_gut, built_environment, dog_gut, mouse_gut, NONE] where NONE is usedto implement Semi-Supervised training (takes longer with more compute).  If MEGAHIT assembly was used then metacoag will fail if the --fasta does not match the --metacoag_graph exactly which will be fixed in future versions. Do not use metacoag with MEGAHIT assembly if  you have performed viral binning prior to prokaryotic binning. [Default: metabat2,semibin2-global,metadecoder,metacoag]')
    parser_binning.add_argument("-m", "--minimum_contig_length", type=int, default=1500, help="Minimum contig length.  Anything under 2500 will default to 2500 for MetaBat2 [Default: 1500] ")
    parser_binning.add_argument("-s", "--minimum_genome_length", type=int, default=200000, help="Minimum genome length.  [Default: 200000]")
    parser_binning.add_argument("--retain_intermediate_bins",action="store_true",help='Retain intermediate bins in fasta.')

    # Metabat2
    parser_metabat2 = parser.add_argument_group('Metabat2 arguments')
    parser_metabat2.add_argument("--metabat2_options", type=str, default="", help="MetaBat2 | More options (e.g. --arg 1 ) [Default: ''] | https://bitbucket.org/berkeleylab/metabat/src/master/")

    # SemiBin2
    parser_semibin2 = parser.add_argument_group('SemiBin2 arguments')
    # parser_semibin2.add_argument("--semibin2_biome", type=str, choices={'ocean', 'wastewater', 'global', 'pig_gut', 'human_oral', 'cat_gut', 'soil', 'chicken_caecum', 'human_gut', 'built_environment', 'dog_gut', 'mouse_gut', 'NONE'}, default="global", help="SemiBin2 | Biome/environment for the built-in model.  Use 'NONE' to implement Semi-Supervised training (takes longer with more compute) [Default: global]")
    parser_semibin2.add_argument("--semibin2_engine", type=str, choices={'auto', 'cpu', 'gpu'}, default="auto", help="SemiBin2 | Device used to train the model [Default: auto]")
    parser_semibin2.add_argument("--semibin2_options", type=str, default="", help="SemiBin2 | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/BigDataBiology/SemiBin")

    # MetaDecoder
    parser_metadecoder = parser.add_argument_group('MetaDecoder arguments')
    parser_metadecoder.add_argument("--metadecoder_coverage_options", type=str, default="", help="MetaDecoder | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/jolespin/metadecoder-nal")
    parser_metadecoder.add_argument("--metadecoder_cluster_options", type=str, default="", help="MetaDecoder | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/jolespin/metadecoder-nal")

    # MetaCoAG
    parser_metacoag = parser.add_argument_group('MetaCoAG arguments')
    parser_metacoag.add_argument("--metacoag_assembler", type=str, choices={"auto", "spades", "megahit", "flye"}, default="auto", help="MeteaCoAG | Assembler used during assembly [Required if MetaCoAG is used]")
    parser_metacoag.add_argument("--metacoag_graph", default="auto", type=str, help="MetaCoAG | de Bruijn graph from SPAdes, MEGAHIT, or metaFlye [Required if MetaCoAG is used, if `auto` then assembly graphs will be looked]")
    parser_metacoag.add_argument("--metacoag_paths", default="auto", type=str, help="MetaCoAG | de Bruijn graph paths from SPAdes or metaFlye [Required if MetaCoAG is used with SPAdes or metaFlye]")
    parser_metacoag.add_argument("--metacoag_options", type=str, default="", help="MetaCoAG | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/jolespin/metacoag-nal")

    # Gene models
    parser_genemodels = parser.add_argument_group('Gene model arguments')
    parser_genemodels.add_argument("--pyrodigal_minimum_gene_length", type=int, default=90, help="Pyrodigal | Minimum gene length [Default: 90]")
    parser_genemodels.add_argument("--pyrodigal_minimum_edge_gene_length", type=int, default=60, help="Pyrodigal | Minimum edge gene length [Default: 60]")
    parser_genemodels.add_argument("--pyrodigal_maximum_gene_overlap_length", type=int, default=60, help="Pyrodigal | Maximum gene overlap length [Default: 60]")
    parser_genemodels.add_argument("--pyrodigal_genetic_code", type=int, default=11, help="Pyrodigal -g translation table [Default: 11]")

    parser_binrefinement = parser.add_argument_group('Bin refinement arguments')
    parser_binrefinement.add_argument("--binnette_contamination_weight", type=float, default=2.0, help="Binette | Contamination weight [Default: 2]")
    parser_binrefinement.add_argument("--binnette_options", type=str, default="", help="Binette | More options (e.g. --arg 1 ) [Default: 'https://github.com/genotoul-bioinfo/Binette']")
    
    parser_quality = parser.add_argument_group('Quality assessment arguments')
    parser_quality.add_argument("--checkm2_completeness", type=float, default=50.0, help="CheckM2 completeness threshold [Default: 50.0]")
    parser_quality.add_argument("--checkm2_contamination", type=float, default=10.0, help="CheckM2 contamination threshold [Default: 10.0]")
    parser_quality.add_argument("--checkm2_low_memory",action='store_true',help="CheckM2 low-memory mode")

    # parser_evaluation.add_argument("--checkm2_options", type=str, default="", help="CheckM lineage_wf | More options (e.g. --arg 1 ) [Default: '']")

    # rRNA
    parser_barrnap = parser.add_argument_group('barrnap arguments')
    parser_barrnap.add_argument("--barrnap_length_cutoff", type=float, default=0.8,  help="barrnap | Proportional length threshold to label as partial [Default: 0.8]")
    parser_barrnap.add_argument("--barrnap_reject", type=float, default=0.25,  help="barrnap | Proportional length threshold to reject prediction [Default: 0.25]")
    parser_barrnap.add_argument("--barrnap_evalue", type=float, default=1e-6,  help="barrnap | Similarity e-value cut-off [Default: 1e-6]")

    # tRNA
    parser_trnascan = parser.add_argument_group('tRNAscan-SE arguments')
    parser_trnascan.add_argument("--trnascan_options", type=str, default="", help="tRNAscan-SE | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/UCSC-LoweLab/tRNAscan-SE")

    # featureCounts
    parser_featurecounts = parser.add_argument_group('featureCounts arguments')
    parser_featurecounts.add_argument("--featurecounts_options", type=str, default="", help="featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/")

    # Tiara
    parser_domain = parser.add_argument_group('Domain classification arguments')
    parser_domain.add_argument("--logit_transform", type=str, default="softmax", help = " Transformation for consensus_domain_classification: {softmax, tss} [Default: softmax]")
    parser_domain.add_argument("--tiara_minimum_length", type=int, default=3000, help="Tiara | Minimum contig length. Anything lower than 3000 is not recommended. [Default: 3000]")
    parser_domain.add_argument("--tiara_options", type=str, default="", help="Tiara | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/ibe-uw/tiara")


    # Options
    opts = parser.parse_args()

    # Directories
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

    directories["project"] = create_directory(opts.project_directory)
    directories["sample"] = create_directory(os.path.join(directories["project"], opts.name))
    directories["output"] = create_directory(os.path.join(directories["sample"], "output"))
    directories["log"] = create_directory(os.path.join(directories["sample"], "log"))
    if not opts.tmpdir:
        opts.tmpdir = os.path.join(directories["sample"], "tmp")
    directories["tmp"] = create_directory(opts.tmpdir)
    directories["checkpoints"] = create_directory(os.path.join(directories["sample"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["sample"], "intermediate"))

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
    if "TMPDIR" in os.environ: print(os.environ["TMPDIR"], file=sys.stdout)
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
