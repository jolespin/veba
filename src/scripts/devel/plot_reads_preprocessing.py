#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, re
from collections import defaultdict, OrderedDict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.8.8"

# Parse Basename
def parse_basename(query: str, naming_scheme: str): # "[ID]_R[DIRECTION]_001.fastq.gz" 
    """
    Adapted from the following source: 
    * @s-b
    * https://stackoverflow.com/questions/73296040/how-to-parse-out-two-fields-out-of-string-using-regex-in-python
    """
    def repl(match_object):
        inside_bracket = match_object.group(1)
        if inside_bracket == "DIRECTION":
            return r"(?P<DIRECTION>[12])"
        if inside_bracket == "ID":
            return r"(?P<ID>[-.\w]+)"
    pattern = re.sub(r"\[(.*?)\]", repl, naming_scheme)
    match = re.match(pattern, query)
    return match["ID"], match["DIRECTION"]
    

def concatenate_seqkit_stats_dataframes(filepaths):
    """
    Input: Path(s) to at least one TSV file.
    Output: Merged input TSV files as a DataFrame.
    """
    
    # Check if at least one TSV file is provided
    if not filepaths:
        raise Excecption("At least 1 seqkit stats TSV file is required")
    # assert len(filepaths) > 0, "At least 1 seqkit stats TSV file is required"

    # Initialize an empty list to store dataframes
    dataframes = []

    for fp in filepaths:
        df = pd.read_csv(fp, sep='\t', index_col=0)
        assert "num_seqs" in df.columns, "Please ensure that input filepaths are from `seqkit stats -a -T`"
        dataframes.append(df)
        # except (FileNotFoundError: # This will be handled by pd.read_csv
        
    # Concatenate all dataframes into a single dataframe
    df_concat = pd.concat(dataframes, axis=0) # This will catch errors on empty lists where there's no dataframes to concatenate.  Also good to be explicit w/ which axis you are concatenating along for readability  

    return df_concat

# Format seqkit stats output
def format_read_preprocessing(
    number_of_sequences:pd.Series,
    naming_scheme,
    sort_by="samples",
    ascending=True,
    ):
    assert sort_by in {"total", "trimmed", "discarded", "samples"}

    output = defaultdict(dict)

    # Format read prepocessing
    for fp, number_of_reads in number_of_sequences.items():
        # From preprocessing module (i.e., fastq_preprocessor)
        if "intermediate" in fp:
            if "_1.fastq.gz" in fp:
                id_sample = fp.split("/")[-4]
                fn = fp.split("/")[-1].split("_1.fastq.gz")[0]
                output[id_sample][fn] = number_of_reads
        else:
            # Unprocessed reads
            basename = fp.split("/")[-1]
            id_sample, direction = parse_basename(basename, naming_scheme=naming_scheme)
            if str(direction) == "1":
                output[id_sample]["unprocessed"] = number_of_reads
                
    df_read_summary = pd.DataFrame(output).T
    df_read_summary.index.name = "id_sample"

    # Add discarded reads
    df_read_summary["discarded"] = df_read_summary["unprocessed"] - df_read_summary["trimmed"]

    if sort_by == "samples":
        df_read_summary = df_read_summary.sort_index(ascending=ascending)
    if sort_by == "total":
        df_read_summary = df_read_summary.sort_values("unprocessed", ascending=ascending)
    if sort_by == "trimmed":
        df_read_summary = df_read_summary.sort_values("trimmed", ascending=ascending)
    if sort_by == "discarded":
        df_read_summary = df_read_summary.sort_values("discarded", ascending=ascending)

    return df_read_summary


# (base) [jespinoz@login01 TestVEBA]$ cat veba_output/preprocess/S1/output/seqkit_stats.concatenated.tsv
# file	format	type	num_seqs	sum_len	min_len	avg_len	max_len
# veba_output/preprocess/S1/intermediate/3__bbduk/kmer_hits_1.fastq.gz	FASTQ	DNA	2903	436935	81	150.5	151
# veba_output/preprocess/S1/intermediate/3__bbduk/kmer_hits_2.fastq.gz	FASTQ	DNA	2903	436849	81	150.5	151
# veba_output/preprocess/S1/intermediate/3__bbduk/non-kmer_hits_1.fastq.gz	FASTQ	DNA	1414483	212887364	75	150.5	151
# veba_output/preprocess/S1/intermediate/3__bbduk/non-kmer_hits_2.fastq.gz	FASTQ	DNA	1414483	212888714	75	150.5	151
# Fastq/S1_1.fastq.gz	FASTQ	DNA	1474203	221880932	75	150.5	151
# Fastq/S1_2.fastq.gz	FASTQ	DNA	1460132	219765860	75	150.5	151
# veba_output/preprocess/S1/intermediate/1__fastp/trimmed_1.fastq.gz	FASTQ	DNA	1417386	213324299	75	150.5	151
# veba_output/preprocess/S1/intermediate/1__fastp/trimmed_2.fastq.gz	FASTQ	DNA	1417386	213325563	75	150.5	151
# veba_output/preprocess/S1/intermediate/1__fastp/trimmed_singletons.fastq.gz	FASTQ	DNA	84846	12748323	75	150.3	151
# veba_output/preprocess/S1/intermediate/2__bowtie2/cleaned_1.fastq.gz	FASTQ	DNA	1417386	213324299	75	150.5	151
# veba_output/preprocess/S1/intermediate/2__bowtie2/cleaned_2.fastq.gz	FASTQ	DNA	1417386	213325563	75	150.5	151
# veba_output/preprocess/S1/intermediate/2__bowtie2/cleaned_singletons.fastq.gz			0	0	0	0.0	0
# veba_output/preprocess/S1/intermediate/2__bowtie2/contaminated_1.fastq.gz			0	0	0	0.0	0
# veba_output/preprocess/S1/intermediate/2__bowtie2/contaminated_2.fastq.gz			0	0	0	0.0	0
# veba_output/preprocess/S1/intermediate/2__bowtie2/contaminated_singletons.fastq.gz			0	0	0	0.0	0



def format_colors(ran_bowtie2=False, ran_bbduk=False, palette="pastel", reverse_color_order=False):
    """
    Scheme: 
    Trimmed + Discarded = Unprocessed
    Cleaned + Contaminated = Trimmed
    Kmer hits + Non-kmer hits = Cleaned
    """
    
    labels = None

    colors = sns.color_palette(palette=palette, n_colors=4)
    
    if reverse_color_order:
        colors = colors[::-1]

    label_to_color = pd.Series(dict(zip(["discarded", "contaminated", "kmer_hits", "non-kmer_hits"], colors)))

    if not any([
        ran_bowtie2,
        ran_bbduk,
        ]):
        labels = ["discarded", "trimmed"]

    if all([
        ran_bowtie2,
        not ran_bbduk,
        ]):
        labels = ["discarded", "contaminated", "cleaned"]


    if all([
        not ran_bowtie2,
        ran_bbduk,
        ]):
        labels = ["discarded", "kmer_hits", "non-kmer_hits"]


    if all([
        ran_bowtie2,
        ran_bbduk,
        ]):
        labels = ["discarded", "contaminated", "kmer_hits", "non-kmer_hits"]

    assert labels is not None

    return label_to_color[labels].to_dict(into=OrderedDict)

def format_proportions(df:pd.DataFrame, percentage=False, ran_bowtie2=False, ran_bbduk=False):
    """
    Scheme: 
    Trimmed + Discarded = Unprocessed
    Cleaned + Contaminated = Trimmed
    Kmer hits + Non-kmer hits = Cleaned
    """
    assert set(df.columns) >= {"discarded", "trimmed", "unprocessed"}, "Please ensure the `seqkit stats -a -T` tables provided are from either VEBA's preprocessing module or fastq_preprocessor.  User provided the following columns in DataFrame: {}".format(", ".join(df.columns))

    df_proportions = pd.DataFrame()

    if not any([
        ran_bowtie2,
        ran_bbduk,
        ]):
        df_proportions["trimmed"] = df["trimmed"]/df["unprocessed"]
        df_proportions["discarded"] = df["discarded"]/df["unprocessed"]

    if all([
        ran_bowtie2,
        not ran_bbduk,
        ]):
        df_proportions["cleaned"] = df["cleaned"]/df["trimmed"] 
        df_proportions["contaminated"] = df["contaminated"]/df["trimmed"] 
        df_proportions["discarded"] = df["contaminated"]/df["unprocessed"] 

    if all([
        not ran_bowtie2,
        ran_bbduk,
        ]):
        df_proportions["discarded"] = df["discarded"]/df["unprocessed"] 
        df_proportions["kmer_hits"] = df["kmer_hits"]/df["trimmed"] 
        df_proportions["non-kmer_hits"] = df["non-kmer_hits"]/df["trimmed"] 

    if all([
        ran_bowtie2,
        ran_bbduk,
        ]):
        df_proportions["discarded"] = df["discarded"]/df["unprocessed"] 
        df_proportions["kmer_hits"] = df["kmer_hits"]/df["cleaned"] 
        df_proportions["non-kmer_hits"] = df["non-kmer_hits"]/df["cleaned"] 
        df_proportions["contaminated"] = df["contaminated"]/df["trimmed"] 

    assert not df_proportions.empty, "Please ensure the `seqkit stats -a -T` tables provided are from either VEBA's preprocessing module or fastq_preprocessor.  User provided the following columns in DataFrame: {}".format(", ".join(df.columns))

    if percentage:
        df_proportions = df_proportions * 100
        total = df_proportions.sum(axis=1)
        assert np.allclose(total.values, np.asarray([100]*total).size), "Percentages don't add up to 100%"
    else:
        total = df_proportions.sum(axis=1)
        assert np.allclose(total.values, np.asarray([1.0]*total).size), "Proportions don't add up to 1.0"

    return df_proportions

def plot_read_preprocessing(
    df:pd.DataFrame,
    percentage=False,
    kmer_hits_relabel=None,
    nonkmer_hits_relabel="other",
    
    palette="pastel",
    border_color="white",
    reverse_color_order=False, 
    figsize=(13,5),
    title=None,
    style="seaborn-white",

    show_xgrid=True,
    show_ygrid=True,
    show_legend=True,
    ylabel="Number of reads",
    xlabel="Samples",
    legend_kws=dict(),
    legend_title=None,
    
    title_kws=dict(),
    barchart_kws=dict(),
    legend_title_kws=dict(),
    axis_label_kws=dict(),
    fig_kws=dict(),
    
    ):
    """
    Input: A DataFrame that is meant to be plotted, as well as a version of it containing proportional data.
    Output: 3 barcharts of the data in the input DataFrame.

    Scheme: 
    Trimmed + Discarded = Unprocessed
    Cleaned + Contaminated = Trimmed
    Kmer hits + Non-kmer hits = Cleaned
    """

    _fig_kws = {"figsize":figsize}
    _fig_kws.update(fig_kws)
    _title_kws = {"fontsize":14, "fontweight":"bold"}
    _title_kws.update(title_kws)
    _barchart_kws ={"edgecolor":border_color}
    _barchart_kws.update(barchart_kws)
    _legend_kws = {'fontsize': 10}#, 'loc': 'center left', 'bbox_to_anchor': (1, 0.5)}
    _legend_kws.update(legend_kws)
    _legend_title_kws = {"size":12, "weight":"bold"}
    _legend_title_kws.update(legend_title_kws)
    _axis_label_kws = {"fontsize":14}
    _axis_label_kws.update(axis_label_kws)
        
    # Checks
    ran_bowtie2 = "contaminated" in df.columns #True/False
    ran_bbduk = all([                          #True/False
        "kmer_hits" in df.columns,
        "non-kmer_hits" in df.columns,
    ])

    # Make a copy so if and/or when we relabel the dataframe it doesn't change the original outside of the function
    df = df.copy()

    # Get a color label dictionary
    label_to_color = format_colors(ran_bowtie2=ran_bowtie2, ran_bbduk=ran_bbduk, palette=palette, reverse_color_order=reverse_color_order)

    # Get proportions
    df_proportions = format_proportions(df, percentage=percentage, ran_bowtie2=ran_bowtie2, ran_bbduk=ran_bbduk)

    if kmer_hits_relabel is not None:
        df.columns = df.columns.map(lambda column_name: kmer_hits_relabel if column_name == "kmer_hits" else column_name)
        df_proportions.columns = df_proportions.columns.map(lambda column_name: kmer_hits_relabel if column_name == "kmer_hits" else column_name)

    if nonkmer_hits_relabel is not None:
        df.columns = df.columns.map(lambda column_name: nonkmer_hits_relabel if column_name == "non-kmer_hits" else column_name)
        df_proportions.columns = df_proportions.columns.map(lambda column_name: nonkmer_hits_relabel if column_name == "non-kmer_hits" else column_name)

    # Create the figure/axes for subplots and set the color palette
    with plt.style.context(style):
        fig, axes = plt.subplots(**_fig_kws, nrows=3)
        
        

            

    
        # # Plot the first subplot (counts data)
        # barchart_df.plot(kind='bar', stacked=True, ax=axes[0])
        # axes[0].set_title('Preprocessing Summary', fontsize=15)
        # axes[0].set_ylabel('Counts', fontsize=13)
        # axes[0].set_xticklabels([])
        
        # # Plot the second subplot (log transformed counts data)
        # barchart_df.plot(kind='bar', stacked=True, ax=axes[1])
        # axes[1].set_ylabel('Log Counts', fontsize=13)
        # axes[1].set_xticklabels([])
        # axes[1].set_yscale(LogScale(axis=0,base=10))
        # axes[1].get_legend().remove()
        
        # # Plot the third subplot (proportions)
        # barchart_df_proportions.plot(kind='bar', stacked=True, ax=axes[2])
        # axes[2].set_ylabel('Proportions', fontsize=13)
        # axes[2].get_legend().remove()
        
        # # Remove the X-axis labels
        # for ax in axes:
        #     ax.set_xlabel('')
        
        # # Adjust the spacing between subplots, legend, and X-label
        # plt.subplots_adjust(hspace=0.1)
        # axes[0].legend(title='Categories', loc='upper left', bbox_to_anchor=(1, 1))
        # plt.xlabel('Samples', fontsize=13)
        
    
        
        return fig, axes, df_proportions
    
#
def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <index.list> -t <table.tsv> -o <output.tsv>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", default="stdin", nargs="+", type=str, help = "path/to/seqkit_stats.tsv. Can be multiple files or one file [Default: stdin]")
    parser.add_argument("-o","--output", required=True, type=str, help = "path/to/figure.[pdf|png|svg]")
    parser.add_argument("-c","--color_palette", default="pastel", help = "Color palette [Default: pastel]")
    parser.add_argument("-s","--figsize", type=str, default="13,5", help = "Comma separated string width,height [Default: 13,5]")
    parser.add_argument("-t","--tick_fontsize", type=str, default="12", help = "xtick fontsize [Default: 12]")
    parser.add_argument("-l","--label_fontsize", type=str, default="14", help = "label fontsize [Default: 14]")
    parser.add_argument("-n","--naming_scheme", default="[ID]_R[DIRECTION]_001.fastq.gz", type=str, help = "Naming scheme.  Use [ID] for identifier name and [DIRECTION] for 1 or 2. [Default: [ID]_R[DIRECTION]_001.fastq.gz]")
    parser.add_argument("-p","--percentage", action = "store_true", help = "Use percentages instead of proportions")            # Set the default naming scheme and help message
    parser.add_argument("--sort_by", type=str, default="samples", help = "Sort samples by {total, trimmed, discarded, samples} [Default: samples]")
    # parser.add_argument("-m", "--mode", type=str, default="short", help = "Mode for fastq {short, long} [Default: short]") # Throw an error if long reads are selected.  Say it will be implemented in the future.

    # Options
    opts = parser.parse_args()

    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # I/O
    if opts.input == "stdin":
        opts.input = sys.stdin 
        df_seqkit_stats = pd.read_csv(opts.input, sep="\t", index_col=0)
        assert "num_seqs" in df_seqkit_stats.columns, "Please ensure that input is from `seqkit stats -a -T`"
    else:
        df_seqkit_stats = concatenate_seqkit_stats_dataframes(opts.input)

    assert df_seqkit_stats.index.nunique() ==  df_seqkit_stats.shape[0]

    # For read preprocessing table
    df_read_summary = format_read_preprocessing(df_seqkit_stats["num_seqs"], naming_scheme=opts.naming_scheme, sort_by=opts.sort_by)
 
    fig, ax, df_proportions = plot_read_preprocessing(
        df_read_summary,
        percentage=opts.percentage,
        # The rest of the kwargs
        )
    fig.savefig(opts.output, dpi=300, bbox_inches="tight")

    # Output table to stdout
    df_read_summary.columns = df_read_summary.columns.map(lambda x: ("counts", x))
    if opts.percentage:
        df_proportions.columns = df_proportions.columns.map(lambda x: ("percentages", x))
    else:
        df_proportions.columns = df_proportions.columns.map(lambda x: ("proportions", x))

    df_output = pd.concat([df_read_summary, df_proportions], axis=1)
    df_output.to_csv(sys.stdout, sep="\t")

if __name__ == "__main__":
    main()
