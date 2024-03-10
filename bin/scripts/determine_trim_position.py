#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse
import numpy as np
import pandas as pd
from collections import OrderedDict

pd.options.display.max_colwidth = 100
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.8.11"

def statistics(fp, retain_index):
    # Build table
    quality_table = list()
    failed_reads = list()
    
    for read in tqdm(pyfastx.Fastq(fp), desc="Reading fastq filepath: {}".format(fp), unit=" read"):
        try:
            quality = np.asarray(read.quali).astype(int)
            quality_table.append(quality)
        except UnicodeDecodeError as e:
            failed_reads.append(read.name)
    if not retain_index:
        os.remove("{}.fxi".format(fp))

    quality_table = np.stack(quality_table)

    df_minmeanmax = pd.DataFrame([
        pd.Series(np.min(quality_table, axis=0), name="min"),
        pd.Series(np.mean(quality_table, axis=0), name="mean"),
        pd.Series(np.max(quality_table, axis=0), name="max"),
    ]).T
    df_quantiles = pd.DataFrame(np.quantile(quality_table, q=[0.25,0.5,0.75], axis=0).T, columns = ["q=0.25", "q=0.5", "q=0.75"])
    df_output = pd.concat([df_minmeanmax, df_quantiles], axis=1)
    df_output.index = df_output.index.values + 1
    df_output.index.name = "position"

    return df_output

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = """
    Single fastq file and output using -o argument: {} -o output.tsv file.fastq.gz\n
    Multiple fastq files, only include basename, and only include a certain field per file: {} -b -f q=0.25 -o output.tsv file_1.fq file_2.fq,...file_n.fq]\n
    Multiple fastq files and output to stdout: {}  file_1.fastq file_2.fq,...file_n.fq.gz > output.tsv
    """.format(__program__, __program__, __program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i", "--fastq_statistics",  type=str, default="stdin", help = "path/to/fastq_statistics.tsv. Rows=Positions (starting at 1), Columns=Samples. [Default: stdin])")
    parser.add_argument("-o","--output_directory", required=False, type=str, help = "Output directory [Required if you did not select --stdout]")
    parser.add_argument("-q","--minimum_quality",  default=30.0,type=float, help = "Minimum quality value")
    parser.add_argument("-m","--minimum_length",  default=100,type=int, help = "Minimum length.  If minimum quality value makes length shorter than this then an error will yield with which samples are responsible [Default: 100]")
    parser.add_argument("-w", "--window_size",  default=4,type=int, help = "Window size [Default: 4]")
    parser.add_argument("-l", "--maximum_average_loss",  type=float, help = "Maximum average loss for window size [Default: --window_size]")
    parser.add_argument("--figure_width",  default=55,type=float, help = "Figure width [Default: 34]")
    parser.add_argument("--plot_quality_label",  default="Quality [q=0.25]", type=str, help = "Quality metric label for plots [Default: Quality [q=0.25]]")
    parser.add_argument("--plot_title",  default="--fastq_statistics", type=str, help = "Title for plots [Default: Filename of --fastq_statistics]")
    parser.add_argument("-s", "--stdout",  action="store_true",  help = "Don't include any plots and just output the position to stdout")

    # parser.add_argument("--no_plots",  action="store_true", help = "Don't produce any plots")
    # parser.add_argument("-s","--sep",  default="\t",type=str, help = "Separator for output values")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Setup
    figsize = (opts.figure_width, opts.figure_width*0.16180339887)

    if opts.maximum_average_loss is None:
        opts.maximum_average_loss = opts.window_size
        
    if opts.fastq_statistics == "stdin":
        opts.fastq_statistics = sys.stdin
        if opts.plot_title == "--fastq_statistics":
            opts.plot_title = "Fastq Statistics [stdin]"
    else:
        if opts.plot_title == "--fastq_statistics":
            opts.plot_title = opts.fastq_statistics
        
    print("* Reading in fastq statistics: {}".format(opts.fastq_statistics), file=sys.stderr)
    df_fastq_statistics = pd.read_csv(opts.fastq_statistics, sep="\t", index_col=0)
    n_positions, n_samples = df_fastq_statistics.shape

    # If no --stdout, then create directory and load matplotlib
    if not opts.stdout:
        import matplotlib.pyplot as plt
        assert opts.output_directory is not None, "If --stdout is not selected then you must specify output directory via -o "
        os.makedirs(opts.output_directory, exist_ok=True)


    # Assertions
    print("* Checking format of fastq statistics", file=sys.stderr)
    assert all(df_fastq_statistics.index.map(lambda x: isinstance(x,int))), "--fastq_statistics index must be integers"
    assert all(np.arange(1, n_positions + 1) == df_fastq_statistics.index), "--fastq_statistics index must be [1,2,...,n] where n is the number of positions"

    # Number of samples above quality threshold per position
    df_boolean = df_fastq_statistics >= opts.minimum_quality
    number_of_samples_above_quality_threshold_per_position = df_boolean.sum(axis=1)
    number_of_samples_above_quality_threshold_per_position.index = number_of_samples_above_quality_threshold_per_position.index.astype(int)

    # Calculate loss per window
    window_to_loss = OrderedDict()
    for i in range(0, n_positions - opts.window_size + 1, 1):
        window = number_of_samples_above_quality_threshold_per_position.iloc[i:i + opts.window_size]
        loss = [n_samples]*opts.window_size - window
        window_to_loss[(window.index[0], window.index[-1])] = loss.mean()
    window_to_loss = pd.Series(window_to_loss)
    window_to_loss.index.names = ["start", "end"]

    start_trim_window, end_trim_window = window_to_loss[window_to_loss >= opts.maximum_average_loss].index[0]
    trim_window = number_of_samples_above_quality_threshold_per_position.loc[start_trim_window:end_trim_window]
    position_to_trim = trim_window.index[np.where(trim_window < n_samples)[0]].values[0]

    # Figures
    if not opts.stdout:
        # df_boolean.to_csv(os.path.join(opts.output_directory, "quality_distributions.boolean.tsv"), sep="\t")
        number_of_samples_above_quality_threshold_per_position.to_frame("number_of_samples_above_quality_threshold").to_csv(os.path.join(opts.output_directory, "number_of_samples_above_quality_threshold_per_position.tsv"), sep="\t")
        window_to_loss.to_frame("average_loss").to_csv(os.path.join(opts.output_directory, "average_loss_per_window.tsv"), sep="\t")

        print("* Plotting distributions of fastq statistics", file=sys.stderr)
        with plt.style.context("seaborn-white"):
            fig, axes = plt.subplots(figsize=figsize, nrows=2, sharex=True)
            ax = axes[0]
            df_fastq_statistics.T.plot(kind="box", ax=ax, color="black")
            ax.axhline(opts.minimum_quality, color="red", label="Minimum quality: {}".format(opts.minimum_quality))
            ax.axvline(position_to_trim, color="red", linewidth=1.618, label="Suggested trim: {}".format(position_to_trim))

            # ax.set_xticklabels(ax.get_xticklabels(), fontsize=12, rotation=90)
            # ax.set_xlabel("Position on read", fontsize=15)
            ax.set_ylabel(opts.plot_quality_label, fontsize=15)
            if opts.plot_title:
                ax.set_title(opts.plot_title, fontsize=15, fontweight="bold")

            print("* Plotting number of samples above quality thresholds for each position", file=sys.stderr)
            ax = axes[1]
            number_of_samples_above_quality_threshold_per_position.plot(marker="o", alpha=0.618,  color="black", label="N", ax=ax)
            for pos, v in number_of_samples_above_quality_threshold_per_position.items():
                if v < n_samples:
                    axes[0].axvline(pos, color="red", linewidth=0.1618, linestyle="--" )
                    axes[1].axvline(pos, color="red", linewidth=0.1618, linestyle="--", )

            ax.axhline(n_samples, color="black", label="Total samples: {}".format(n_samples))
            ax.axvline(position_to_trim, color="red", linewidth=1.618, label="Suggested trim: {}".format(position_to_trim))

            ax.set_xticks(range(1, n_positions+1))
            ax.set_xticklabels(ax.get_xticks(), fontsize=12, rotation=90)
            ax.set_xlabel("Position on read", fontsize=15)
            ax.set_ylabel("Number of samples\nabove quality threshold", fontsize=15)

            ax.set_xlim((1-0.5, n_positions+0.5))
            # ax.legend()
            # fig.savefig(os.path.join(opts.output_directory, "number_of_samples_above_quality_threshold_per_position.pdf"), dpi=300, format="pdf", bbox_inches="tight")

            ax.legend()
            fig.savefig(os.path.join(opts.output_directory, "quality_distributions.pdf"), dpi=300, format="pdf", bbox_inches="tight")


    # Position to trim 
    if position_to_trim < opts.minimum_length:
        print("! Position to trim (pos={}) is less than --minimum_length {}".format(position_to_trim, opts.minimum_length), file=sys.stderr)
        if not opts.stdout:
            print("! Review the low quality samples listed in the following file: {}".format(os.path.join(opts.output_directory, "low_quality_samples.tsv")), file=sys.stderr)
            df_fastq_statistics.loc[position_to_trim][lambda x: x < opts.minimum_quality].sort_values().to_frame("position_{}".format(position_to_trim)).to_csv(os.path.join(opts.output_directory, "low_quality_samples.tsv"), sep="\t")
        else:
            print("! Review the low quality samples listed below:", file=sys.stderr)
            df_fastq_statistics.loc[position_to_trim][lambda x: x < opts.minimum_quality].sort_values().to_frame("position_{}".format(position_to_trim)).to_csv(sys.stderr, sep="\t")

        sys.exit(1)

    print("* Position to trim: {}".format(position_to_trim), file=sys.stderr)
    if not opts.stdout:
        with open(os.path.join(opts.output_directory,"position_to_trim.txt"), "w") as f:
            print(position_to_trim, file=f)
    else:
        print(position_to_trim, file=sys.stdout)
    

if __name__ == "__main__":
    main()
