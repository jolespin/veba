### CRISPR-Cas system screening of *de novo* genomes
If you have prokaryotic sequences, it may be beneficial for you to screen for CRISPR-Cas systems as these may provide insight into evolutionary history of CRISPR-Cas systems and also may be useful in biotechnology.  This tutorial will show you how use `CRISPRCasTyper` as a post hoc analysis for screening genomes.

What you'll end up with at the end of this is are CRISPR-Cas system candidates and spacer seqeunces.

Please refer to the [end-to-end metagenomics](end-to-end_metagenomics.md) or [recovering viruses from metatranscriptomics](recovering_viruses_from_metatranscriptomics.md) workflows for details on binning, clustering, and annotation.

_____________________________________________________

#### Steps:
0. Create a `conda` environment for `CRISPRCasTyper`
1. Screen using `CRISPRCasTyper`
2. Interpret the results

_______________________________________________________

#### 0. Create and activate a `conda` environment for `CRISPRCasTyper`

```
mamba create -n cctyper_env -c conda-forge -c bioconda -c russel88 cctyper -y

conda activate cctyper_env
```

#### 1. Screen using `CRISPRCasTyper`

At this point, it's assumed you have the following: 

* Genome assemblies.  These can either be MAGs binned with VEBA, binned elsewhere, or even reference genomes you downloaded. 


For this dataset, we are going to use the all of the MAGs from the Siberian permafrost case study in the manuscript (BioProject: PRJNA596250). 


```
# Concatenate the assemblies
cat veba_output/binning/prokaryotic/*/output/genomes/*.fa > veba_output/misc/all_genomes.fa

# Run CCTyper
cctyper --seed 0 --prodigal meta veba_output/misc/all_genomes.fa veba_output/misc/cctyper_output
```

This generates several files with the most useful shown below: 

* cas_operons.tab - All certain Cas operons
* CRISPR_Cas.tab - CRISPR_Cas loci, with consensus subtype prediction
* crisprs.gff - CRISPR-Cas system features in GFF format
* plot.svg - Plot of CRISPR-Cas system features in SVG format
* spacers/*.fa - Fasta files with all spacer sequences

The full output file list and descriptions are available here: https://github.com/Russel88/CRISPRCasTyper
_____________________________________________________

#### Next steps:

Try to identify the source organisms of the spacer sequences or see if your CRISPR-associated proteins are unique.