### Phylogenetic inference
This walkthrough goes through phylogenetic inference of eukaryotic diatoms.

What you'll end up with at the end of this are phylogenetic trees for recovered eukaryotic diatoms.

It is assumed you've completed either the [end-to-end metagenomics](end-to-end_metagenomics.md) or [recovering viruses from metatranscriptomics](recovering_viruses_from_metatranscriptomics.md) walkthroughs.

_____________________________________________________

#### Steps:

1. Download the proteomes of similar organisms
2. Perform phylogenetic inference on proteomes

_____________________________________________________


#### 1. Download the proteomes of similar organisms

We are going to placing the following diatom genomes recovered from the *Plastisphere* dataset:

```
SRR17458614__METABAT2__E.1__bin.2
SRR17458615__METABAT2__E.1__bin.2
SRR17458638__METABAT2__E.1__bin.2
SRR17458638__METABAT2__E.1__bin.3
```

**You can find download them on FigShare**:

Espinoza, Josh (2022): Genome assemblies, gene models, gene annotations, taxonomy classifications, clusters, and counts tables. figshare. Dataset. https://doi.org/10.6084/m9.figshare.20263974.v1 


We also need a set of similar organisms to place the recovered diatoms.  For this, we are pulling from [*MMETSP*](https://zenodo.org/record/3247846/files/mmetsp_dib_trinity2.2.0_pep_zenodo.tar.gz) and [*NCBI*](https://www.ncbi.nlm.nih.gov/assembly/?term=Bacillariophyceae).

Once we have all of the proteomes downloaded, we need a list file that gives the path to where these proteomes are located.  I've copied over all of the proteomes into a directory called `diatoms`.  Once I did that I was able to get all of the paths easily: `ls diatoms > veba_output/misc/diatoms.proteins.list`

`diatoms.proteins.list` looks like this:

```
cat veba_output/misc/diatoms.proteins.list

diatoms/GCA_000149405.2.faa
diatoms/GCA_000150955.2.faa
diatoms/GCA_000296195.2.faa
diatoms/GCA_001750085.1.faa
diatoms/GCA_002217885.1.faa
diatoms/GCA_900660405.1.faa
diatoms/MMETSP0009.faa
diatoms/MMETSP0010.faa
diatoms/MMETSP0013.faa
diatoms/MMETSP0014.faa
diatoms/MMETSP0015.faa
diatoms/MMETSP0017.faa
diatoms/MMETSP0088.faa
diatoms/MMETSP0090.faa
diatoms/MMETSP0091.faa
diatoms/MMETSP0092.faa
diatoms/MMETSP0139.faa
diatoms/MMETSP0140.faa
diatoms/MMETSP0141.faa
diatoms/MMETSP0142.faa
diatoms/MMETSP0149.faa
diatoms/MMETSP0150.faa
diatoms/MMETSP0152.faa
diatoms/MMETSP0154.faa
diatoms/MMETSP0156.faa
diatoms/MMETSP0158.faa
diatoms/MMETSP0160.faa
diatoms/MMETSP0169.faa
diatoms/MMETSP0171.faa
diatoms/MMETSP0174.faa
diatoms/MMETSP0176.faa
diatoms/MMETSP0200.faa
diatoms/MMETSP0316.faa
diatoms/MMETSP0317.faa
diatoms/MMETSP0318.faa
diatoms/MMETSP0319.faa
diatoms/MMETSP0320.faa
diatoms/MMETSP0321.faa
diatoms/MMETSP0322.faa
diatoms/MMETSP0327.faa
diatoms/MMETSP0329.faa
diatoms/MMETSP0397.faa
diatoms/MMETSP0403.faa
diatoms/MMETSP0404.faa
diatoms/MMETSP0418.faa
diatoms/MMETSP0492.faa
diatoms/MMETSP0493.faa
diatoms/MMETSP0494.faa
diatoms/MMETSP0562.faa
diatoms/MMETSP0563.faa
diatoms/MMETSP0578.faa
diatoms/MMETSP0580.faa
diatoms/MMETSP0593.faa
diatoms/MMETSP0603.faa
diatoms/MMETSP0604.faa
diatoms/MMETSP0693.faa
diatoms/MMETSP0696.faa
diatoms/MMETSP0697.faa
diatoms/MMETSP0698.faa
diatoms/MMETSP0699.faa
diatoms/MMETSP0705.faa
diatoms/MMETSP0706.faa
diatoms/MMETSP0707.faa
diatoms/MMETSP0708.faa
diatoms/MMETSP0713.faa
diatoms/MMETSP0716.faa
diatoms/MMETSP0717.faa
diatoms/MMETSP0718.faa
diatoms/MMETSP0719.faa
diatoms/MMETSP0724.faa
diatoms/MMETSP0725.faa
diatoms/MMETSP0726.faa
diatoms/MMETSP0727.faa
diatoms/MMETSP0733.faa
diatoms/MMETSP0734.faa
diatoms/MMETSP0735.faa
diatoms/MMETSP0736.faa
diatoms/MMETSP0737.faa
diatoms/MMETSP0738.faa
diatoms/MMETSP0739.faa
diatoms/MMETSP0740.faa
diatoms/MMETSP0744.faa
diatoms/MMETSP0745.faa
diatoms/MMETSP0746.faa
diatoms/MMETSP0747.faa
diatoms/MMETSP0751.faa
diatoms/MMETSP0752.faa
diatoms/MMETSP0753.faa
diatoms/MMETSP0754.faa
diatoms/MMETSP0786.faa
diatoms/MMETSP0789.faa
diatoms/MMETSP0794.faa
diatoms/MMETSP0800.faa
diatoms/MMETSP0816.faa
diatoms/MMETSP0850.faa
diatoms/MMETSP0851.faa
diatoms/MMETSP0852.faa
diatoms/MMETSP0853.faa
diatoms/MMETSP0878.faa
diatoms/MMETSP0879.faa
diatoms/MMETSP0880.faa
diatoms/MMETSP0881.faa
diatoms/MMETSP0898.faa
diatoms/MMETSP0899.faa
diatoms/MMETSP0900.faa
diatoms/MMETSP0901.faa
diatoms/MMETSP0902.faa
diatoms/MMETSP0903.faa
diatoms/MMETSP0904.faa
diatoms/MMETSP0905.faa
diatoms/MMETSP0906.faa
diatoms/MMETSP0907.faa
diatoms/MMETSP0908.faa
diatoms/MMETSP0909.faa
diatoms/MMETSP0910.faa
diatoms/MMETSP0911.faa
diatoms/MMETSP0912.faa
diatoms/MMETSP0913.faa
diatoms/MMETSP0918.faa
diatoms/MMETSP0920.faa
diatoms/MMETSP0970.faa
diatoms/MMETSP0971.faa
diatoms/MMETSP0972.faa
diatoms/MMETSP0973.faa
diatoms/MMETSP0998.faa
diatoms/MMETSP1001.faa
diatoms/MMETSP1002.faa
diatoms/MMETSP1005.faa
diatoms/MMETSP1010.faa
diatoms/MMETSP1012.faa
diatoms/MMETSP1013.faa
diatoms/MMETSP1039.faa
diatoms/MMETSP1040.faa
diatoms/MMETSP1057.faa
diatoms/MMETSP1058.faa
diatoms/MMETSP1059.faa
diatoms/MMETSP1060.faa
diatoms/MMETSP1061.faa
diatoms/MMETSP1062.faa
diatoms/MMETSP1063.faa
diatoms/MMETSP1064.faa
diatoms/MMETSP1065.faa
diatoms/MMETSP1066.faa
diatoms/MMETSP1067.faa
diatoms/MMETSP1070.faa
diatoms/MMETSP1071.faa
diatoms/MMETSP1171.faa
diatoms/MMETSP1175.faa
diatoms/MMETSP1176.faa
diatoms/MMETSP1322.faa
diatoms/MMETSP1336.faa
diatoms/MMETSP1352.faa
diatoms/MMETSP1360.faa
diatoms/MMETSP1361.faa
diatoms/MMETSP1362.faa
diatoms/MMETSP1394.faa
diatoms/MMETSP1405.faa
diatoms/MMETSP1406.faa
diatoms/MMETSP1407.faa
diatoms/MMETSP1408.faa
diatoms/MMETSP1409.faa
diatoms/MMETSP1410.faa
diatoms/MMETSP1411.faa
diatoms/MMETSP1412.faa
diatoms/MMETSP1413.faa
diatoms/MMETSP1414.faa
diatoms/MMETSP1415.faa
diatoms/MMETSP1416.faa
diatoms/MMETSP1417.faa
diatoms/MMETSP1418.faa
diatoms/MMETSP1419.faa
diatoms/MMETSP1420.faa
diatoms/MMETSP1421.faa
diatoms/MMETSP1422.faa
diatoms/MMETSP1423.faa
diatoms/MMETSP1428.faa
diatoms/MMETSP1429.faa
diatoms/MMETSP1432.faa
diatoms/MMETSP1434.faa
diatoms/MMETSP1435.faa
diatoms/MMETSP1437.faa
diatoms/MMETSP1442.faa
diatoms/MMETSP1443.faa
diatoms/MMETSP1447.faa
diatoms/MMETSP1449.faa
diatoms/SRR17458614__METABAT2__E.1__bin.2.faa
diatoms/SRR17458615__METABAT2__E.1__bin.2.faa
diatoms/SRR17458638__METABAT2__E.1__bin.2.faa
diatoms/SRR17458638__METABAT2__E.1__bin.3.faa
```

#### 2. Perform phylogenetic inference on proteomes

Now that we have all of the files we need, we can perform phylogenetic inference using BUSCO's eukaryota_odb10 markers and score cutoffs.  For eukaryotes, it's advised that you use the eukaryota_odb10 marker set because this is the core marker set used for classification.  This isn't the case for prokaryotes and viruses.  If you don't have enough resources to run maximum likelihood trees via *IQTREE2* then use `--no_iqtree`.

**Conda Environment:** `conda activate VEBA-phylogeny_env`


```
# Set the number of threads
N_JOBS=16

N=phylogeny-diatoms
rm -f logs/${N}.*

# Proteome file paths
PROTEINS=veba_output/misc/diatoms.proteins.list

# Marker set
HMM=${VEBA_DATABASE}/MarkerSets/eukaryota_odb10.hmm
SCORES={VEBA_DATABASE}/MarkerSets/eukaryota_odb10.scores_cutoff.tsv.gz

# Minimum ratio of genomes include in alignment. 
# This removes markers that are under represented. Default is 0.95
MINIMUM_GENOMES_ALIGNED_RATIO=0.95

# Output directory
OUT_DIR=veba_output/phylogeny/diatoms

# Directory
CMD="source activate VEBA-phylogeny_env && phylogeny.py -a ${PROTEINS} -o ${OUT_DIR} -p ${N_JOBS}  -f name --no_iqtree -d ${HMM} -s ${SCORES} --minimum_genomes_aligned_ratio ${MINIMUM_GENOMES_ALIGNED_RATIO} 

# Either run this command or use SunGridEnginge/SLURM
```

The following output files will produced: 

* alignment_table.boolean.tsv.gz - Alignment table of (n = genomes, m = markers, ij=fasta pass qc)
* concatenated_alignment.fasta - Concatenated protein alignment of all marker hits
* concatenated_alignment.fasttree.nw - FastTree newick format based on concatenated alignment
* concatenated_alignment.fasttree.nw.pdf - PDF visualization of newick tree by ETE3
* prefiltered_alignment_table.tsv.gz - Prefiltered alignment table of (n = genomes, m = markers, ij=fasta alignment)
* output.treefile - IQTREE2 newick format based on concatenated alignment (if --no_iqtree is not selected)
* output.treefile.pdf - PDF visualization of newick tree by ETE3

#### Next steps:

Now it's time to visualize the tree.  We recommend either using [*ETE 3*](http://etetoolkit.org/docs/latest/tutorial/index.html) (either the web interface or Python) or [iTOL](https://itol.embl.de/) for a powerful web interface.