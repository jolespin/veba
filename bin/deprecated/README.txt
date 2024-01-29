* binning-prokaryotic.py - Deprecated on 2023.1.24.  This version uses CheckM (v1) which does not have direct support for CPR.  To get around this, VEBA ran CheckM and then GTDBTk and then CheckM once more.  Updated binning-prokaryotic.py uses CheckM2 which has native support for CPR.

* preprocess-kneaddata.py - Deprecated 2021.12.1.  This version uses Kneaddata while the later versions use fastq_preprocessor.

* cluster.py - Deprecated 2023.1.26.  This version uses the old clustering methodology based on OrthoFinder.

* classify-eukaryotic.py - Deprecated 2023.2.1.  This version can only work on genomes that have been binned throuhg VEBA.  The new version can handle de novo classifications.