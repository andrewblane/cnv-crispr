**CNV Dataset used**
Searching the [GEO](https://www.ncbi.nlm.nih.gov/geo/), a few suitable A375 CNV
datasets were found. Of three candidate datasets, one was generated using a
17.5k feature array (GSM218051), one using a 105k feature array (GSM935710) and,
the one I chose, using a 2x 400k feature array [GSM952580]
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM952580). I reasoned that
having more features would increase the mappability (coverage and accuracy) of
the CNV dataset onto the supplied guide-seq dataset.

**Requirements:**
python3
  - Will not run in python 2.x

python-tk (apt-get install python-tk)
  - Generates a warning when loaded that can be safely ignored on Python 3.

*Packages:*
pandas
numpy
matplotlib
numpy
scipy


**Instructions:**
1. Set the three input paths in guide-seq_cnv_analysis.py:
    These refer to the supplied guide-seq data, a chosen CNV data source and
    a supplied Affymetrix array key that relates numeric gene ID keys in the
    CNV data to gene names. The defaults work with the pull as submitted.

2. I tried to minimize exotic module requirements. Provided all are in place
    and installed correctly, the script can be run as a regular Python script:
    python3 guide-seq_cnv_analysis.py

**Expected output:**
    a. A CSV file showing CNV at each gene in the guide-seq dataset. A subset
    of these are also printed to stdout when the script is run.
    b. A PDF with a 2d histogram of the relationship between CNV and guide-seq
    abundance. The strategy for calculating this is described in more detail as
    script comments. A linear fit is plotted through the data and statistics
    shown on the plot.

**Summary of findings:**
The data show that CNV affects guide selection, with regions of higher copy
number associated with lessened guide selection. This relationship is weakened
at later timepoints, presumably because other selective pressures come to
dominate by then.
