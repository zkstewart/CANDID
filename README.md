# CANDID
Cluster-based Alignment-free Novel Domain Iterative Discovery

CANDID is a program that is in its final stage of development and otherwise fully functional. Its intent is to discover novel domains by excluding known domain regions from proteins and subsequently utilising sensitive MMseqs2 search to find potential domain regions which are clustered based on alignment-free distance measures and iteratively used as HMMER queries until convergence is reached.

## Prerequisites
CANDID operates as a pipeline which calls quite a few external programs and resources, as well as packages within Python and R. Below is a list detailing all the things you'll need to successfully run CANDID.

### Mandatory
* If running on Windows, Cygwin is required; testing is necessary to identify the minimum required set of modules
* Python 3.X
  * biopython
  * numpy
* Python 2.7
* CD-Hit
* HMMER 3.1+ [Needs testing with 3.2 to ensure everything is alright]
* Pre-constructed HMM database of known domains [Put hmm_db_download.py script in this repository after updates]
* SignalP 4.1
* seg (ftp://ftp.ncbi.nih.gov/pub/seg/seg)
* psCoils (http://www.biocomp.unibo.it/piero/PS-COILS/download)
* MMseqs2 (https://github.com/soedinglab/MMseqs2)
* MAFFT
* R (version shouldn't matter so long as the below two packages are available for it)
  * msa (https://bioconductor.org/packages/release/bioc/html/msa.html)
  * odseq (https://bioconductor.org/packages/release/bioc/html/odseq.html)

### Optional / choices
* Python 3.X
  * If NOT running Hammock for clustering, alfpy is required https://github.com/aziele/alfpy
  * If NOT running Hammock for clustering, HDBSCAN is required https://github.com/scikit-learn-contrib/hdbscan
* Hammock can optionally replace ALF-HDBSCAN clustering (https://github.com/krejciadam/hammock)

Note that, while Hammock is offered as an option for clustering, from initial testing alignment-free HDBSCAN (ALF-HDBSCAN) clustering seems like a better option in most scenarios. It is a lot faster which is important for the iterative nature of CANDID to converge upon a solution quicker. Moreover, it is more likely to detect clusters in a data set where Hammock may struggle. If you have the time and computational resources, there's no reason why you shouldn't be able to try out both solutions and see what works best.

**IMPORTANT**: I haven't figured out how to get Hammock to run on Windows using Cygwin. If you can do this, let me know how. Otherwise, you may be limited to ALF-HDBSCAN on this operating system.

## Development progress

Changes that have occurred recently:

1. MSA curation functions are completely implemented. ODseq's outlier predictions require additional validation with my own (rough) outlier detection method. This combination helps to temper these results to be less strict in certain situations.
2. Optimal default parameters are likely to be minsize = (2 or 3) and minsample = 2. minsize depends on the amount of data that is input; with small datasets use 2, and with larger ones use 3. minsample is tricky; in most cases 2 works best but in at least one case with a smaller dataset only minsample == 1 was able to recover domains.
3. Testing code is being stripped out since we're largely done with developing individual functions.

Directions for improvement:

1. Improvement to CANDID is likely only possible now through the assistant function detailed below for combining outputs from multiple runs. The core functionality of the program is complete and I do not have any ideas or plans for adding new things in.
2. Do a bit more testing to see if I should change the convergence detection to stop in the first instance that no changes are detected.
3. Change how output files are presented/hide the processing files elsewhere.
3. Continued updating to this README. Include subheadings like 'How to use', 'Parameter selection', 'HMM database creation', 'CANDID output combining', and anything else that might be relevant.

## Extra programs under development for assisting CANDID operations
1. HMM database formatting code needs to be updated and included here.
2. CANDID output clustering/combining program. Currently, CANDID's weakness is that HDBSCAN can give variable results depending on specified settings. Combining multiple runs is a logical solution, but doing this isn't entirely easily. I think that developing something based on Hammock to combine the HMMs from multiple runs would be sensible. This would need to be able to combine clusters together without redundant/overlapping sequence, and identify novel clusters present in only one/some settings configurations. In short, difficult, but not as big of an undertaking as CANDID itself was.
