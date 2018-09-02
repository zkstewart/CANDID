# CANDID
Cluster-based Alignment-free Novel Domain Iterative Discovery

This is a program still undergoing development - it's getting pretty close to being "finished", though! Its intent is to discover novel domains by excluding known domain regions from proteins and subsequently utilising sensitive MMseqs2 search to find potential domain regions which are clustered based on alignment-free distance measures and iteratively used as HMMER queries until convergence is reached.

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
* Hammock can optionally replace alfpy-HDBSCAN clustering (https://github.com/krejciadam/hammock)

Note that, while Hammock is offered as an option for clustering, from initial testing alignment-free HDBSCAN (ALF-HDBSCAN) clustering seems like a better option in most scenarios. It is a lot faster which is important for the iterative nature of CANDID to converge upon a solution quicker. Moreover, it is more likely to detect clusters in a data set where Hammock may struggle. If you have the time and computational resources, there's no reason why you shouldn't be able to try out both solutions and see what works best.

**IMPORTANT**: I haven't figured out how to get Hammock to run on Windows using Cygwin. If you can do this, let me know how. Otherwise, you may be limited to ALF-HDBSCAN on this operating system.

## Progress

Changes that have occurred recently:

1. MSA curation functions are being included to remove bad quality clusters and remove outliers from clusters. Unfortunately, I do not know of a way to get HDBSCAN to consistently output "perfect" clusters without having the odd outliers here and there. I don't think this is possible, and it is just part of the beast we're working with - really short protein sequence clustering is hard. Hammock struggles in the exact same ways which is reassuring and disappointing that no perfect solution exists.

Nonetheless, I have included msa_trim to remove extensions surrounding the putative domain region and remove excessively gappy sequences. msa_score acts as a simple heuristic measure to remove obviously poor quality clusters. odseq_outlier_detect is a WIP function which will call ODseq on each cluster and remove outliers. I would have liked to implement something similar to this directly within Python for convenience/speed of use/to not require R as another external program to the already long list, but I can't imitate it well enough. Combined, these functions allow us to use relaxed HDBSCAN clustering (maximising our ability to detect clusters) while culling the poor quality clusters that inevitably result from this.
2. coord_lists_overlap_cull has been added to ensure that we don't rediscover non-novel domain regions. This resulted in an immediate and drastic improvement to the results.

Directions for improvement:

1. Finish msa curation system
2. Update benchparse system to work with the newer version of the code.
3. Figure out optimal default parameters - more testing!
4. Strip out testing code once I have figured out the optimal default parameters.
