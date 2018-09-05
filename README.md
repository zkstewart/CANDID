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
  * If NOT running Hammock for clustering, alfpy is required https://github.com/aziele/alfpy (also known as Alfree)
  * If NOT running Hammock for clustering, HDBSCAN is required https://github.com/scikit-learn-contrib/hdbscan
* Hammock can optionally replace ALF-HDBSCAN clustering (https://github.com/krejciadam/hammock)

Note that, while Hammock is offered as an option for clustering, from initial testing alignment-free HDBSCAN (ALF-HDBSCAN) clustering seems like a better option in most scenarios. It is a lot faster which is important for the iterative nature of CANDID to converge upon a solution quicker. Moreover, it is more likely to detect clusters in a data set where Hammock may struggle. If you have the time and computational resources, there's no reason why you shouldn't be able to try out both solutions and see what works best.

**IMPORTANT**: I haven't figured out how to get Hammock to run on Windows using Cygwin. If you can do this, let me know how. Otherwise, you may be limited to ALF-HDBSCAN on this operating system.

## How to use

Calling CANDID.py on command-line like ```CANDID.py -h``` or ```CANDID.py -help-long``` will print details relevant to program operation. Here, I'll provide a more thorough description of the necessary file inputs, what parameters can be altered, and their impacts.

### Preparing the HMM database

Prior to runnning CANDID it is necessary to obtain a HMM database to exclude known domains from discovery. It's important that this database be as comprehensive as possible so as to reduce the chance we "discover" a domain that has already been described and modelled in a database, such as those curated by Pfam, SMART, NCBI's CDD, CATH, or SUPERFAMILY. The program generate_hmm_db.py is designed to assist in this task, converting the CDD and optionally (as a strong recommendation) the databases of SUPERFAMILY and CATH. These resources can be obtained from the following locations.

* NCBI's CDD includes models from Pfam, SMART, and a variety of other sources. Its download is automatically coordinated by the generate_hmm_db.py program, but it can be downloaded beforehand from ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/fasta.tar.gz.
* SUPERFAMILY requires you to sign-up and obtain a licence for its use. You can do this at http://supfam.org/SUPERFAMILY/downloads.html then follow the instructions to download the HMM library.
* CATH provides a file with naming format "cath-S35-hmm3-v#\_#\_#.lib.gz" within its releases ftp site (ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/all-releases). As of writing this (05-09-2018) the latest version is 4_1_0, and its full link is ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/all-releases/v4_1_0/sequence-data/cath-S35-hmm3-v4_1_0.lib.gz.

Once complete, you will need to download and install the above-mentioned prerequisites and you will be ready to begin using CANDID.

### Configuring CANDID

CANDID is designed to be fully user-configurable, which means there are a lot of parameters that can be specified. By default most such parameters are hidden as their defaults should be used in most circumstances, but you can reveal these with the ```-help-long``` program call.

All arguments can be specified on command-line. Since this is cumbersome, CANDID will automatically format a .config file which can be used for subsequent runs. On the command-line, provide the program call ```-generate_config``` in addition to those arguments which are mandatory or configurable and this file will be generated if the chosen parameters pass a variety of checks to ensure their correctness. After this, the config file within the specified output directory will be automatically parsed upon resuming the program without ```-generate_config```, or you can keep a basic config file which you can point CANDID to with the ```-config``` argument. When CANDID parses a config file, any arguments provided on the command-line will overrule those in the .config file; this can let you make minor adjustments when pointing CANDID to a basic config file.

CANDID needs to know the location of external programs that it relies upon. Many of these do not need to be explicitly specified if they are discoverable within your PATH, but certain ones (which are pointed out in the program help text) must be explicitly defined.

Hammock is an optional choice for performing the protein sequence clustering step. If you choose to use Hammock, you do not need the alfpy or HDBSCAN Python packages, but you must provide the location of the Hammock.jar file and have a working installation of Java.

### CANDID pipeline and parameter choices

The first step of CANDID attempts to remove highly similar sequences from consideration with strict CD-HIT settings. This is necessary since high identity sequences will align end-to-end and may include alignments of multiple domains; CANDID instead wants to find partial alignments of parts of sequences which theoretically correspond to single domain regions shared between sequences. It is not recommended to change these parameters since they are designed to strike a balance between strictly removing highly similar sequences without reducing our chances of finding real domains excessively.

The second step of CANDID will run HMMER search to find known domain models (as per the HMM database possibly created by generate_hmm_db.py) within the CD-HIT clustered sequences. An E-value cutoff of 0.1 is recommended by default, and as such, this parameter is hidden by default. Any predictions made will result in the CD-HIT clustered sequences having these portions of their sequence masked.

The third step of CANDID involves the prediction of signal peptides using SignalP. These features may share high similarity with other sequences but are not considered domains and, as such, any predictions will be masked from the protein sequences. The only relevant parameter here is the organism type which is not hidden from a normal ```-h``` call and should be specified to be whatever is most appropriate.

The fourth step of CANDID involves seg and PSCOILS prediction of low-complexity regions and coiled coils, respectively. There are no relevant parameters for tuning this operation. Once again, predictions will be masked from the sequence since these regions, despite often sharing sequence similarity with other proteins, are not considered domains.

The fifth, sixth, and seventh steps of CANDID involve all-against-all sequence comparison using MMseqs2. The masked sequences will be converted into a MMseqs2 database and queried against each other to find local regions of similarity. Putative domain regions will be identified on the basis of these alignments by looking for 'peaks' in the amount of sequences that align against residues within a sequence. Theoretically, we would expect homologous sequences with weak similarity that made it through CD-HIT clustering to align end-to-end resulting in no noticeable peaks or troughs when plotting a line where X = sequence position and Y = alignment coverage. Residues which correspond to domain regions should have a higher amount of aligned sequences than inter-domain regions, resulting in peaks and troughs in our imaginary graph. The E-value used for returning significant alignments (0.1) is hidden from normal use since it is important that we be relaxed enough to find weakly conserved domains, but not so relaxed that we alter our signal:noise ratio and make later clustering more difficult and less powerful.

The eigth and final processing step of CANDID involves the iterative clustering loop which gives the program its name. Although Hammock can be substituted instead, by default an alignment-free clustering process (ALF-HDBSCAN) is used to compare putative domain regions obtained through parsing MMseqs2 results and to cluster these on the basis of similarity rather than some measure of identity. This provides us the ability to find domain regions that are only weakly defined, but is also the cause of some difficulty. All of the parameters behind this process are hidden since their defaults should work in most cases, but this is where user tuning should be employed to 1) find a set of parameters that works best for your data, and/or 2) to obtain the results from multiple parameter configurations and combine these using CANDID_combine [Forthcoming].

Without going into excessive detail, ALF-HDBSCAN clustering can be altered by changing these factors:

* Alignment-free scoring parameters:
  * "Word size" of similarity scoring can be configured. The recommended value range is 1-6 according to the alfree publication (doi:10.1186/s13059-017-1319-7) and is specified as 2 by default here.
  * Protein alphabet can be reduced to 15 or 11 characters. This can be useful in certain situations to identify domains with low identity but high conservation of residue-type (e.g. polarity characteristic); this is turned off by default.
  * The algorithm can be chosen out of normalised Google distance (default, alfree benchmarks suggests it is the best) as well as Bray-Curtis dissimilarity or Canberra distance. While alfree does provide (many) more distance measures, these three are the "best" according to alfree benchmark as well as my own internal testing.
* HDBSCAN clustering parameters:
  * Minimum cluster size can be specified according to the size of your dataset / the required amount of occurrences of a putative "domain" for it to be considered a true domain. A size of 3 is a recommended minimum since, while two occurrences of a conserved sequence may indicate that it is a domain, it likely isn't a domain in the sense that it is found in multiple and possibly unrelated gene families resulting from exon shuffling or other evolutionary processes.
  * Minimum sample size is somewhat confusing and its effects on clustering performance are far-reaching. By default I recommend you leave this at 2 since increasing or decreasing this seems to result in worse clusters for poorly understood reasons, but your own dataset might behave differently.
  * The algorithm used to cluster points can be changed to 'leaf' as opposed to 'excess of mass' which is the recommended default by HDBSCAN's creators. In some scenarios this results in improvement, but in others it might lead to what would be a single good domain being split into multiple smaller but highly homogenous clusters.
  * HDBSCAN by default does not allow a single cluster to be obtained for various reasons. You can change this on the command-line, but this should only be done if you ran HDBSCAN previously and didn't find any "good" domains or if you realistically expect only one novel domain to be present in your dataset (which is an unjustifiable assertion for the most part - we don't know how many novel domains exist unless the data is simulated).
* Iteration control parameter:
  * CANDID will, by default, continue to iteratively discover domains and modify their boundaries until it reaches a point where convergence occurs - no further modifications occur when we perform the same process more than once. This can take a long time in practice, possibly dozens of iterations might occur. For extremely large datasets this can be problematic, especially since in later iterations we typically only see small changes to domain boundaries and the removal of outliers rather than the discovery of new domain regions. You can specify a limit on the amount of iterations that can occur if CANDID would otherwise take "too long".
  * Iteration can be halted without killing the program mid-way through an iteration by creating a file called "CANDID_exit_marker" (no file extension) within the specified output directory. At the end of an iteration, if this file exists, CANDID will act as if we reached a limit on the number of iterations and will format the file outputs appropriately. Thus, you are able to run the program without a limit on convergence and decide when you think it has run long enough.

Within this iterative loop is an additional outlier detection method which utilises ODseq as well as reclustering with HDBSCAN using a sum-of-pairs scoring scheme. There are no relevant parameters for altering this behaviour. Clusters which do not score well and individual outlier sequences are removed during this process. These steps are important since HDBSCAN does tend to discover some poor quality clusters in addition to including outliers in clusters on occasion, especially when the input dataset is somewhat noisy. 

Finally, this program can be resumed by re-running the program and it will attempt to pick up where it left off. It uses a system of files with suffix '.complete' to note when a specific action has completed successfully. As a result, you can 'hack' this program to use your own inputs at any stage along the process so long as it is named appropriately. This isn't entirely recommended, but it's something that can be done if you understand what CANDID is doing and have some reason to believe you can do a specific step better.

## Development progress

Changes that have occurred recently:

1. MSA curation functions are completely implemented. ODseq's outlier predictions require additional validation with my own (rough) outlier detection method. This combination helps to temper these results to be less strict in certain situations.
2. Optimal default parameters are likely to be minsize = (2 or 3) and minsample = 2. minsize depends on the amount of data that is input; with small datasets use 2, and with larger ones use 3. minsample is tricky; in most cases 2 works best but in at least one case with a smaller dataset only minsample == 1 was able to recover domains.
3. Testing code is being stripped out since we're largely done with developing individual functions.
4. generate_hmm_db.py code update is complete. This acts as a helper program for making a database for excluding known domains from sequences necessary for CANDID to function properly.

Directions for improvement:

1. Do a bit more testing to see if I should change the convergence detection to stop in the first instance that no changes are detected.
2. Change how output files are presented/hide the processing files elsewhere.
3. Continued updating to this README. Include subheading 'CANDID output combining' and anything else that might be relevant.

## Extra programs under development for assisting CANDID operations
1. CANDID output clustering/combining program. Currently, CANDID's weakness is that HDBSCAN can give variable results depending on specified settings. Combining multiple runs is a logical solution, but doing this isn't entirely easily. I think that developing something based on Hammock to combine the HMMs from multiple runs would be sensible. This would need to be able to combine clusters together without redundant/overlapping sequence, and identify novel clusters present in only one/some settings configurations. In short, difficult, but not as big of an undertaking as CANDID itself was.
