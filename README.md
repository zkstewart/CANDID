# CANDID
Cluster-based Alignment-free Novel Domain Iterative Discovery

This is a program undergoing sporadic development. Its intent is to discover novel domains by excluding known domain regions from proteins and subsequently utilising sensitive MMseqs2 search to find potential domain regions which are clustered based on alignment-free distance measures and iteratively used as HMMER queries until convergence is reached.

Changes that have occurred recently:

1. MMseqs2 is integrated definitively as the searching program. Ideally this program needs to be able to find alternative alignments during profile search, but for the time being its speed and sensitivity can't be matched by PSI-BLAST.
2. Old, poorly written code has been improved where possible. Functions are now more "drag-and-droppable" and the code should overall be more readable.
3. Domain identification from MMseqs2 tabular output is modified to use an improved peak finding approach. This should work better the NC Check system I was previously using.

Directions for improvement:
1. More tidying to code, especially in domclust. Determine whether I'm going to use the parse_joiner function (or something similar).
2. Trial different approaches for the hmm_grow function. I'd like to "unlock" the system to populate entirely new sequences on each iteration, but I still need a way of interrogating the data to find out when we've reached convergence. Unsure how to handle this. I might just opt to modernise hmm_grow to be more friendly and readable.
3. Test different approaches for domain clustering. Alignment-free clustering appears to give variable results which are difficult to control. Clustering needs to produce consistently accurate results. I will need to see what is currently available, alignment-free or otherwise. I wouldn't mind finding something that is more computationally intensive; if this proved to be a bottleneck, I could try other ways of growing the hmm model (such as using jackhmmer).
4. See how other programs compare. HC-HMM, Sub-HMMs, and Hammock might contain useful ideas; Hammock might actually be useful as a replacement for the clustering step in this pipeline? PSCAN also looks like it uses a very similar idea to my own (it uses DBSCAN*, I use HDBSCAN) but with a different scoring stage (unsure how it replaces my alignment-free similarity scoring).
