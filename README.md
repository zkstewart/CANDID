# CANDID
Cluster-based Alignment-free Novel Domain Iterative Discovery

This is a program undergoing sporadic development. Its intent is to discover novel domains by excluding known domain regions from proteins and subsequently utilising sensitive MMseqs2 search to find potential domain regions which are clustered based on alignment-free distance measures and iteratively used as HMMER queries until convergence is reached.

Directions for improvement:
1 - Replace PSI-BLAST with MMseqs2 by default more cleanly while keeping the ability to use PSI-BLAST as part of the program.
2 - Potentially scrap alignment-free clustering in favour of single-linkage clustering via MMseqs2. The reasons for this is that alignment-free clustering appears to give varying results depending on the number of results found. I want the program to be able to find short domains/motifs; relying on E-values for this is not ideal as they are sensitive to the alignment length whereas alignment-free distance measures should not be as sensitive to this. That said, MMseqs2 is already responsible for finding the domains/motifs, so there might not be any benefit to using a separate program for alignment-free clustering. Needs testing. If I do change this, the program name will need to change too. It will also make PSI-BLAST useless to keep as a function in this program if we use MMseqs for the clustering.
3 - Potentially remove HMMER3 iterative queries. From initial tests, MMseqs2 appears to be sensitive enough that this step no longer offers benefit whereas it used to discover new alignments when this program used PSI-BLAST. Also another thing that will require the program's name to change if removed.
