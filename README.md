# CANDID
Cluster-based Alignment-free Novel Domain Iterative Discovery

This is a program undergoing sporadic development. Its intent is to discover novel domains by excluding known domain regions from proteins and subsequently utilising sensitive MMseqs2 search to find potential domain regions which are clustered based on alignment-free distance measures and iteratively used as HMMER queries until convergence is reached.

Changes that have occurred recently:

1. MMseqs2 is integrated definitively as the searching program.
2. Old, poorly written code has been improved where possible.
3. Domain identification from MMseqs2 tabular output is modified to use an improved peak finding approach.
4. System now self-populates new domains each iteration on the basis of the HMMER domtblout file while still being able to keep track of changes. I am doing this with a simple coordinate comparison system.
5. Hammock is integrated as a possible option for clustering. From testing my original approach seems better and is certainly a lot faster.
6. MSA scoring function means we shouldn't get low-quality clusters in our output.

Directions for improvement:
1. More code tidying - requires me to actually finalise the program so I can strip out testing code.
2. Update benchparse system to work with the newer version of the code.
3. Make changes to user input to be less confusing.
