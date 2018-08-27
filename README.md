# CANDID
Cluster-based Alignment-free Novel Domain Iterative Discovery

This is a program undergoing sporadic development. Its intent is to discover novel domains by excluding known domain regions from proteins and subsequently utilising sensitive MMseqs2 search to find potential domain regions which are clustered based on alignment-free distance measures and iteratively used as HMMER queries until convergence is reached.

Changes that have occurred recently:

1. System now self-populates new domains each iteration on the basis of the HMMER domtblout file while still being able to keep track of changes. I am doing this with a simple coordinate comparison system.
2. Hammock is integrated as a possible option for clustering. From testing my original approach seems better and is certainly a lot faster.
3. MSA scoring function means we shouldn't get low-quality clusters in our output.
4. User input should be less confusing. 'help' and 'help-long' help to reduce the overwhelming amount of tunable parameters to just the set that are relevant to most users who will be using _what is hopefully_ my optimised parameters.

Directions for improvement:

1. Update benchparse system to work with the newer version of the code.
2. Figure out optimal default parameters - more testing!
3. Strip out testing code once I have figured out the optimal default parameters.
