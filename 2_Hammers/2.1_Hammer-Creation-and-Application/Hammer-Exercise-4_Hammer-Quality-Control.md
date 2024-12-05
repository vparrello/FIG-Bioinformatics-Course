# Hammer Exercise 4 - Hammer Quality-Control

* Objective: Filter out candidate hammers that occur multiple times
within the complete sequences of the RepGen genomes.

So far we have only looked at the Kmers within the PheS SOUR sequences;
however, the PheS is only one gene out of several thousand genes
within a typical bacterial genome, and it is only about 
1000 basepairs long, whereas a typical bacterial genome is 
several million basepairs long. Hence, it is possible
that a Kmer that was found exactly once in exactly one
representative genome when only looking at the sequences
for a set of SOURs within a RepGen set
might stilll occur more than once when looking at the full-length 
contig sequences for that RepGen set. 
Hence, once a set of hammer-candidates has been found,
we must doublecheck them to filter out candidates
that also occur outside of the selected SOUR genes in any of the representative genomes.