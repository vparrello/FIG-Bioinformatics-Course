# Hammer Exercise 2 - Using Hammers to find the Nearest Representative

Objective: Given a set of representative genomes,  a set of associated "Hammers", and an unknown genome, determine which representative the unkonwn genome is most similar to.

In this exercise we will generate code that will search the "Mystery Genome" to see if it contains any "hammers". If the "Mystery Genome" contains many of the same "hammer" barcodes as one of the representative genomes, and very few "hammers" for any other representative genome, then this indicates that there is a high probability that the "Mystery Genome" and the representative genome that contributed the largest number of "hammer hits" are related.

## Materials

```
FIG-Bioinformatics-Course/
├── 2_Hammers 
│   └── 2.1_Hammers-as-Genetic-Barcodes/
│       └── Genetic-Barcodes-Exercise-2-New-Genome.md(you are here)
└── Data/
    ├── MysteryGenome1.fna
    ├── rep10.list.tbl
    └── rep10.seqs.tbl
```

## Exercises

1. We need to compare the DNA sequences in `MysteryGenome1.fna` to the set of hammers
you generated in the previous exercise. Ask Grimoire to write a program that will:

    * Accept a mandatory TSV "hammers" filename as an argument.

    * Skip the header-line of the hammer-file and then read the first and second columns into a dictionary as a `hammer` and `genome_id`, respectively.

    * Determine the Kmer-length `K` of the hammers from the hammers dictionary.
        (NOTE: all the hammers in the file will have the same length.)

    * Use BioPython to read the sequences of the genome from `STDIN`.

    * For each sequence, extract all possible Kmers, and if a Kmer is a hammer, increment the score for its associated `genome_id`; then repeat this operation on the reverse-complement of that sequence, since a gene can face in either direction.

    * Print out a TSV file of the `genome_ids` found and their associated scores, sorted by decreasing score.
    
    * Finally, please translate the python code into pseudocode.

2. Once Grimoire is done, please paste its pseudocode and code into `Code/hammer_compare.py` and save it as usual.

3. Run your script on the `MysteryGenome1.fna`. The result should be a table of representative genome-IDs followed by how many hammers from that representative hit the "Mystery Genome". Take the genome-ID with the top number of hits as the most probable representative for the "Mystery Genome".


Congratulations! You have constructed the genetic barcodes for the PheS SOUR for the set of `rep10` representative genomes, and then used them to identify the `rep10` genome `511145.12` as the most likely representative for your "Mystery Genome"!

## Caveats

In this exercise we have presented a greatly simplified "toy version" of the the hammer-creation process to illustrate the basic concepts behind hammers. The "production" version of the code is significantly more complex; we use the sequences for the 20 SOURs that we have found to yield the most reliably informative hammers, and we also cross-check the entire genome in addition to the SOUR sequences to make certain that the hammer does not occur somewhere else within the genome, since by definition a hammer can only occur once within one representative genome and must be contained within a SOUR. We have avoided these additional complexities within this exercise in order to keep this exercise simple, but in a later exercise, you will refine your prompt to take these additional complexities into account.


## Self Check

Below are the list of genomes in the representative genome set. Use the counts for the mystery genome to check your work.

The "Mystery Genome" should have `511145.12` "Escherichia coli str. K-12 substr. MG1655" as its highest-scoring representative genome.

### NOTE: Table of counts is still missing !!!