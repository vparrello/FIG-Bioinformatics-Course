# Hammers as Genetic Barcodes Exercise 2 - New Genome and its Representative

Objective: Given a set of represetative genomes, find its associated "Genetic Barcodes" ("Hammers"), and use them to classify the genomes within a sample.

Normally when working with bioinformatics data, we are given a sample containing one or more genomes that have been annotated. By "annotated", we mean that the sample has been compared to the rest of genome space, it has been sorted into "bins" based on the most likely reference or representative genome or genomes to be compared to, the sample may have been "assembled" into contigs, and the gene-containing regions on the contigs may have been located and had roles assigned to them. If a sample has not yet been annotated, then we can use our hammers to try and find which genomes it probably contains.

In this exercise you will be given a sample that has not been annotated. We will then walk you through the steps of finding the hammers associated with the genes for a single SOUR for a representative set of genomes, and finally we will search a "Mystery Genome" to see if it contains any hammers. If the "Mystery Genome" contains many of the same barcodes as a representative genome, and very few barcodes for any other representative genome, then this indicates that there is a high probability that the "Mystery Genome" and the representative genome that contributed the largest number of "hammer hits" are related.

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

1. First we need to make a set of hammers that can be used as "barcodes" for at least one SOUR. We will do this by extracting all the singly-occuring 20-mers from the DNA-sequences of a particular gene with the role-abbreviation of "PheS" (short for "Phenylalanine tRNA Synthase"). How do you think one might do that?
* Hint: Can you remember which datatype would be appropriate
for associating a string with the number of times that it occurs?

    Load in the file `Definitions.html` as in previous exercises. Then, ask Grimoire to write a python program named `hammer_creator.py` that will:

        * Accept as '-K' a mandatory Kmer-length command-line argument;
        
        * Read a tab-separated-value file from `STDIN`, skip the header-line, parse the first two columns into (genome_id, sequence) pairs, and strip off any newlines that might have been attached to the sequence

        * Find all of the Kmers that occur exactly once in exactly one genome; these Kmers are the "Hammers"

        * Print a two column tab-separated table of the hammers and the `genome_id` that the hammer was found in, with column-headers of "hammer" and "genome_id".

    Paste Grimoire's program into the code-template `Code/hammer_creator.py`.
    The template also has a "block comment" section reserved for the program's pseudocode, so paste in Grimoire's pseudocde if it generated it, else ask it to generate pseudocode for the program if it didn't. Once you are done copying and pasting, save the code using the "Save" menu-item under the "File" menu.

2. Next, we need to extract the (genome_id, sequences) pairs from  `rep10.seqs.tbl`. You already have a program in your "toolkit" that performs this function, it is called `Code/cmd_tsv_select_columns.py`. The relevant columns in `rep10.seqs.tbl` are `genome_id` and `seed_dna`. In `FASTA-Exercise-2`, you learned how to "pipe" the output from one command to the input of another command, so you have all the tools you need:

```
python Code/cmd_tsv_select_columns.py genome_id seed_dna < Data/rep10.seqs.tbl | Code/hammer_creator.py -K 20 > Data/rep10.hammers.tbl
``` 

3. We now need to compare the DNA sequences in `MysteryGenome1.fna` to your set of hammers.
Ask Grimoire to write a program that will:

        * Accept a mandatory TSV "hammers" filename as an argument.

        * Skip the header-line of the hammer-file and then read the first and second columns into a dictionary as a `hammer` and `genome_id`, respectively.

        * Determine the Kmer-length `K` of the hammers from the hammers dictionary.
        (NOTE: all the hammers in the file will have the same length.)

        * Use BioPython to read the sequences of the genome from `STDIN`.

        * For each sequence, extract all possible Kmers, and if a Kmer is a hammer, increment the score for its associated `genome_id`; then repeat this operation on the reverse-complement of that sequence, since a gene can face in either direction.

        * Print out a TSV file of the `genome_ids` found and their associated scores, sorted by decreasing score.
    
        * Finally, please translate the python code into pseudocode.

6. Once Grimoire is done, please paste its pseudocode and code into `Code/hammer_compare.py` and save it as usual.

7. Run your script on the `MysteryGenome1.fna`. The result should be a table of representative genome-IDs followed by how many hammers from that representative hit the "Mystery Genome". Take the genome-ID with the top number of hits as the most probable representative for the "Mystery Genome".


Congratulations! You have constructed the genetic barcodes for the PheS SOUR for the set of `rep10` representative genomes, and then used them to identify the `rep10` genome `511145.12` as the most likely representative for your "Mystery Genome"!




### Self Check

Below are the list of genomes in the representative genome set. Use the counts for the mystery genome to check your work.

The "Mystery Genome" should have `511145.12` "Escherichia coli str. K-12 substr. MG1655" as its highest-scoring representative genome.