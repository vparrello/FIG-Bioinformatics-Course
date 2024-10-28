# Hammers as Genetic Barcodes Exercise 2 - New Genome and its Representative

Objective: Given a set of representative genomes, find its associated "Genetic Barcodes" ("Hammers"), and use them to determine which representatives are most similar to the genomes found within a genetic sample.

Normally when working with bioinformatics data, we are given a sample containing one or more genomes that have been annotated. By "annotated", we mean that the sequence-data within the sample has been compared to the rest of the known universe of genomes, the data may have been been sorted into "bins" based on the most likely reference or representative genome or set of genomes that the user is interested in comparing them to, the sequence-data may have been "assembled" into longer sequences called "contigs", and the gene-containing regions on these contigs may have been located and had roles assigned to them. However, complete annotation of a sample can be an expensive and time-consuming process. Often, we are only interested in performing a "census" of which genera are present in a sample. We shall see that a set of "Hammers" will allow us to quickly perform such a "census", sidestepping the more complex and expensive complete annotation process.

In this exercise, we will work through a simplified "toy problem" that will illustrate the steps involved in constructing a set of hammers, and using these hammers to analysze a sample. You will be given a sample that contains sequences from a single "Mystery Genome" that have not been annotated, and a set of sequences implementing a single SOUR within a set of Representative Genomes. We will then walk you through a simplified version of the steps required to find the hammers associated with that single SOUR. Finally, we will search the "Mystery Genome" to see if it contains any "hammers". If the "Mystery Genome" contains many of the same "hammer" barcodes as one of the representative genomes, and very few "hammers" for any other representative genome, then this indicates that there is a high probability that the "Mystery Genome" and the representative genome that contributed the largest number of "hammer hits" are related.

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

1. First we need to make a set of hammers that can be used as "barcodes" for the selected SOUR. We will do this by extracting all the singly-occuring 20-mers from the DNA-sequences of a particular gene with the role-abbreviation of "PheS" (short for "Phenylalanine tRNA Synthase"). How do you think one might do that?
* Hint: Can you remember which datatype would be appropriate for associating a string with the number of times that it occurs?

    Load in the file `Definitions.html` as in previous exercises. Then, ask Grimoire to write a python program named `hammer_creator.py` that will:

        * Accept as '-K' a mandatory Kmer-length command-line argument;
        
        * Read a tab-separated-value file from `STDIN`, skip the header-line, parse the first two columns into (genome_id, sequence) pairs, and strip off any newlines that might have been attached to the sequence

        * Find all of the Kmers that occur exactly once in exactly one genome; these Kmers are the "Hammers"

        * Print a two column tab-separated table of the hammers and the `genome_id` that the hammer was found in, with column-headers of "hammer" and "genome_id".

    Paste Grimoire's program into the code-template `Code/hammer_creator.py`.
    The template also has a "block comment" section reserved for the program's pseudocode, so paste in Grimoire's pseudocde if it generated it, else ask it to generate pseudocode for the program if it didn't. Once you are done copying and pasting, save the code using the "Save" menu-item under the "File" menu.

2. Next, we need to extract the (genome_id, sequence) pairs from  `rep10.seqs.tbl`. You already have a program in your "toolkit" that performs this function, it is called `Code/cmd_tsv_select_columns.py`. The relevant columns in `rep10.seqs.tbl` are `genome_id` and `seed_dna`. In `FASTA-Exercise-2`, you learned how to "pipe" the output from one command to the input of another command, so you have all the tools you need:

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

## Caveats

In this exercise we have presented a greatly simplified "toy version" of the the hammer-creation process to illustrate the basic concepts behind hammers. The "production" version of the code is significantly more complex; we use the sequences for the 20 SOURs that we have found to yield the most reliably informative hammers, and we also cross-check the entire genome in addition to the SOUR sequences to make certain that the hammer does not occur somewhere else within the genome, since by definition a hammer can only occur once within one representative genome and must be contained within a SOUR. We have avoided these additional complexities within this exercise in order to keep this exercise simple, but in a later exercise, you will refine your prompt to take these additional complexities into account.


## Self Check

Below are the list of genomes in the representative genome set. Use the counts for the mystery genome to check your work.

The "Mystery Genome" should have `511145.12` "Escherichia coli str. K-12 substr. MG1655" as its highest-scoring representative genome.

### NOTE: Table of counts is still missing !!!