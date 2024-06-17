# Hammers as Genetic Barcodes Exercise 2 - New Genome and its Representative

Objective: Find a "Genetic Barcode" for a genome in order to classify one inside a sample.

Normally when working with bioinformatics data, we are given a sample containing one or more genomes that has been annotated. This means that the sample has been compared to the rest of genome space, it has been sorted into "bins" based on the most likely reference or representative genome or genomes to be compared to, the sample may have been "assembled" into contigs, and the probably gene-regions on the contigs may have be located and had roles assigned to them. If a sample has not been annotated, then we need to use our hammers to try and find which genomes it probably contains.

In this exercise you will be given a sample that has not been annotated. We will then walk you through the steps of finding the hammers from a single representative genome to see if any of those hammers are present in the mystery genome. If they are, that tells us that they also have the same barcodes in their genome and therefore are more likely to be the same genome.

## Materials

TODO insert a sample data here with planted genomes in it. -  MysteryGenome.fasta
TODO have a sample representative genome sample for creating hammers. - Rep10.list.tbl, Rep10.seqs.tbl, genomeid.fasta
TODO program that will grab 20-mers without killing the memory in most machines

## Exercises

1. First we need to make a set of hammers that can be used as "barcodes" for at least 1 representative genome. We will do this by extracting all the singly-occuring 20-mers from the DNA-sequences of a particular gene with the role-abbreviation of "PheS" (short for "Phenylalanine tRNA Synthase").
How do you think one might we do that?
* Hint: Can you remember which python datatype can only contain a single copy of a string?

Ask Grimoire to write a python program that will:

* Accept a Kmer-length `K` as a mandatory command-line argument;
* Read a tab-separated list of (genome_id, sequence) pairs from `STDIN`, skipping the header-line;
* Find all of the Kmers that occur exactly once in exactly one sequence (these Kmers are the "Hammers"), and the `genome_id` that the Kmer occurred in;
* Print a two column tab-separated table of the hammers and their  `genome_id`, with column-headers of "hammer" and "genome_id".

Paste Grimire's program into the code-template `Code/hammer_creator.py`.
The template also has a "block comment" section reserved for the program's pseudocode, so paste in Grimoire's pseudocde if it generated it, and ask it to generate pseudocode for the program if it didn't. Once you are done copying and pasting, save the code using the "Save" menu-item under the "File" menu.

2. Next, we need to extract the (genome_id, sequences) pairs from  `rep10.seqs.tbl`. You already have a program in your "toolkit" that performs this function, it is called `cmd_tsv_select_columns.py`. The relevant columns in `rep10.seqs.tbl` are `genome_id` and `seed_dna`.
In `FASTA-Exercise-2`, you learned how to "pipe" the output
from one command to the input of another command, so you have all the tools you need:

```
python Code/cmd_tsv_select_columns.py genome_id seed_dna | Code/hammer_creator.py -K 20 > Data/rep10.h20.tbl
``` 

3. We now need to compare the MysteryGenome to your set of hammers. Ask Grimoire to make another program. One that will read a sequence, make a 20mer, compare it with the set, and then add to a counter every time it has a match. Paste this code into hammer_compare.py.

6. Add hammer_compare.py to your framework from Step 3. If you need help, ask Grimoire to adjust the format to include all the programs. 

7. Run your script on the MysteryGenome.fasta. You should have a count print out of how many hammers were hit. 


Congratulations! You have created the genetic barcodes from the PheS SOUR for the TODO REPGEN NUMBER HERE genome. We need to do the same for the rest of the representative genome set and see how the counts differ. Please complete the hammers for all of the representative genomes so that you are prepared for the next exercise.


### Self Check

Below are the list of genomes in the representative genome set. Use the counts for the mystery genome to check your work.

TODO put answers here