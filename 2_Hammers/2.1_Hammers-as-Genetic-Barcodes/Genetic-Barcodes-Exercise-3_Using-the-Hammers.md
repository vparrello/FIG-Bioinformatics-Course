# Hammers as Genetic Barcodes Exercise 3 - Using the Hammer Set

Objective: Practice using the hammers against two new genomes and determine if they are close

In the last lesson, we created a set of hammers(genetic barcodes) for a single genome and determined how close it was to your representative. You then found out which genome gets hit the most in terms of hammers and labelled the mystery genome as a neighbor of your representative. 

Follow the same pattern as the lesson before and find which representative genome is closest to your mystery genome. Then we give you a sample with multiple genomes in it. Compare your hammer results of what is in the sample to what the solution says to see how accurate your hammers are.

## Materials

MysteryGenome2.fasta
Rep10.list.tbl
Rep10.seqs.tbl 
rep10.hammers.tbl - output from the last lesson

## Exercises

1. Recall your `Code/cmd_tsv_select_columns.py`, `Code/hammer_compare.py`, and `Code/hammer_creator.py` tools. Review them for their use case and inputs. 

2. Since we have already created the hammer set, we do not necessarily need to use it again. Using `Code/hammer_compare.py` exclusively allows you to save time in the process. Use this program to get an output for `MysteryGenome2.fasta` to help review your process from the previous lesson.

3. Next we will use this same process to find the hammer counts for a mystery sample from BV-BRC. Review how to use their commmand line tool to download the sample "SRR19064441". Be sure to move it into the `Data` folder so that your programs can more easily access it.

4. Run your `Code/hammer_compare.py` program against the contents of the sample `SRR19064441.fasta`. This will have multiple genomes inside of it so do not eliminate any genomes from consideration when looking at your output.

4. Currently we have a few counts of how many hammers we are expecting to be hit from a sample. But without doing this hundreds of times, it would be very hard to have an accurate picture of how many hammer hits you need to verify its presence. Ask Grimoire to explain to you how a genetic barcode might be inaccurate or inconclusive.

5. Tuning our barcodes so that they only hit our target is crucial to having a useful tool. Without this, we might be casting a wider area than necessary, hitting genomes that are not close, or classifying not enough genomes to be accurate. Ask Grimoire to describe some ways that we can tune our hammers to accurately depict our genome space.

In the next unit, we will explore tuning the hammers to be more precise and accurate by determining how worthy a hammer is.

### Self Check

TODO put answers here
2. Mystery Genome: 
3. Table of genomes in the sample
