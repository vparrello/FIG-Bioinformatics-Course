# Hammers

Hammers are "Genetic Barcodes" that we can use to try and create targeted antibiotics or as a faster way to classify and annotate genomes. In this module, the goal is to create a set of hammers that can be used to quickly identify genomes and use that data to project possible diseases in a human metagenomic sample.

## 2.1 Hammers as "Genetic Barcodes"

### Genetic Barcodes-1 Why do we need hammers and what are they?
Discover what "hammers" are, and how they can be used to classify a unknown genome

### Genetic Barcodes-2 

* Why do we care about hammers?

    * Problem: We have a new genome, and we need to find the closest representative genome (RepGen) within a given RepGen Set.

    * We could find the new genome's PheS and compare it to all the other PheSs - Problem comes in if our new genome is not yet annotated (takes a lot of time)

    * Suppose I have a given RepGen; is our new genome close to this Repgen?

    * We can avoid first locating the PheS by instead looking for a "hammer",
    which is a type of "genetic barcode".
    A "Hammer" is a 20-mer that is found:

        * within a PheS gene (or in general within a gene implementing a member from a set of SOURs)

        * occurs exactly once in exactly one representative genome,

        * is found nowhere else in any other representative genome.

    If we find a hammer in a new genome, then it provides strong evidence that we have found a PheS gene that is similar to the PheS gene in the representative genome that the hammer came from, because a 20-mer is a highly specific DNA pattern.

    * We find these hammers and see how many are in the new one - We don't really know what this number means

    * Grab a couple other genomes and see how many hammers hit. We should grab some that we know are close and some that we know are far. 

### 2.2 Hammer Tuning
    * Precision - Make sure the hammer doesn't occur in other representative genomes

        * Do some fancy dancy stuff so it doesn't take too long (subset of repgens)

    * Accuracy - Picking the correct roles

        * PHEs is our go to. What do the others do?

        * Run the hammer program to look for hammers inside the roles

        * We know they are bad if they find everyone or noone - check against jaccard similarities of the repgens

        * We made the list for you because we're awesome

    * Worthiness - How fraction of of a regen's peers are hit by the hammer

        * Different use cases for different hammer counts. 

        * Hitting a wide area is good for classification - lots of hits means lots of peers are hit which is good for classifying a genome

        * hitting a small area is good for antibiotics - few hits means less collateral damage from the antibiotic in your system


### 2.3 Creation of the Hammer Set
    * take everything you've learned, and create the hammer set for all the repgens. 

    * Use case for classification

    * Use case for antibiotics

    * Use case for closest genome
