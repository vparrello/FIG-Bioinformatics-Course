# Hammers

### 2.1 Hammers as "Genetic Barcodes"
    * Why do we care about hammers?
        * Problem: We have a new genome and we need to find his closest representative
        * Find PheS and compare it to all the other PheS - Problem comes in if our genome is not annotated (takes a lot of time)
        * Take a single repgen, is our genome close to this one?
        *We need to look for a "hammer" which is a genetic bar code. This occurs in the PheS but nowhere else in the genome. If we find this hammer, then we know we have a PheS like our original one.
        * We find these hammers and see how many are in the new one - We don't really know what this number means
        *Grab a couple other genomes and see how many hammers hit. We should grab some that we know are close and some that we know are far. 

### 2.2 Hammer Tuning
    * Precision - Make sure the hammer doesn't occur in other representative genomes
        Do some fancy dancy stuff so it doesn't take too long (subset of repgens)
    * Accuracy - Picking the correct roles
        PHEs is our go to. What do the others do?
        Run the hammer program to look for hammers inside the roles
        We know they are bad if they find everyone or noone - check against jaccard similarities of the repgens
        We made the list for you because we're awesome
    * Worthiness - How many peers are hit by the hammer
        Different use cases for different hammer counts. 
        Hitting a wide area is good for classification - lots of hits means lots of peers are hit which is good for classifying a genome
        hitting a small area is good for antibiotics - few hits means less collateral damage from the antibiotic in your system


### 2.3 Creation of the Hammer Set
    * take everything you've learned, and create the hammer set for all the repgens. 
    * Use case for classification
    * Use case for antibiotics
    * Use case for closest genome
