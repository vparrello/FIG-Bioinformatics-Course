# Hammers

A "hammer" is a 20-character DNA sequence that occurs precisely once in precisely one representative genome within a given RepGen set. It must also occur within a gene implementing a SOUR. Thus, a hammer acts as a genomic signature. One may think of a hammer as a type of "Genetic Barcode" that we can use to quickly classify an unknown genome, or to try to create a targeted antibiotic that will kill a pathogen while leaving other microbes alone. In this module, the goal is to create a set of hammers that can be used to quickly identify single genomes, and to determine which genomes are present within a metagenomic sample.

## 2.1 Hammers Creation and Application

### Hammer Exercise 1 - What are "Hammers", why do we need them, and how do we build them?
Discover what "hammers" are, how to construct them, and how they can be used to classify an unknown genome

### Hammer Exercise 2 - Using Hammers to find the Nearest Representative
Use your hammer tool to attempt to identify a "mystery genome".

### Hammer Exercise 3 - Using Hammers on a Metagenomic Sample
Use your hammer tool to identify which genomes are present within a "metagenomic sample", i.e., a sample that contains several genomes.

## 2.2 Hammer Tuning

### Hammer Tuning-1 Precision
Make sure the hammer doesn't occur in other representative genomes

### Hammer Tuning-2 Accuracy
Picking the correct roles mean that the hammers are more accurate for your purpose.  

### Hammer Tuning-3 Worthiness
How many "peer genomes" of a representative are hit by the hammer you chose?

* Hitting a wide area is good for classification - lots of hits means lots of peers are hit which is good for classifying a genome

* Hitting a small area is good for antibiotics - few hits means less collateral damage from the antibiotic in your system

## 2.3 Creation of the Hammer Set

### Creation-1 Create a hammer set
Create a set for all available repgen genomes to use for both classification, and antibiotics

### Creation-2 Use Case: Classification
Use the hammer set to classify genomes inside of a sample.

### Creation-3 Use Case: Antibiotics
Use the hammer set to target a group for antibiotics

### Creation-4 Use Case: Closest Genome
Use the hammer set to find the genome closest to the mystery genome
