# Hammers

A "hammer" is a 20-character DNA sequence that occurs precisely once in precisely one representative genome within a given RepGen set. It must also occur within a gene implementing a SOUR. Thus, a hammer acts as a genomic signature. One may think of a hammer as a type of "Genetic Barcode" that we can use to quickly classify an unknown genome, or to try to create a targeted antibiotic that will kill a pathogen while leaving other microbes alone. In this module, the goal is to create a set of hammers that can be used to quickly identify the closest representative to a genome, and to determine a set of representatives that characterize a metagenomic sample.

## 2.1 Hammers Creation and Application

### Hammer Exercise 1 - What are "Hammers", why do we need them, and how do we build them?
Discover what "hammers" are, how to construct them, and how they can be used to classify an unknown genome

### Hammer Exercise 2 - Using Hammers to find the Nearest Representative
Use your hammer tool to attempt to identify a "mystery genome".

### Hammer Exercise 3 - Using Hammers on a Metagenomic Sample
Use your hammer tool to identify which genomes are present within a "metagenomic sample", i.e., a sample that contains several genomes.

### Hammer Exercise 4 - Hammer Quality-Control

So far we have only looked at the Kmers in the PheS SOUR sequences.
Now we must make sure that a hammer doesn't occur outside of a SOUR in any of the representative genomes.

<<<<<<< HEAD
# Hammer Exercise 5 - Building Hammers Using More Than One SOUR
=======
## 2.2 Hammers Using More Than One Role
>>>>>>> 05ec586 (Intermediate save-point. I've split off "Multirole Hammers" into a new unit, "2.2_Hammers-Using-More-Than-One-Role", with corresponding revisions to 'SyllabusUnit2.md'.)

So far we have built hammers for a RepGenSet
using only its PheS SOUR sequences.
We can obtain more reliable genome assignments by requiring a consensus on which genomes are present
within a sample, by allowing a "jury" of several different SOURs
to "vote" on which RepGen genomes are closest
to a sample.

<<<<<<< HEAD
=======
### Multirole Hammer Sets Exercise 1 - Building Hammers Using More Than One SOUR

Enhancing the Hammer-Creation algorithm to support hammers
that are signatures of particular roles within a RepGen.

### Multirole Hammer Sets Exercise 2 - Applying Multirole Hammers

Enhancing the Hammer-Application algorithm to implement
"Jury Voting" for whether a RepGen is present within a sample.
>>>>>>> 05ec586 (Intermediate save-point. I've split off "Multirole Hammers" into a new unit, "2.2_Hammers-Using-More-Than-One-Role", with corresponding revisions to 'SyllabusUnit2.md'.)
