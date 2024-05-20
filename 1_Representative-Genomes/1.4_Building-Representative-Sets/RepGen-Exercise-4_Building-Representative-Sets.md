# RepGen Exercise 4 - Building a Representative Set

Objective: Combining what we have learned to create a program that builds a set of Representative Sequences from an input FAST file.

## Materials

(Fill this in later.)

1. Attach the `Definitions.html` file as in previous exercises.

2. Reformulate the following prompt as a set of program specifications

```
I am now about to give you a description of the command-line interface for a program to compute a "set of representative sequences" (RepGen set) using "Stingy Addition" as defined in the uploaded definitions-file.

The script should accept the following mandatory command-line arguments, in both long-form and the specified short-form:
* Kmer-length  (integer, short argument-name '-k')
* MinSim       (percent, short argument-name '-m')
* input FASTA-file   (filename, short argument-name '-f')
* output RepSeq-file  (filename, short argument-name '-r')

MinSim should be converted from a percent to a fraction.

The program should use BioPython to read the FASTA-file into a list, and then sort the sequence-objects by order of decreasing sequence-length, resolving length-ties by sequence, and then by sequence-ID.

The main body of the program should construct a subset of the input sequences that satisfies the provided definition of a "Representative Set" (RepGen set), using the Kmer jaccard-similarity as your measure of similarity.

When you are done, please write out the representative set to the RepSeq-file in FASTA format.
```