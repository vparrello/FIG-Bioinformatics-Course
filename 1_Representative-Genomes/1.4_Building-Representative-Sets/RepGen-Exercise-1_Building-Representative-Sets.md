# RepGen Exercise 4 - Building a Representative Set

Objective: Combining what we have learned to create a program that builds a set of Representative Sequences from an input FASTA file.

## Materials

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

```
FIG-Bioinformatics-Course/
├── 1_Representative-Genomes/
│   └── 1.4_Building-Representative-Sets/
│       └── RepGen-Exercise-4_Building-Representative-Sets.md (you are here)
└── Data/
    └── Universe.fasta
```

## Overview:

As a result of the rapidly decreasing cost of sequencing genetic data, we are now at the stage where on the order of a million  bacterial genomes have have been sequenced and deposited in the public databases, and several billion associated gene-sequences, leading to an abundances of riches that have become difficult to manage.
There has thus developed a need to construct manageably small sets of "representatives" that capture the diversity of the full set of genomes and genes.

We define a "Set of Representatives" as follow:

* Let U be the Universe of entities (in our case, genomes or sequences) being considered.

* Let there be some measure of similarity between two entities.

* Then a "Set of Representatives" or "RepSet" is a subset of U such that:

    - no two members of the RepSet are more similar than some threshold similarity `Sim`, and

    - every member of U is similar to at least one member of the RepSet.

When the set of entities being represented is a set of genomes, we usually refer to it as a "RepGen Set" for short, and a member of the "RepGen Set" is called a "RepGen". In this exercise, we will look at the problem of constructing a set of "Representative Sequences".

Many measures of similarity can be used, but the measure we will use in this exercise is "Number of Kmers in Common". 

* NOTE: In previous exercises, we have worked with data-files having names like `rep10` or `rep200`; each of these files were obtained by constructing a set of representative genomes, and the number in the name was the threshold number of Kmers in common. For example, no two sequences in `rep10` have 10 or more short Kmers in common; no two members of `rep200` have 200 or more short Kmers in common, etc. In our "Production RepGen Sets" such as `rep20` or `rep200`, we have selected the protein-sequence that implements the role "Phenylalanine tRNA Synthetase" (abbreviated as "PheS") as the sequence to base our RepGen Sets on, and "number of protein 8-mers in common between two sequences" as our measure of similarity.
"PheS" is an example of a "Singly-Occurring Universal Role" (SOUR), which is a role implemented by exactly one gene (it is "Singly Occuring"), and that can be expected to be found in every genome (it is "Universal").

There are several algorithms for building RepSets. In this course, we will be using the "Stingy Addition" algorithm. The pseudocode for "Stingy Addition" is as follows:
```
  RepGenSet = ()   #... i.e. the "empty list"
  ToCheck   = U    #... The "Universe of Genomes"
  Foreach G1 in ToCheck:
    If no G in RepGenSet is more similar to G1 than MinSim:
      Add G1 to RepGenSet
            	
  Return RepGenSet
```

## Exercises:

1. Attach the `Definitions.html` file as in previous exercises.

2. The following prompt provides the program specification for the representative-set algorithm:

```
I will now give you a description of the command-line interface for the program `build_representative_set.py` to compute a "set of representative sequences" (RepGen set) using the "Stingy Addition" algorithm as defined in the uploaded definitions-file.

The script should accept the following mandatory command-line arguments, in both long-form and the specified short-form:
* Kmer-length  (integer, short argument-name '-k')
* Sim          (similarity threshold, short argument-name '-s')
* input FASTA-file   (filename, short argument-name '-f')
* output RepSeq-file  (filename, short argument-name '-r')

The measure of similarity to be used is "Number of Kmers in common".

The program should use BioPython to read the FASTA-file into a list, and then sort the sequence-objects by order of decreasing sequence-length, resolving length-ties by sequence, and then by sequence-ID.

The main body of the program should construct a subset of the input sequences that satisfies the provided definition of a "Representative Set" (RepGen set).

When you are done, please write out the representative set to the RepSeq-file in FASTA format.
```

Save the program that Grimoire generates as `build_representative_set.py`.

3. Run the program as follows:

```
python Code/build_representative_set.py -k 8 -s 10 -f Data/Universe.fasta -r Data/myrep10.faa
```

Check your result by running:
```
python Code/fasta_reader.py < Data/myrep10.faa
```
The result should end with:
```
Number of sequences: 155
Average sequence length: 363.73
```

## COMMENT: Why the Sort?

In the above exercise, we asked Grimoire to sort the sequences into
a particular order before processing them; you may be wondering,
why did we do that? The short answer is "reproducibility".

The representative-set algorithm processes the sequences in the order
that they were read in from the input-file, and at each step decides
whether to keep or reject a sequence based on how similar it is
to the sequences that were already accepted as representatives.
We would like the set of representatives selected to be "stable",
in the sense that we do not want the set of representatives selected
to depend on the order that they were seen within the input-file.
By sorting the sequences before processing them, we ensure that
the same set of representatives will be selected regardless of the order
that they appear within the input file.
