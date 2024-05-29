# Kmer Exercise 3 - Find Nearest Reference-Sequence

Objective: Determine which representative sequence is closest to a query sequence.

A common task in bioinformatics is to determine which sequence in a reference database is "most similar" to a query-sequence according to one of several measures of "similarity". For example, you might have a newly-sequenced genome of unknown species, and you'd like to know what its "closest relative" is. In this exercise, we will implement a simplified version of this task: We will read in a set of reference-sequences from a file (or as we call them, "representative sequences" for reasons that will become obvious in Section 1.4), and then for each entry in a set of test-sequences, report the nearest reference-sequence based on jaccard-similaritity.

## Materials: 

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

```
FIG-Bioinformatics-Course
├── Definitions.html
├── Data/
│   └── rep10.seqs.faa
└── 1_Representative-Genomes
    └── 1.3_Kmers-and-Jaccard-Similarities
        ├── Kmer-Exercise-3_Finding-Nearest-Neighbors.md (you are here)
        ├── Data/
        │   └── test-sequences.faa
        └── Solutions
            └── nearest-reference-sequence_solution.py
```

## Exercises

1. Attach the `Definitions.html` file as in provious exercises, again using the prompt:
    ```
    I am going to give you some definitions in an attached file; you don't need to respond to them, just learn them. I'm then going to ask you some questions.
    ```

2. Ask Grimoire to write a program that implements the following specifications:
    * Accepts a mandatory integer argument, short-name `-K`, long-name `--kmer-length`
    * Accepts a mandatory FASTA filename argument, short-name `-R`, long-name `--RepSet`
    * Reads in the FASTA RepSet file using BioPython, and stores the sequence-ID, sequence-description, and the sequence
    * Reads a set of FASTA sequences from STDIN, and for each input sequence reports as a TSV the input ID, and the most similar RepSet sequence-ID, sequence-description, and Kmer jaccard-similarity


## NOTES:

1.) Need to create the file `rep10.seqs.faa` and the test-sequences.

2.) What is the best way to facillitate the student checking their results? Do we include the sequence-description in the test-file and have Grimoire output the description, so that the student can compare the input description and the RepGenSet description? (In a real-world problem, the sequences will not come "pre-labeled".)