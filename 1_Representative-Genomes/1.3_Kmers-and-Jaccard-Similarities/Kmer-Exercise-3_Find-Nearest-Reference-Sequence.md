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

1. Attach the `Definitions.html` file as in previous exercises, again using the prompt:
    ```
    I am going to give you some definitions in an attached file; you don't need to respond to them, just learn them. I'm then going to ask you some questions.
    ```

2. Ask Grimoire to write a program named `nearest-reference-sequence.py` that implements the following specifications:
    * Accepts a mandatory integer argument, short-name `-K`, long-name `--kmer-length`
    * Accepts a mandatory FASTA filename argument, short-name `-R`, long-name `--RepSet`
    * Reads in the FASTA RepSet file using BioPython, then splits the "sequence-description" into "sequence-ID" and "sequence-comment" fields, and stores the sequence-ID, sequence-comment, and sequence.
    * Reads a set of FASTA query-sequences from STDIN, and extracts the query-ID, query-comment, and query-sequence
    * For each input query-sequence, reports as a TSV the query-ID and the query-description, the most similar RepSet sequence-ID and RepSet-ID's description, and the Kmer jaccard-similarity.

3. Save Grimoire's program as in previous exercises.

4. Run the program as follows:
```
python3 nearest-reference-sequence.py -K 8 -R ../Data/rep10.seed_proteins.faa Data/test-nearest.faa > Data/test-nearest.jaccard.tab
```
5. You can compare your output to the provided solution with the following command:
```
diff Data/test-nearest.jaccard.tab Solutions/test-nearest.jaccard.solutions.tab
```
`diff` is a command that compares two files. If the files being compared are identical, then `diff` should produce no output; however, if the files are different, then `diff` will list the differences between the two files. So in this case, no result is a good result! :-)

6. Open the output-file `Data/test-nearest.jaccard.tab` in VScode and scroll through the output. Notice that in many cases, the description of the query-sequence (which is the name of the genome the query came from) is similar to the description of the representative sequence --- and that in the cases where the two descriptions are very different, the jaccard-similarity is often small, indicating that the two sequences are only weakly similar. Bioinformaticians use more sophisticated versions of the tool you have just written to determine whether a sequence or genome is one that they have seen before, or something that is different enough that it is possibly a new species.

7. BONUS: Re-run step (4.) with smaller and larger values of `K` to observe how the nearest reference sequence and the jaccard-similarity changes. Note that if `K` is too large, most of the jaccard-similarities will collapse to `0.0`, indicating that large `K` is so "specific" that it is not very "sensitive", as discussed in Kmer-Exercise-2.


## NOTES:

1.) Need to create the file `rep10.seqs.faa` and the test-sequences.

2.) What is the best way to facilitate the student checking their results? Do we include the sequence-description in the test-file and have Grimoire output the description, so that the student can compare the input description and the RepGenSet description? (In a real-world problem, the sequences will not come "pre-labeled".)