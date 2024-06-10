# Kmer Exercise 4 - Protein vs DNA Kmers

Objectives:
1. Understand how the Jaccard-similarity computed using a DNA-sequence's protein translation correlates with the Jaccard-similarity computed using the DNA itself

2. Get a better feel for the tradeoffs between "sensitivity" and "specificity"



## Materials: 

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

```
FIG-Bioinformatics-Course
├── Definitions.html
├── Data/
│   └── rep10.seqs.tbl
└── 1_Representative-Genomes
    └── 1.3_Kmers-and-Jaccard-Similarities
        ├── Kmer-Exercise-4_Proteinvs-DNA-Kmers.md (you are here)
        └── Solutions
            └── (TBD)
```

## Exercises

1. Attach the `Definitions.html` file as in previous exercises, again using the prompt:
    ```
    I am going to give you some definitions in an attached file; you don't need to respond to them, just learn them. I'm then going to ask you some questions.
    ```

    2. Ask Grimoire to write a program named `protein_vs_dna_jaccard.py` meeting the following specifications:
    * Mandatory integer-valued protein Kmer-length argument, short-name `-P`, long-name `--protK`.
    * Mandatory integer-valued DNA Kmer-length argument, short-name `-D`, long-name `--dnaK`.
    * Mandatory TSV-filename argument, short-name `-d`, long-name `--data`.
    * Each line of the TSV-file contains data from a genome. The program should skip the header-line, then read the `genome_id`, `seed_protein`, and `seed_dna` columns.
    * Foreach pair of genomes, calculate the protein jaccard-similarity using the data in column `seed_protein` with protein Kmer-length`protK`, and the DNA jaccard-similarity using the data in column `seed_dna` with DNA Kmer-length `dnaK`; count the number of pairs processed, the number of nonzero protein similarities, and the number of nonzero DNA similarities.
    * Plot the similarity-pairs with the DNA jaccard-similarity on the X-axis and the protein jaccard-similarity on the Y-axis.
    * Print the number of nonzero protein-similarities and the number of nonzero DNA-similarities to STDERR, then exit.

