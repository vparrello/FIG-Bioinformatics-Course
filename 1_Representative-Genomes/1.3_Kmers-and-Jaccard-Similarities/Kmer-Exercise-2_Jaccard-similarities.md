# Kmer Exercise 2 - Jaccard Similarities

Objective: Learn about the concept of "Jaccard similarities", and how they are used to compare sequences.

Comparing and contrasting sequences is the main way that scientists find patterns within DNA and protein sequences. The Jaccard Similarity is a formula that we can use to determine how similar two sets (or groups of Kmers) are from each other. If two sequences that perform then same role are very similar, there is good chance that they from two closely-related species or genera. But if they are very different, the chances are high that they they are from two distantly-related species or genera. This lesson allows you to explore using the Jaccard Similarities Formula to estimate sequence-similarity using kmers.

## Materials: 

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

```
FIG-Bioinformatics-Course/
├── Definitions.html
└── 1_Representative-Genomes/
    ├── Code/
    │   └── kmer_jaccard_similarities.py
    ├── Data/
    │   └── Sample1.fasta
    └── 1.3_Kmers-and-Jaccard-Similarities/
        ├── Kmer-Exercise-2_Jaccard-similarities.md  (you are here)
        └── Solutions/
            └── kmer_jaccard_similarities_solution.py
```

## Exercises:

1. Prepare Grimoire for this session as in the previous exercise by once again entering the prompt:

```
I am going to give you some definitions in an attached file;
you don't need to respond to them, just learn them.
I'm then going to ask you some questions.
```

and again attaching the file `Definitions.html` using the "paperclip" icon before clicking the "Send message" ("Up-arrow") icon.

2. Ask Grimoire "What is a 'Jaccard Similarity', and how is it used?"

3. Ask Grimoire "How would I compute the jaccard-similarity of two sequences using Kmers?"

4. Ask Grimoire to write a program that:

* Takes a Kmer-length as a command-line argument,

* Uses BioPython to read a FASTA-formatted file from STDIN,

* Computes the number of Kmers in common and the jaccard-similarities for all pairs of sequences in the FASTA file

* Prints the sequence-IDs, number of Kmers in common, and the jaccard-similarities to STDOUT in tab-separated format if the number of Kmers in common is not zero.

* Prints the number of pairs and the number of pairs with nonzero similarity to STDERR.

5. Ask Grimoire to translate the program into pseudocode.

6. Use VScode to save the generated pseudocode and code into the stub-program `Code/kmer_jaccard_similarities.py`,

7. Use VScode to open a "terminal" window and run the program on the data-file `Sample1.fasta` using '20' as the value of `K`:

```
python Code/kmer_jaccard_similarities.py 20 < Data/Sample1.fasta > Sample1.out
```

* Notice that the number of pairs with nopnzero similarity is less than the number of pairs of sequences; this is not an error, it means that many of the sequences do not have even a single 20-mer in common, which implies that our sample is diverse.

* Run the program `1_Representative-Genomes/1.3_Kmers-and-Jaccard-Similarities/Solutions/kmer_jaccard_similarities_solution.py/kmer_jaccard_similarities_solution.py` using the same command-line arguments, and check your program's output against the output from the solution program using the `diff` command.

8. Repeat the exercise with `K=10`, and note that now most of the jaccard-similarities are non-zero.
    * This exercise illustrates that there is a tradeoff between "sensitivity" and "specificity": Larger values of `K` yield results that are far more specific (i.e., fewer "hits" will be found, and that the pairs that are found will be "closer" to each other, i.e. "more similar"), but the price to be paid for a more "specific" value of `K` is that it is also less "sensitive", so that if `K` is too large, the search may not find anything at all.

9. Ask Grimoire to explain the difference between the concepts of "sensitivity" and "specificity", and how they are relevant in bioinformatics and in medicine.

10. In bioinformatics, we will often specify similarity-thresholds in terms of the number of Kmers that two sequences must have in common rather than in terms of their Jaccard-similarity.

11. Ask Grimoire to write a program `num_kmers_vs_jaccard.py` that takes a Kmer-length and a FASTA filename as command-line arguments, and then for all pairs of sequences plots the number of Kmers in common between vs their jaccard-similarity. Save this program as in exercise (5.), and then run it for various values of `K` to get a feel for how the number of Kmers in common and the jaccard-similarities vary with `K`. As in the previous exercise, you should see that larger `K` will result in fewer pairs of sequences with Kmers in common (also referred to as fewer "hits"), and also that fewer Kmers in common correlates with a smaller jaccard-similarity, so that larger `K` means smaller sensitivity (smaller scores) but also larger specificity (fewer hits).

## Solution Check instructions:

(TBD)


# Notes:
    Step 4: Can we expand that?
        Add pseudocode step before making the program
        Add bullets to this step for students to understand what is happening behind the scenes
    Step 5: This zooms in to a very specific detail. What is the purpose of explaining this detail?
        This is to tell Grimoire to use the fstring
        Need to add a formating exercise to the TSV unit for output
    Step 6: We have always provided the command line prompt before. Is there a reason we are not including it here?
        
    Exercise after this needs to give the students a use case for the program. 
        Can we use it to determine taxonomy? define all close genomes to a single genome?
        Exploration into why we chose 8mers for protiens and 20mers for DNA
        Define a use case in which a 0.0 jaccard similarity is good or in which it is bad. Same for 100% similar cases

Steps from here to Hammers:
    Create a Repgen Set
    Separate the PHEs sequence of each repgen (Do we want to have more than 1 role?)
    
     
