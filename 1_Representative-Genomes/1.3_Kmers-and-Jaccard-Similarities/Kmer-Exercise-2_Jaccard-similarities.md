#### Kmer Exercise 2 - Jaccard Similarities

Objective: Learn about the concept of "Jaccard similarities", and how they are used to compare sequences.

Comparing and contrasting sequences is the main way that scientists find patterns within DNA and protien sequences. The Jaccard Similarity is a formula that we can use to determine how similar two sets (or groups of kmers) are from each other. If two sequences are very similar, they might be from the same species or genus. If they are very different, chances are they are not the same species or genus. This lesson allows you to explore the Jaccard Similarities Formula within the context of kmers.

#### Materials: 

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)
``` 
FIG-Bioinformatics-Course/
    Basics.html
    1 - Representative Genomes/
        Data/Sample1.fasta
        1.3 - Kmers and Jaccard Similarities/
            Code/kmer_jaccard_similarities.py
            Kmer Exercise 2 - Jaccard similarities.md
            Solutions/kmer_jaccard_similarities_solution.py
```

#### Exercise:

0. Prepare Grimoire for this session as in the previous exercise by once again entering the prompt:

```
I am going to give you some definitions in an attached file;
you don't need to respond to them, just learn them.
I'm then going to ask you some questions.
```

and then attaching the file "Basics.html" using the "paperclip" icon before clicking the "Send message" ("Up-arrow") icon.

1. Ask Grimoire "What is a 'Jaccard Similarity', and how is it used?"

2. Ask Grimoire "How would I compute the jaccard-similarity of two sequences using Kmers?"

3. Ask Grimoire to write a program that takes a Kmer-length and a FASTA filename as command-line arguments, and then prints the jaccard-similarities for all pairs of sequences in the FASTA file to STDOUT.
    * The program should use the output-format "f'{id1}\t{id2}\t{jaccard_sim}'", where 'id1' is the identifier of the first sequence in the comparison, 'id2' is the identifier of the second sequence, and 'jaccard_sim' is the jaccard-similarity between the sequences.

4. Ask Grimoire to give you a detailed explanation of what the output "f-string" format means.

5. Use VScode to save the generated code into the stub-program `Code/kmer_jaccard_similarities.py`, and then use VScode to open a "terminal" window and run the program on the data-file `Sample1.fasta` using '20' as the value of `K`.
    * Grimoire should have shown you how to run the program from a terminal-window, but if it didn't, please ask it to show you how to run the program.

6. Run the same command for the program 'Code/kmer_jaccard_similarities_solution and check your program's output against the output from the solution program.
    * Notice that most of the 20-mer jaccard-similarities are "0.0"; this is not an error --- it means that most of the sequences do not have even a single 20-mer in common and means our sample is diverse.

7. Ask Grimoire to explain the difference between the concepts of "sensitivity" and "specificity", and how they are relevant in bioinformatics and in medicine.


8. Repeat the exercise with `K=10`, and note that now most of the jaccard-similarities are non-zero.
    * This illustrates that there is a tradeoff between "sensitivity" and "specificity": Larger values of `K` yield results that are for more specific (i.e., fewer "hits" with be found, and the pairs of sequences will be "closer" to each other, i.e. "more similar"), but the price to be paid is that a more "specific" value of `K` is less "sensitive", i.e., it may not find anything at all.

## Solution Check instructions:


# Notes:
    Step 3: Can we expand that?
        Add pseudocode step before making the program
        Add bullets to this step for students to understand what is happening behind the scenes
    Step 4: This zooms in to a very specific detail. What is the purpose of explaining this detail?
        This is to tell Grimoire to use the fstring
        Need to add a formating exercise to the TSV unit for output
    Step 5: We have always provided the command line prompt before. Is there a reason we are not including it here?
    Step 8: Can we get away from the i.e. piece? 
        Add to step 7?
        
    Exercise after this needs to give the students a use case for the program. 
        Can we use it to determine taxonomy? define all close genomes to a single genome?
        Exploration into why we chose 8mers for protiens and 20mers for DNA
        Define a use case in which a 0.0 jaccard similarity is good or in which it is bad. Same for 100% similar cases


     
