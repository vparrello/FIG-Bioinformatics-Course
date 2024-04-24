#### FASTA Exercise 2 - Working with FASTA files.

Objective: Use Grimoire to write a program that reads and operates on FASTA-formatted sequence-data.

FASTA files are the bulk of the data that scientists use to explore the bioinformatic scientific space. This means that we need to not only understand the data that is inside of the file, but also have an extremely efficient way to sort through that data and pull the needed information from it. This exercise creates a rudimentary program that can read a FASTA file efficiently and quickly.

#### Materials: 

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)
``` 
FIG-Bioinformatics-Course/
    1 - Representative Genomes/
        Data/Sample1.fasta
        Code/fasta_reader.py
        Code/fasta_reader_solution.py
        1.2 - FASTA Files/
            FASTA Exercise 1 - What is FASTA format.md
            
```

#### Exercise:

1. Ask Grimoire to write a python program that will read a named FASTA file without using BioPython, then print the tab-separated ID and length for each sequence to STDOUT, and finally the number of sequences read and the average length of the set of sequences to STDERR. 
    *WARNING: Grimoire sometimes forgets that the "sequence-ID" stops at the first "whitespace character", and that the rest of the record-header after the sequence-ID is considered a "description" of the sequence, not part of its ID.

2. Copy the pseudocode to the clipboard, and then use VScode to save the pseudocode to the pseudocode section of the program template `Code/fasta_reader.py` 

3. Copy the program itself to the clipboard, and then use VScode to save the program to the code section of the template `Code/fasta_reader.py`.

4. Ask Grimoire to explain line-by-line how the program works.

5. Use the VScode terminal-window to run the program on the file `Data/Sample1.fasta`.
    * Grimoire should have given you an example of how to run the program, but if it didn't, please ask it to show you how. 
    * If you are working within a Windows environment, you should specify that in your questions.

## Solution Check instructions:
Use the solution code provided to check your results.


# NOTES: 
    Bullet list for directory tree
    Add file extension why the file extensions matter to the user and not code 
    How the two alphabets are different between protein and DNA
    More samples for Exercise 3
    Needs solution code for fasta_reader - in the bin directory on Gordon's computer
    Have a protein and dna sample for data

    Exercises following:
        Read the rep200.seqs.tbl and create a fasta file from that list
        Extract sequences from tsv_header
        Need an exercise that explores BioPython and its uses

        