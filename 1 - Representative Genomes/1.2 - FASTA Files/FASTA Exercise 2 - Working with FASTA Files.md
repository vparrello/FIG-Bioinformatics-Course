#### FASTA Exercise 2 - Working with FASTA files.

Objective: Use Grimoire to write a program that reads and operates on FASTA-formatted sequence-data.

#### Materials: 

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)
``` 
FIG-Bioinformatics-Course/
    1 - Representative Genomes/
        1.2 - FastaA Files/
            FASTA Exercise 1 - What is FASTA format.md
            bin/fasta_reader.py
            bin/fasta_reader_solution.py
            Data/Sample1.fasta
```

#### Exercise:

1. Ask Grimoire to write a python program that will read a named FASTA file without using BioPython, then print the tab-separated ID and length for each sequence to STDOUT, and finally the number of sequences read and the average length of the set of sequences to STDERR. 
(WARNING: Grimoire sometimes forgets that the "sequence-ID" stops at the first 
"whitespace character", and that the rest of the record-header after the
sequence-ID is considered a "description" of the sequence, not part of its ID.)

2. Copy the pseudocode to the clipboard, and then use VScode to save the pseudocode to the pseudocode section of the program template `bin/fasta_reader.py` within the `1.2 - FASTA Files` directory.

3. Copy the program itself to the clipboard, and then use VScode to save the program to the code section of the template `bin/fasta_reader.py`.

4. Ask Grimoire to explain line-by-line how the program works.

5. Use the VScode terminal-window to run the program on the file `Data/Sample1.fasta` within the `1.2 - FASTA Files` directory. (Grimoire should have given you an example of how to run the program, but if it didn't, please ask it to show you how. (If you are working within a Windows environment, you should specify that in your questions.))

6. Use the solution code provided to check your results.
