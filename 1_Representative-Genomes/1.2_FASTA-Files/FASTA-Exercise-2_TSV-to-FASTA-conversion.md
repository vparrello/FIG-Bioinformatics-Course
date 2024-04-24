#### FASTA Exercise 2 - TSV to FASTA conversion

Objectives: 
1. Create a program that converts a 3-column tab-separated file into a FASTA-formatted file.

2. Learn how to "pipe" the output of one commmand-line program to the input of another program, thereby chaining together their functions.

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

* FIG-Bioinformatics-Course/
    * 1_Representative-Genomes/
        * Data/
            * rep10.seqs.tbl 
        * 1.2_FASTA-Files/
            * FASTA-Exercise-2_TSV-to-FASTA-translation.md

#### Exercise: 

In bioinformatics we are often faced with the need to extract data from files and convert them from the format produced by one program to a different format needed by a second program.
In this exercise, you will first create a program that will reformat a 3-column tab-separated file into a FASTA-formatted file. You will then use the program `cmd_tsv_select_columns.py` that you created in `TSV Ex. 3` to extract 3 selected columns from a tab-separated file containing genome and sequence-data, and "pipe" them to the program that converts 3-column data into FASTA format. 