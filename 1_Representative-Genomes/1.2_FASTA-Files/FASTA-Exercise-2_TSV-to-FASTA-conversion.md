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

1. In bioinformatics we are often faced with the need to extract data from files and convert them from the format produced by one program to a different format needed by some other program.
In this exercise, you will first create a program that will reformat a 3-column tab-separated file into a FASTA-formatted file. You will then use the program `cmd_tsv_select_columns.py` that you created in `TSV Ex. 2` to extract 3 selected columns from a tab-separated file containing genome and sequence-data, and "pipe" them to the program that converts 3-column data into FASTA format.

If you use your script `tsv_headers.py` on the file `Data/rep10.seqs.tbl`, you should see something like this:
```
% python3 ~/bin/tsv_headers.py ../Data/rep10.seqs.tbl 
Field names in the TSV file are:
genome_id
genome_name
seed_protein
ssu_rna
seed_dna
```
The above means that this is a 5-column TSV file; the first two columns are a genome-ID and its human-readable biological name, while the last three columns contain sequence-data.

2. Ask Grimoire to write a program named `3col_to_fasta.py` that:
    * accepts TSV-data from STDIN,
    * skips the header-line,
    * treats the first column as a sequence-ID, the second column as a sequence-description, and the third column as sequence-data,
    * writes the data to STDOUT in FASTA format.

Use VScode to save the program into the `Code/` directory as usual.

3. Use `cmd_tsv_select_columns.py` from TSV-Ex-2 to select the columns 'genome_id', 'genome_name', and 'seed_protein', and "pipe" the output to `3col_to_fasta.py`:
```
python3 ../Code/cmd_tsv_select_columns.py genome_id genome_name seed_protein < ../Data/rep10.seqs.tbl | python3 3col_to_fasta > rep10.seed_protein.faa
```

* NOTE: The above should all be entered on a single command-line, even though your browser has probably split it across several lines.

The results should match the corresponding file in the `Solutions/` directory.

* NOTE: If Grimoire did not use `BioPython`, you should instead check the alternate solution file.