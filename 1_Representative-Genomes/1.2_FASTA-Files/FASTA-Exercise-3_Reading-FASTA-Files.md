#### FASTA Exercise 3 - Reading FASTA files.

Objective: Use Grimoire to write a program that reads and operates on FASTA-formatted sequence-data.

FASTA files are a common format used to store and transmit bioinformatic sequence data. This means that we need to not only understand how these data are represented inside the file, but also have an effective way reading FASTA data and extracting information from it. This exercise creates a rudimentary program that can read a FASTA file, and report a set of summary statistics on the file's contents.

#### Materials: 

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

* FIG-Bioinformatics-Course/
    * 1_Representative-Genomes/
        * Data/
            * Sample1.fasta
        * Code/
            * fasta_reader.py
        * Solutions/
            * fasta_reader_solution.py
        * 1.2_FASTA-Files/
            * FASTA-Exercise-3_Reading-FASTA-Files.md
            

#### Exercise:

1. Ask Grimoire to explain in pseudocode how to read and parse the sequence of records within a FASTA-format file into a list containing "ID", "Description", and "Sequence" data.

2. Ask Grimoire to convert its pseudocode into python code, and then explain to you how the code works line-by-line if it did not already do so.

3. Ask Grimoire to give pseudocode for and then write a python program that will accept a FASTA filename as a command-line argument, read the file, print a tab-separated ID and sequence-length for each record to STDOUT, and finally print the number of sequences read and the average length of the set of sequences to STDERR.

    * WARNING: Grimoire sometimes forgets that the "sequence-ID" stops at the first "whitespace character", and that the rest of the record-header after the sequence-ID is considered a "description" of the sequence, not a part of its ID, so if you see additional text after the sequence-ID in the TSV output-table, you will need to remind Grimoire that it should not include the sequence-description in its TSV output.

4. Copy the pseudocode to the clipboard, and then use VScode to save the pseudocode to the pseudocode section of the program template `Code/fasta_reader.py` 

5. Copy the program itself to the clipboard, and then use VScode to save the program to the code section of the template `Code/fasta_reader.py`.

6. Ask Grimoire to explain line-by-line how the program works.

7. Use the VScode terminal-window to run the program `Code/fasta_reader.py` on the file `Data/Sample1.fasta`.
    * Grimoire should have given you an example of how to run the program, but if it didn't, please ask it to show you how. 
    * If you are working within a Windows environment, you should specify that in your questions.

* Bonus Exercise 1: In FASTA-Ex-2, you extracted a FASTA-file from a TSV-file, which you saved as `rep10.seed_protein.faa`. Run `fasta_reader.py` on `rep10.seed_protein.faa`, and save the TSV-output as follows:
```
python3 Code/fasta_reader.py < rep10.seed_protein.faa > rep10.seed_protein.lengths.tab
```

* Bonus Exercise 2: Repeat the procedure from FASTA-Ex-2, but this time extract the field `seed_dna` instead of `seed_protein`:
```
python3 ../Code/cmd_tsv_select_columns.py genome_id genome_name seed_dna < ../Data/rep10.seqs.tbl | python3 3col_to_fasta > rep10.seed_dna.fna
```
(Note that this time we have used the file-extension ".fna" instead of ".faa", because we are extracting "nucleic-acid" data instead of "amino-acid" data.)

Now, run `fasta_reader.py` on `rep10.seed_dna.fna`:
```
python3 Code/fasta_reader.py rep10.seed_dna.fna > rep10.seed_dna.lengths.tab
```
Compare the files `rep10.seed_protein.lengths.tab` to `rep10.seed_dna.lengths.tab`. You should notice that for each protein-sequence, the corresponding DNA-sequence is 3 times longer; this is because 3 DNA characters translate to a single amino-acid character. We will explore the concept of "sequence translation" further in the next exercise.

## Solution Check instructions:
Use the solution code provided to check your results.
