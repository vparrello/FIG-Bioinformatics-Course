# FASTA Exercise 4 - DNA to Protein Translation

Objective: Learn how to translate DNA sequences into protein sequences.

We have mention that there are three fundamental types of biological sequences (DNA, RNA, and protein), and that each type of sequence has its own "alphabet". We will now learn how cells translate genes encoded in DNA to proteins that are composed of amino-acids.


## Materials: 

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

```
FIG-Bioinformatics-Course/
├── 1_Representative-Genomes/
│   ├── 1.2_FASTA-Files/
│   │   └── FASTA-Exercise-4_DNA-to-Protein-Translation.md (you are here)
│   └── Solutions/
│       └── translate-DNA_solution.py        
├── Code/
│   └── translate-DNA.py
└── Data/
    └── Sample1.fasta
```

## Exercises:

1. With only a few exceptions (most notably during the replication of "retroviruses" such as HIV), information in cells flows in one direction: From DNA, to RNA, to Proteins. This one-way flow of molecular information has become known as the "Fundamental Dogma of Molecular Biology". Ask Grimoire to explain this "Fundamental Dogma" to you in more detail. If Grimoire uses terms that are unfamiliar to you, ask it to explain those terms as well.

2. A protein-encoding DNA sequence is structured as a sequence of 3-character groups known as "codons". Each codon encodes a single amino-acid. Since DNA uses a 4-character alphabet, there are 4x4x4 = 64 possible codons, which map into 20 possible amino-acids plus a "STOP" signal indicating the end of the gene. Ask Grimoire to show you the mapping from DNA to amino-acids in Bacteria.

3. Ask Grimoire to explain in pseudocode how to translate from a DNA sequence into a protein sequence.

4. Ask Grimoire to write a program named `translate_DNA.py` that reads a set of FASTA-encoded DNA sequences from STDIN, translates each sequence from DNA to protein, and writes the translated sequences in FASTA format to STDOUT. Use VScode to save the program to the `Code/` directory for this exercise.

5. Run the program on the file `Data/Sample1.fasta`:
    ```
    python3 Code/translate_DNA.py < Data/Sample1.fasta > Data/Sample1.faa
    ```
Note that we are recommending a file-extension of `.faa` to remind you that the translation is an Amino-Acid FASTA file.

6. Check your output translation against the file `Solutions/Sample1.faa`.

7. Bonus: Real biology contains many complications, and one of them is that some lifeforms use slightly different versions of the genetic code, and slightly different methods of transcription and translation.
    * Ask Grimoire how transcription differs in Bacteria, Archaea, and Eukaryotes.
    * Ask Grimoire how translation differs in Bacteria, Archaea, and Eukaryotes.
    * Ask Grimoire to explain the differences between the various different genetic codes.
