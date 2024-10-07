########################################################################
#...Prompt that generated this program:
"""
Please write a program named `translate_DNA.py`
that reads a set of FASTA-encoded DNA sequences from STDIN,
standardizes their case, translates each sequence from DNA to protein
using the "Bacterial" genetic code (NCBI transl_table=11),
and writes the translated sequences in FASTA format to STDOUT.
Please be careful to ensure that all three possible START codons
are translated as 'M' when they are in the initial position.
Embedded STOPs should be translated as '*'.
Terminal STOPs should not be translated.
Invalid codons should be translated to 'X'.
"""
########################################################################
#...Pseudocode for this program:
"""
# Define the bacterial genetic code as a dictionary
genetic_code = {
    "TTT": "F", "TTC": "F", ..., "TGG": "W"
}

# Define the set of start codons
start_codons = {"ATG", "GTG", "TTG"}

# Function to translate a DNA sequence into a protein sequence
function translate_sequence(dna_seq):
    # Initialize an empty list for the protein sequence
    protein_seq = []

    # Convert the DNA sequence to uppercase
    dna_seq = dna_seq.upper()

    # Loop over the DNA sequence in steps of 3 (codon length)
    for i from 0 to length(dna_seq) - 2, step 3:
        codon = dna_seq[i:i+3]  # Extract the codon

        # Check if the codon is a start codon and is at the initial position
        if i == 0 and codon in start_codons:
            append "M" to protein_seq
        else if codon in genetic_code:
            amino_acid = genetic_code[codon]

            # If the codon is a stop codon and not terminal, append "*"
            if amino_acid == "*" and i + 3 < length(dna_seq):
                append "*" to protein_seq
            else if amino_acid != "*":
                append amino_acid to protein_seq
        else:
            # For invalid codons, append "X"
            append "X" to protein_seq

    # Return the translated protein sequence as a string
    return join(protein_seq)

# Function to read FASTA input from STDIN, translate, and write output to STDOUT
function read_fasta_and_translate():
    current_sequence = []
    current_header = None

    # Loop over each line from standard input
    for each line in STDIN:
        remove whitespace from line

        # If the line starts with ">", it's a FASTA header
        if line starts with ">":
            # If there's a previous sequence, translate it
            if current_header is not None:
                dna_seq = join(current_sequence)
                protein_seq = translate_sequence(dna_seq)
                print current_header
                print protein_seq

            # Set the new header and reset the sequence accumulator
            current_header = line
            current_sequence = []
        else:
            # Accumulate the sequence lines
            append line to current_sequence

    # Translate and print the last sequence
    if current_header is not None and current_sequence is not empty:
        dna_seq = join(current_sequence)
        protein_seq = translate_sequence(dna_seq)
        print current_header
        print protein_seq

# Main program entry point
if __name__ == "__main__":
    call read_fasta_and_translate()
"""
########################################################################
#...Grimoire-generated code:
import sys
import re

# Bacterial Genetic Code NCBI transl_table=11
genetic_code = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',  # Start codon
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',  # Start codon
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',  # Stop codons
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*',  # Stop codon
    'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W',
}

# Define the start codons for NCBI transl_table=11
start_codons = {'ATG', 'GTG', 'TTG'}

def translate_sequence(dna_seq):
    protein_seq = []
    dna_seq = dna_seq.upper()
    
    # Translate the DNA sequence in triplets (codons)
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        
        if i == 0 and codon in start_codons:  # Translate start codons as 'M'
            protein_seq.append('M')
        elif codon in genetic_code:
            amino_acid = genetic_code[codon]
            if amino_acid == '*' and i + 3 < len(dna_seq):  # Skip terminal stop codons
                protein_seq.append('*')
            elif amino_acid != '*':  # Non-terminal stops
                protein_seq.append(amino_acid)
        else:
            protein_seq.append('X')  # Invalid codon
        
    return ''.join(protein_seq)

def read_fasta_and_translate():
    current_sequence = []
    current_header = None
    
    for line in sys.stdin:
        line = line.strip()
        
        if line.startswith(">"):  # Header line
            if current_header:  # If there is a previous sequence, translate it
                dna_seq = ''.join(current_sequence)
                protein_seq = translate_sequence(dna_seq)
                print(current_header)
                print(protein_seq)
            current_header = line
            current_sequence = []
        else:  # Sequence line
            current_sequence.append(line)
    
    # Translate the last sequence in the file
    if current_header and current_sequence:
        dna_seq = ''.join(current_sequence)
        protein_seq = translate_sequence(dna_seq)
        print(current_header)
        print(protein_seq)

if __name__ == "__main__":
    read_fasta_and_translate()
