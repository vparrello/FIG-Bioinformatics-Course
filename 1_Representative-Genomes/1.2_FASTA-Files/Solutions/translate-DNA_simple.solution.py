########################################################################
#...Prompt that generated this program:
"""
Please write a program named translate_DNA.py
that reads a set of FASTA-encoded DNA sequences from STDIN,
standardizes sequences to upper-case, translates each sequence
from DNA to protein, and writes the translated sequences in FASTA format
to STDOUT. Invalid codons should be translated to 'X'.
"""
########################################################################
#...Pseudocode for this program:
"""
# Define the genetic code dictionary mapping DNA codons to amino acids
GENETIC_CODE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*', 'TGA': '*'
}

# Function: translate_dna
# Purpose: Translate a DNA sequence to a protein sequence
Function translate_dna(sequence):
    Convert sequence to uppercase
    Initialize an empty list for protein sequence
    For each triplet of nucleotides (codon) in the sequence:
        If codon is found in GENETIC_CODE dictionary:
            Append corresponding amino acid to protein list
        Else:
            Append 'X' (for invalid codon) to protein list
    Return the protein sequence as a string

# Function: process_fasta
# Purpose: Parse FASTA input and process sequences
Function process_fasta(input_stream):
    Initialize an empty string for sequence
    Initialize an empty string for header
    For each line in input_stream:
        Strip whitespace from the line
        If line starts with '>':  # FASTA header line
            If sequence is not empty:
                Print current header
                Print translated protein sequence using translate_dna(sequence)
            Update header with new line
            Reset sequence to empty
        Else:  # Sequence line
            Append line to sequence
    If sequence is not empty after the loop:
        Print last header
        Print translated protein using translate_dna(sequence)

# Main logic
If the script is executed as the main program:
    Call process_fasta with input from standard input (STDIN)
"""
########################################################################
#...Program generated by Grimoire:
import sys

# Define the standard genetic code as a dictionary
GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGA':'*'
}

# Function to translate a DNA sequence into a protein sequence
def translate_dna(sequence):
    sequence = sequence.upper()  # Convert sequence to uppercase
    protein = []
    # Process the sequence 3 nucleotides (1 codon) at a time
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon in GENETIC_CODE:
            protein.append(GENETIC_CODE[codon])
        else:
            protein.append('X')  # Invalid codons are translated to 'X'
    return ''.join(protein)

# Function to parse FASTA input and translate DNA sequences
def process_fasta(stdin):
    sequence = ''
    header = ''
    for line in stdin:
        line = line.strip()
        if line.startswith('>'):  # FASTA header line
            if sequence:
                # Translate the previous sequence
                print(header)
                print(translate_dna(sequence))
            header = line  # Save new header
            sequence = ''  # Reset sequence
        else:
            sequence += line  # Append sequence lines
    # Translate the last sequence in the file
    if sequence:
        print(header)
        print(translate_dna(sequence))

if __name__ == "__main__":
    process_fasta(sys.stdin)
