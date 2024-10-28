########################################################################
#...Prompt that generated this program:
"""
Please write a program named `translate-DNA_intermediate.py`
that reads a set of FASTA-encoded DNA sequences from STDIN,
standardizes their case, translates each sequence from DNA to protein
using the "Bacterial" genetic code (NCBI transl_table=11),
and writes the translated sequences in FASTA format to STDOUT.
Please be careful to ensure that all three possible START codons
are translated as 'M' when they are in the initial position.
STOP codons should terminate translation.
Invalid codons should be translated to 'X'.
"""
########################################################################
#...Pseudocode for this program:
"""
# Define the bacterial genetic code using a dictionary
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'M', 
    ... # other codon mappings
    'TAA': '*', 'TAG': '*', 'TGA': '*',  # Stop codons
}

# Function to translate a DNA sequence to a protein sequence
function TRANSLATE_DNA(DNA_SEQUENCE):
    PROTEIN_SEQUENCE = empty list
    SEQUENCE_LENGTH = length of DNA_SEQUENCE
    
    for every 3 characters in DNA_SEQUENCE (index i = 0, step by 3):
        CODON = substring of DNA_SEQUENCE from i to i+3
        if CODON length is less than 3:
            break
        
        if CODON is a start codon (ATG, GTG, TTG) and index i == 0:
            append 'M' to PROTEIN_SEQUENCE
        else if CODON exists in GENETIC_CODE:
            AMINO_ACID = GENETIC_CODE[CODON]
            if AMINO_ACID is '*':  # Stop codon
                break
            append AMINO_ACID to PROTEIN_SEQUENCE
        else:
            append 'X' to PROTEIN_SEQUENCE  # Invalid codon
            
    return PROTEIN_SEQUENCE as string

# Function to read FASTA sequences from stdin
function READ_FASTA():
    SEQUENCES = empty dictionary
    CURRENT_SEQUENCE = empty list
    CURRENT_HEADER = None
    
    for each LINE from standard input:
        if LINE starts with '>':  # Header line
            if CURRENT_HEADER is not None:
                add CURRENT_HEADER: UPPERCASE(CURRENT_SEQUENCE) to SEQUENCES
            CURRENT_HEADER = LINE
            CURRENT_SEQUENCE = empty list
        else:
            append LINE to CURRENT_SEQUENCE
    
    if CURRENT_HEADER is not None:
        add CURRENT_HEADER: UPPERCASE(CURRENT_SEQUENCE) to SEQUENCES
    
    return SEQUENCES

# Function to write FASTA sequences to stdout
function WRITE_FASTA(SEQUENCES):
    for each HEADER and PROTEIN_SEQUENCE in SEQUENCES:
        print HEADER
        for each 60 characters in PROTEIN_SEQUENCE:
            print the substring

# Main function to handle the translation workflow
function MAIN():
    # Read DNA sequences from input
    DNA_SEQUENCES = READ_FASTA()
    
    TRANSLATED_SEQUENCES = empty dictionary
    
    # Translate each sequence
    for each HEADER, DNA_SEQUENCE in DNA_SEQUENCES:
        TRANSLATED_SEQUENCES[HEADER] = TRANSLATE_DNA(DNA_SEQUENCE)
    
    # Output translated sequences
    WRITE_FASTA(TRANSLATED_SEQUENCES)

# Run the main function if script is executed
if __name__ == '__main__':
    MAIN()
"""
########################################################################
#...Grimoire-generated code:
import sys

# Define the bacterial genetic code (NCBI transl_table=11)
genetic_code = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',  # Start codon ATG
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'M',  # GTG start codon
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*',  # Stop codon
    'TAG': '*',  # Stop codon
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*',  # Stop codon
    'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# Function to translate a DNA sequence to a protein sequence
def translate_dna(dna_seq):
    protein_seq = []
    seq_len = len(dna_seq)
    
    for i in range(0, seq_len, 3):
        codon = dna_seq[i:i+3]
        if len(codon) < 3:
            break
        
        # Translate the start codons ATG, GTG, TTG to 'M' if in the initial position
        if i == 0 and codon in ['ATG', 'GTG', 'TTG']:
            protein_seq.append('M')
        elif codon in genetic_code:
            aa = genetic_code[codon]
            if aa == '*':  # Stop codon
                break
            protein_seq.append(aa)
        else:
            protein_seq.append('X')  # Invalid codon
        
    return ''.join(protein_seq)

# Function to read FASTA sequences from stdin
def read_fasta():
    seqs = {}
    current_seq = []
    current_header = None

    for line in sys.stdin:
        line = line.strip()

        if line.startswith('>'):  # Header line
            if current_header:
                seqs[current_header] = ''.join(current_seq).upper()
            current_header = line
            current_seq = []
        else:
            current_seq.append(line)
    
    if current_header:  # Add the last sequence
        seqs[current_header] = ''.join(current_seq).upper()

    return seqs

# Function to write FASTA sequences to stdout
def write_fasta(seqs):
    for header, protein_seq in seqs.items():
        print(header)
        for i in range(0, len(protein_seq), 60):
            print(protein_seq[i:i+60])

# Main function to read, translate, and output sequences
def main():
    dna_seqs = read_fasta()
    translated_seqs = {}

    for header, dna_seq in dna_seqs.items():
        translated_seqs[header] = translate_dna(dna_seq)
    
    write_fasta(translated_seqs)

if __name__ == '__main__':
    main()
