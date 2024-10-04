########################################################################
#...Prompt that generated this program:
"""
Please write a program named `translate_DNA.py`
that reads a set of FASTA-encoded DNA sequences from STDIN,
translates each sequence from DNA to protein using the "Bacterial"
genetic code (NCBI transl_table=11), and writes the translated sequences
in FASTA format to STDOUT.
Please be careful to ensure that all three possible START codons
are translated as 'M' when they are in the initial position.
Embedded STOPs should be translated as '*'.
Terminal STOPs should not be translated.
"""
########################################################################
#...Pseudocode for this program:
"""
# Define a dictionary `GENETIC_CODE` containing codons mapped to their corresponding amino acids
# Include start codons and stop codons as per the bacterial genetic code (NCBI Transl Table 11)
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', ... , 'TGA': '*', 'ATG': 'M', 'GTG': 'V', 'TTG': 'L', ...
}

# Define sets for start and stop codons
START_CODONS = {'ATG', 'GTG', 'TTG'}
STOP_CODONS = {'TAA', 'TAG', 'TGA'}

# Define a function `parse_fasta(input_stream)` that:
#   - Reads input line by line
#   - Identifies headers (lines starting with '>')
#   - Aggregates the DNA sequence corresponding to each header
#   - Returns a dictionary with headers as keys and their respective sequences as values
function parse_fasta(input_stream):
    Initialize an empty dictionary `sequences`
    Set `current_header` to None
    Initialize an empty list `current_seq`
    
    for each line in `input_stream`:
        Strip the line of whitespace
        if the line starts with '>':
            if `current_header` is not None:
                Add `current_header` and `current_seq` (joined) to `sequences`
            Set `current_header` to the line
            Reset `current_seq` to an empty list
        else:
            Append the line to `current_seq`
    
    if `current_header` is not None:
        Add `current_header` and `current_seq` (joined) to `sequences`
    
    return `sequences`

# Define a function `translate_dna_to_protein(dna_sequence)` that:
#   - Reads the DNA sequence in groups of 3 (codons)
#   - Translates each codon to its corresponding amino acid using the `GENETIC_CODE`
#   - Translates the first codon if it is a start codon to 'M'
#   - Stops translation at the first stop codon but translates embedded stop codons as '*'
#   - Returns the translated protein sequence, excluding the terminal stop codon
function translate_dna_to_protein(dna_sequence):
    Initialize an empty list `protein_sequence`
    
    for each codon (i.e., 3-nucleotide chunk) in `dna_sequence`:
        if codon length is less than 3, break
        
        if codon is in `START_CODONS` and it is the first codon:
            Append 'M' to `protein_sequence`
        elif codon is in `STOP_CODONS`:
            Append '*' to `protein_sequence`
            break  # Stop translation at the first stop codon
        else:
            Find the corresponding amino acid in `GENETIC_CODE` for the `codon`
            Append the amino acid to `protein_sequence`
    
    Join the `protein_sequence` list into a string
    Remove any trailing '*' from the sequence (i.e., terminal stop codon)
    return the joined `protein_sequence`

# Define the main function that:
#   - Reads the DNA sequences in FASTA format from stdin using `parse_fasta()`
#   - Translates each sequence using `translate_dna_to_protein()`
#   - Outputs the translated protein sequences in FASTA format, wrapping the output at 60 characters
function main():
    Call `parse_fasta()` to get `dna_sequences` from stdin
    
    for each sequence in `dna_sequences`:
        Call `translate_dna_to_protein(sequence)` to get the protein sequence
        Print the corresponding header
        Print the protein sequence, wrapped to 60 characters per line

# If this script is executed as the main program, call the `main()` function
if __name__ == "__main__":
    Call `main()`
"""
########################################################################
#...Grimoire-generated code:
import sys

# Genetic code for the bacterial code (NCBI Transl Table 11)
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',  # Start codon
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',  # GTG also start codon
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',  # Stop codons
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*',  # Stop codon
    'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

START_CODONS = {'ATG', 'GTG', 'TTG'}
STOP_CODONS = {'TAA', 'TAG', 'TGA'}

def parse_fasta(input_stream):
    """Parse FASTA format sequences from input stream."""
    sequences = {}
    current_header = None
    current_seq = []

    for line in input_stream:
        line = line.strip()
        if line.startswith(">"):
            if current_header:
                sequences[current_header] = ''.join(current_seq)
            current_header = line
            current_seq = []
        else:
            current_seq.append(line)
    
    if current_header:
        sequences[current_header] = ''.join(current_seq)

    return sequences

def translate_dna_to_protein(dna_sequence):
    """Translate a DNA sequence into a protein sequence using the bacterial genetic code."""
    protein_sequence = []
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        if len(codon) < 3:
            break  # Skip incomplete codons at the end
        
        if codon in START_CODONS and i == 0:
            protein_sequence.append('M')
        elif codon in STOP_CODONS:
            protein_sequence.append('*')
        else:
            amino_acid = GENETIC_CODE.get(codon, 'X')  # X for unknown codons
            protein_sequence.append(amino_acid)

        if codon in STOP_CODONS:
            break  # Stop translation at the first stop codon

    return ''.join(protein_sequence).rstrip('*')  # Remove terminal stop

def main():
    # Read DNA sequences from stdin (FASTA format)
    dna_sequences = parse_fasta(sys.stdin)
    
    # Translate each sequence and write the result to stdout in FASTA format
    for header, sequence in dna_sequences.items():
        protein = translate_dna_to_protein(sequence)
        print(header)
        for i in range(0, len(protein), 60):
            print(protein[i:i+60])

if __name__ == "__main__":
    main()
