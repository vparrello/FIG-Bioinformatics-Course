import sys
import argparse
from Bio import SeqIO

def is_valid_dna(kmer):
    valid_dna_chars = set("acgt")
    return all(char in valid_dna_chars for char in kmer)

def is_valid_protein(kmer):
    valid_protein_chars = set("ACDEFGHIKLMNPQRSTVWY")
    return all(char in valid_protein_chars for char in kmer)

def process_fasta(Kmer_length, sequence_type):
    num_sequences = 0
    num_valid_kmers = 0
    num_invalid_kmers = 0

    for record in SeqIO.parse(sys.stdin, "fasta"):
        num_sequences += 1
        sequence = str(record.seq)
        seq_id = record.id
        
        if sequence_type == 'dna':
            sequence = sequence.lower()
        elif sequence_type == 'protein':
            sequence = sequence.upper()
        
        for i in range(len(sequence) - Kmer_length + 1):
            kmer = sequence[i:i + Kmer_length]
            if sequence_type == 'dna' and not is_valid_dna(kmer):
                print(f"{kmer}\t{seq_id}", file=sys.stderr)
                num_invalid_kmers += 1
            elif sequence_type == 'protein' and not is_valid_protein(kmer):
                print(f"{kmer}\t{seq_id}", file=sys.stderr)
                num_invalid_kmers += 1
            else:
                print(f"{kmer}\t{seq_id}")
                num_valid_kmers += 1

    print(f"Number of sequences read: {num_sequences}", file=sys.stderr)
    print(f"Number of valid kmers written to STDOUT: {num_valid_kmers}", file=sys.stderr)
    print(f"Number of invalid kmers found: {num_invalid_kmers}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description="Process FASTA file and extract Kmers")
    parser.add_argument('-K', '--Kmer', type=int, required=True, help="Kmer length")
    parser.add_argument('-t', '--type', choices=['dna', 'protein'], help="Type of sequence (dna or protein)")
    args = parser.parse_args()

    process_fasta(args.Kmer, args.type)

if __name__ == "__main__":
    main()
