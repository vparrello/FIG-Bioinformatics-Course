import sys
import gzip
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


def read_hammer_file(file_path):
    """
    Reads a tab-separated-value file, possibly gzipped, and loads the first and second columns into a dictionary.
    """
    hammer_dict = {}
    with (gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r')) as file:
        header = next(file)  # Skip the header
        for line in file:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            hammer, feature_id = parts[0], parts[1]
            hammer_dict[hammer] = feature_id
    hammer_length = len(next(iter(hammer_dict))) if hammer_dict else 0
    return hammer_length, hammer_dict


def find_kmers(sequence, k):
    """
    Finds all kmers of length k in a sequence and its reverse complement.
    """
    kmers = []
    sequence_rc = str(Seq(sequence).reverse_complement())
    for i in range(len(sequence) - k + 1):
        kmers.append(sequence[i:i + k])
        kmers.append(sequence_rc[i:i + k])
    return kmers


def count_hammers(fasta_file, hammer_length, hammer_dict):
    """
    Counts occurrences of hammers (kmers) in the DNA sequences from a FASTA file.
    """
    hammer_counts = defaultdict(int)
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        kmers = find_kmers(sequence, hammer_length)
        for kmer in kmers:
            if kmer in hammer_dict:
                hammer_counts[kmer] += 1
    return hammer_counts


def main():
    if len(sys.argv) != 2:
        print("Usage: python report_hammer_counts.py <hammer_file>", file=sys.stderr)
        sys.exit(1)

    hammer_file = sys.argv[1]
    hammer_length, hammer_dict = read_hammer_file(hammer_file)
    
    if hammer_length == 0:
        print("Error: No hammers found in the input file.", file=sys.stderr)
        sys.exit(1)
    
    hammer_counts = count_hammers(sys.stdin, hammer_length, hammer_dict)
    
    sorted_hammers = sorted(
        ((hammer, count, hammer_dict[hammer]) for hammer, count in hammer_counts.items() if count > 1),
        key=lambda x: -x[1]
    )
    
    for hammer, count, feature_id in sorted_hammers:
        print(f"{hammer}\t{count}\t{feature_id}")


if __name__ == "__main__":
    main()
