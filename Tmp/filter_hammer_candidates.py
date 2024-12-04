import os
import sys
import gzip
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


def read_hammer_file(file_path):
    """
    Reads a tab-separated file, possibly gzipped, and loads the first column as keys and the entire line as values.
    """
    hammer_dict = {}
    with (gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r')) as file:
        header = next(file).strip()  # Save the header line
        for line in file:
            parts = line.strip().split("\t")
            if len(parts) < 1:
                continue
            hammer = parts[0]
            hammer_dict[hammer] = line.strip()
    hammer_length = len(next(iter(hammer_dict))) if hammer_dict else 0
    return hammer_length, header, hammer_dict


def find_kmers(sequence, k):
    """
    Finds all kmers of length k in a sequence and its reverse complement, converting sequence to lowercase.
    """
    sequence = sequence.lower()
    sequence_rc = str(Seq(sequence).reverse_complement())
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmers.append(sequence[i:i + k])
        kmers.append(sequence_rc[i:i + k])
    return kmers


def count_kmers_in_contigs(contigs_directory, hammer_length):
    """
    Counts occurrences of kmers in all FASTA files within the given directory.
    """
    kmer_counts = defaultdict(int)
    for file_name in os.listdir(contigs_directory):
        file_path = os.path.join(contigs_directory, file_name)
        if not os.path.isfile(file_path):
            continue
        for record in SeqIO.parse(file_path, "fasta"):
            sequence = str(record.seq)
            kmers = find_kmers(sequence, hammer_length)
            for kmer in kmers:
                kmer_counts[kmer] += 1
    return kmer_counts


def main():
    parser = argparse.ArgumentParser(description="Filter hammer candidates based on contig sequences.")
    parser.add_argument("-H", "--hammer-file", required=True, help="Path to the hammer file (can be gzipped).")
    parser.add_argument("-C", "--contigs-directory", required=True, help="Path to the directory containing FASTA contig files.")
    args = parser.parse_args()

    hammer_length, header_line, hammer_dict = read_hammer_file(args.hammer_file)
    if hammer_length == 0:
        print("Error: No hammers found in the input file.", file=sys.stderr)
        sys.exit(1)

    kmer_counts = count_kmers_in_contigs(args.contigs_directory, hammer_length)

    # Print header line
    print(header_line)

    # Filter and print hammers
    printed_count = 0
    for hammer, value in hammer_dict.items():
        if kmer_counts[hammer] == 1:
            print(value)
            printed_count += 1

    # Print summary to STDERR
    print(f"Candidates read:  {len(hammer_dict)}", file=sys.stderr)
    print(f"Hammers accepted: {printed_count}", file=sys.stderr)


if __name__ == "__main__":
    main()
