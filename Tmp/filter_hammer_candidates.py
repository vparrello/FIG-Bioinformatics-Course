import os
import sys
import argparse
from Bio import SeqIO
from collections import defaultdict

def read_hammer_candidates(input_stream):
    """Reads the hammer candidate file from STDIN."""
    header = input_stream.readline().strip()
    hammer_dict = {}
    hammer_length = None

    for line in input_stream:
        fields = line.strip().split("\t")
        hammer = fields[0]
        if hammer_length is None:
            hammer_length = len(hammer)
        elif len(hammer) != hammer_length:
            raise ValueError("Inconsistent hammer lengths in input")
        hammer_dict[hammer] = line.strip()

    return len(hammer_dict), hammer_length, header, hammer_dict

def process_contigs(directory, hammer_length, hammer_dict):
    """Process contigs and filter hammer candidates."""
    candidate_counts = defaultdict(int)
    files = [os.path.join(directory, f) for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

    for file in files:
        for record in SeqIO.parse(file, "fasta"):
            seq = str(record.seq).lower()
            reverse_complement = str(record.seq.reverse_complement()).lower()

            for i in range(len(seq) - hammer_length + 1):
                kmer = seq[i:i + hammer_length]
                if kmer in hammer_dict:
                    candidate_counts[kmer] += 1
                    if candidate_counts[kmer] > 1:
                        del hammer_dict[kmer]

            for i in range(len(reverse_complement) - hammer_length + 1):
                kmer = reverse_complement[i:i + hammer_length]
                if kmer in hammer_dict:
                    candidate_counts[kmer] += 1
                    if candidate_counts[kmer] > 1:
                        del hammer_dict[kmer]

def main():
    parser = argparse.ArgumentParser(description="Filter hammer candidates.")
    parser.add_argument('-C', '--contigs-directory', required=True, help='Directory containing contigs files')
    args = parser.parse_args()

    # Read hammer candidates from STDIN
    num_candidates, hammer_length, header, hammer_dict = read_hammer_candidates(sys.stdin)

    # Process contigs directory
    process_contigs(args.contigs_directory, hammer_length, hammer_dict)

    # Output accepted hammers
    print(header)
    for hammer in sorted(hammer_dict.keys()):
        print(hammer_dict[hammer])

    # Output statistics to STDERR
    print(f"Candidates read: {num_candidates}", file=sys.stderr)
    print(f"Hammers accepted: {len(hammer_dict)}", file=sys.stderr)

if __name__ == "__main__":
    main()
