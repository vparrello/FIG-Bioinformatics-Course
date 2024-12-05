import os
import sys
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


def read_hammer_candidates(stdin):
    """Reads the hammer candidate file from STDIN."""
    hammer_dict = {}
    header_line = stdin.readline().strip()
    hammer_length = None

    for line in stdin:
        fields = line.strip().split("\t")
        key = fields[0]
        if hammer_length is None:
            hammer_length = len(key)
        hammer_dict[key] = line.strip()

    return len(hammer_dict), hammer_length, header_line, hammer_dict


def process_contigs(contigs_directory, hammer_length, hammer_dict):
    """Processes the contig files and identifies hammer candidates."""
    candidate_count = defaultdict(int)
    eliminated_count = 0

    for filename in os.listdir(contigs_directory):
        file_path = os.path.join(contigs_directory, filename)

        if not os.path.isfile(file_path):
            continue
        else:
            print(f"Processing file '{filename}'", file=sys.stderr)
        
        for record in SeqIO.parse(file_path, "fasta"):
            seq = str(record.seq).lower()
            rev_comp = str(Seq(seq).reverse_complement()).lower()

            # Generate Kmers for forward and reverse complement
            kmers = set(seq[i:i + hammer_length] for i in range(len(seq) - hammer_length + 1))
            kmers.update(rev_comp[i:i + hammer_length] for i in range(len(rev_comp) - hammer_length + 1))

            for kmer in kmers:
                if kmer in hammer_dict:
                    candidate_count[kmer] += 1
                    if candidate_count[kmer] > 1:
                        del hammer_dict[kmer]
                        eliminated_count += 1

    return eliminated_count


def main():
    parser = argparse.ArgumentParser(description="Filter hammer candidates.")
    parser.add_argument("-C", "--contigs-directory", required=True, help="Directory containing FASTA contig files.")
    args = parser.parse_args()

    # Read hammer candidates from STDIN
    candidates_read, hammer_length, header_line, hammer_dict = read_hammer_candidates(sys.stdin)

    # Process the contigs directory
    eliminated_count = process_contigs(args.contigs_directory, hammer_length, hammer_dict)

    # Print the header line to STDOUT
    print(header_line)
    # Sort and print the accepted hammer candidates
    for key in sorted(hammer_dict.keys()):
        print(hammer_dict[key])

    # Print statistics to STDERR
    accepted_count = len(hammer_dict)
    print(f"Candidates read: {candidates_read}", file=sys.stderr)
    print(f"Candidates eliminated: {eliminated_count}", file=sys.stderr)
    print(f"Hammers accepted: {accepted_count}", file=sys.stderr)


if __name__ == "__main__":
    main()
