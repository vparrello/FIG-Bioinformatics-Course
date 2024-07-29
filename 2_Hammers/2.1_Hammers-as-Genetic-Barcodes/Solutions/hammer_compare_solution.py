#########################################################
#       Insert your Pseudo Code Below                   #
# Use the '#' character at the beginning of a           #
# 'single-line' comment to ensure that it will          #
# not register as code.                                 #
# The three quotes below will protect a "block comment" #                           #
"""
import sys
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

function main()
    if length of sys.argv is not 2
        print "Usage: python hammer_compare.py <hammers_file>"
        exit program with status 1

    hammers_file = sys.argv[1]
    hammers = empty dictionary

    # Read the hammers file
    open hammers_file as file
        create csv reader with tab delimiter
        skip the header line
        for each row in reader
            hammer = row[0]
            genome_id = row[1]
            hammers[hammer] = genome_id

    # Determine Kmer length K
    k = length of any key in hammers

    # Initialize scores dictionary
    scores = defaultdict(int)

    # Read sequences from STDIN
    for each record in SeqIO.parse from stdin in fasta format
        sequence = string representation of record.seq
        rev_sequence = string representation of record.seq.reverse_complement()

        # Scan the sequence and its reverse complement
        for each seq in [sequence, rev_sequence]
            for i from 0 to length of seq - k + 1
                kmer = substring of seq from index i to i + k
                if kmer exists in hammers
                    genome_id = hammers[kmer]
                    increment scores[genome_id] by 1

    # Print results sorted by score in descending order
    sorted_scores = sort scores.items by value in descending order
    for each genome_id, score in sorted_scores
        print genome_id and score separated by tab

if __name__ == "__main__":
    call main()
"""

#########################################################
# Replace the 'pass' statement' below with your code

import sys
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

def main():
    if len(sys.argv) != 2:
        print("Usage: python hammer_compare.py hammers_file < query_genome.fna")
        sys.exit(1)

    hammers_file = sys.argv[1]
    hammers = {}

    # Read hammers file
    with open(hammers_file, newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader)  # Skip header line
        for line_number, row in enumerate(reader, start=1):
            try:
                hammer = row[0]
                genome_id = row[1]
                hammers[hammer] = genome_id
            except IndexError as e:
                print(f"Error '{e}' processing row {line_number}:\n{row}", file=sys.stderr)
                sys.exit(1)

    # Determine Kmer length K
    k = len(next(iter(hammers.keys())))

    # Initialize scores dictionary
    scores = defaultdict(int)

    # Read sequences from STDIN
    for record in SeqIO.parse(sys.stdin, "fasta"):
        sequence = str(record.seq)
        rev_sequence = str(record.seq.reverse_complement())

        # Scan the sequence and its reverse complement
        for seq in [sequence, rev_sequence]:
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i + k]
                if kmer in hammers:
                    genome_id = hammers[kmer]
                    scores[genome_id] += 1

    # Print results sorted by score in descending order
    sorted_scores = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    for genome_id, score in sorted_scores:
        print(f"{genome_id}\t{score}")

if __name__ == "__main__":
    main()

