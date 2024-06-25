#########################################################
#       Insert your Pseudo Code Below                   #
# Use the '#' character at the beginning of a           #
# 'single-line' comment to ensure that it will          #
# not register as code.                                 #
# The three quotes below will protect a "block comment" #                           #
"""
FUNCTION main():
    PARSE command-line argument K AS integer
    INITIALIZE kmer_counts AS empty dictionary
    INITIALIZE kmer_genomes AS empty dictionary

    READ input data from STDIN
    SPLIT input data by newline into lines
    SKIP the first line (header)

    FOR each line IN lines:
        SPLIT line by tab INTO genome_id AND sequence
        FOR each possible Kmer of length K IN sequence:
            EXTRACT Kmer
            INCREMENT kmer_counts[Kmer]
            ADD genome_id TO kmer_genomes[Kmer]

    INITIALIZE hammers AS empty list
    FOR each Kmer IN kmer_counts:
        IF kmer_counts[Kmer] == 1 AND LENGTH(kmer_genomes[Kmer]) == 1:
            ADD (Kmer, FIRST_ELEMENT(kmer_genomes[Kmer])) TO hammers

    PRINT "hammer\tgenome_id"
    FOR each hammer, genome_id IN hammers:
        PRINT hammer + "\t" + genome_id

CALL main()
   
"""
#########################################################
# Replace the 'pass' statement' below with your code

import sys
import argparse
from collections import defaultdict

def main():
    # Parse command-line argument
    parser = argparse.ArgumentParser(description='Find hammers in genome sequences.')
    parser.add_argument('K', type=int, help='Length of the Kmer')
    args = parser.parse_args()
    K = args.K
    
    # Dictionary to store Kmer occurrences
    kmer_counts = defaultdict(int)
    kmer_genomes = defaultdict(set)
    
    # Read input data from STDIN
    input_data = sys.stdin.read().strip().split('\n')
    
    # Skip the header line
    input_data = input_data[1:]
    
    # Process each line in the input data
    for line in input_data:
        genome_id, sequence = line.split('\t')
        for i in range(len(sequence) - K + 1):
            kmer = sequence[i:i+K]
            kmer_counts[kmer] += 1
            kmer_genomes[kmer].add(genome_id)
    
    # Find hammers
    hammers = [(kmer, list(genomes)[0]) for kmer, count in kmer_counts.items() if count == 1 and len(kmer_genomes[kmer]) == 1]
    
    # Print results
    print("hammer\tgenome_id")
    for hammer, genome_id in hammers:
        print(f"{hammer}\t{genome_id}")

if __name__ == "__main__":
    main()

