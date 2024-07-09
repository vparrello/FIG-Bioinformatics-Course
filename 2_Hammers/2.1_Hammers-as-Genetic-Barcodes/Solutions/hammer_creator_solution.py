#########################################################
#       Insert your Pseudo Code Below                   #
# Use the '#' character at the beginning of a           #
# 'single-line' comment to ensure that it will          #
# not register as code.                                 #
# The three quotes below will protect a "block comment" #                           #
"""
# define function to parse arguments
function parse_arguments():
    create an ArgumentParser object
    add argument '-K' with type int and required=True
    return parsed arguments

# define function to read and parse sequences
function read_input():
    genome_sequences = empty list
    input_lines = read all input from stdin and split by newline
    for each line in input_lines starting from the second line:
        split line by tab
        genome_id = first part
        sequence = second part, stripped of any whitespace
        append (genome_id, sequence) tuple to genome_sequences
    return genome_sequences

# define function to map Kmers to genomes
function extract_kmers(genome_sequences, K):
    kmer_genomes = dictionary with default value as an empty set
    for each (genome_id, sequence) in genome_sequences:
        for i from 0 to length of sequence minus K + 1:
            kmer = substring of sequence from i to i + K
            add genome_id to the set in kmer_genomes[kmer]
    return kmer_genomes

# define function to find the hammers
function find_hammers(kmer_genomes):
    hammers = empty dictionary
    for each kmer in kmer_genomes:
        if length of kmer_genomes[kmer] is 1:
            genome_id = first element in list(kmer_genomes[kmer])
            hammers[kmer] = genome_id
    return hammers

# define function for main program
function main():
    args = parse_arguments()
    K = args.K
    genome_sequences = read_input()
    kmer_genomes = extract_kmers(genome_sequences, K)
    hammers = find_hammers(kmer_genomes)
    
    print "hammer\tgenome_id"
    for each hammer in hammers:
        print hammer + "\t" + hammers[hammer]

# entry point of program
if __name__ == "__main__":
    main()

"""
#########################################################
# Replace the 'pass' statement' below with your code

import sys
import argparse
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description="Find Hammers in DNA sequences")
    parser.add_argument('-K', type=int, required=True, help="Kmer length")
    return parser.parse_args()

def read_input():
    genome_sequences = []
    input_lines = sys.stdin.read().strip().split('\n')
    for line in input_lines[1:]:  # Skip header line
        parts = line.strip().split('\t')
        genome_id = parts[0]
        sequence = parts[1].strip()
        genome_sequences.append((genome_id, sequence))
    return genome_sequences

def extract_kmers(genome_sequences, K):
    kmer_genomes = defaultdict(set)
    for genome_id, sequence in genome_sequences:
        for i in range(len(sequence) - K + 1):
            kmer = sequence[i:i + K]
            kmer_genomes[kmer].add(genome_id)
    return kmer_genomes

def find_hammers(kmer_genomes):
    hammers = {}
    for kmer in kmer_genomes:
        if len(kmer_genomes[kmer]) == 1:
            genome_id = list(kmer_genomes[kmer])[0]  # Use list() and [0] to get the single genome_id
            hammers[kmer] = genome_id
    return hammers

def main():
    args = parse_arguments()
    K = args.K
    genome_sequences = read_input()
    kmer_genomes = extract_kmers(genome_sequences, K)
    hammers = find_hammers(kmer_genomes)
    
    print("hammer\tgenome_id")
    for hammer in hammers:
        print(f"{hammer}\t{hammers[hammer]}")

if __name__ == "__main__":
    main()
