#########################################################
#       Insert your Pseudo Code Below                   #
# Use the '#' character at the beginning of a           #
# 'single-line' comment to ensure that it will          #
# not register as code.                                 #
# The three quotes below will protect a "block comment" #                           #
"""
# Define function to parse command line arguments
function parse_arguments():
    Create argument parser
    Add '-K' argument with integer type and required flag
    Return parsed arguments

# Define function to read input from STDIN
function read_input():
    Read all input data from STDIN and split by newlines
    Skip the header line
    Initialize empty list for genome sequences
    For each line in input data:
        Split line by tab character
        If line has at least two columns:
            Extract genome_id and sequence from first two columns
            Append (genome_id, sequence) pair to genome sequences list
    Return genome sequences list

# Define function to find hammers in genome sequences
function find_hammers(genome_sequences, k):
    Initialize empty dictionary for kmer counts
    Initialize empty dictionary for kmer genomes
    For each (genome_id, sequence) pair in genome sequences:
        For each position in sequence from 0 to length(sequence) - k:
            Extract kmer of length k from current position
            Increment kmer count in kmer counts dictionary
            Add genome_id to kmer's genome set in kmer genomes dictionary
    Initialize empty list for hammers
    For each kmer, genome_ids pair in kmer genomes dictionary:
        If kmer count is 1 and genome_ids set contains exactly one element:
            Add (kmer, genome_id) to hammers list
    Return hammers list

# Main function
function main():
    Parse arguments and get kmer length k
    Read genome sequences from input
    Find hammers using genome sequences and kmer length
    Print header "hammer\tgenome_id"
    For each (hammer, genome_id) in hammers:
        Print hammer and genome_id separated by tab

# Entry point of the program
if __name__ == '__main__':
    Call main function
   
"""
#########################################################
# Replace the 'pass' statement' below with your code

import sys
import argparse
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description='Find hammers in genome sequences.')
    parser.add_argument('-K', type=int, required=True, help='Length of Kmer')
    return parser.parse_args()

def read_input():
    input_data = sys.stdin.read().strip().split('\n')
    header = input_data[0]  # Skip the header
    genome_sequences = []
    for line in input_data[1:]:
        parts = line.split('\t')
        if len(parts) >= 2:
            genome_id = parts[0]
            sequence = parts[1]
            genome_sequences.append((genome_id, sequence))
    return genome_sequences

def find_hammers(genome_sequences, k):
    kmer_counts = defaultdict(int)
    kmer_genomes = defaultdict(set)
    
    for genome_id, sequence in genome_sequences:
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            kmer_counts[kmer] += 1
            kmer_genomes[kmer].add(genome_id)

    hammers = [(kmer, list(genome_ids)[0]) for kmer, genome_ids in kmer_genomes.items() if kmer_counts[kmer] == 1 and len(genome_ids) == 1]
    
    return hammers

def main():
    args = parse_arguments()
    k = args.K
    genome_sequences = read_input()
    hammers = find_hammers(genome_sequences, k)
    
    print("hammer\tgenome_id")
    for hammer, genome_id in hammers:
        print(f"{hammer}\t{genome_id}")

if __name__ == '__main__':
    main()
