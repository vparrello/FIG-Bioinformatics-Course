#########################################################
#  Prompt that generated this solution:
"""
Please write a python program named 'hammer_creator.py' that will:

        * Accept as '-K' or '--kmer_length' a mandatory Kmer-length command-line argument.
	
        * Accept as '-i' or '--input' the name of a mandatory FASTA-formated sequence file.

	* Accept as '-o' or '--output' the name of a mandatory tab-separated-value hammer-file.

	* The sequence-identifiers of the input FASTA file are "feature_IDs",
	or "fids" for short. A "feature_ID" has the form 'fig|X.Y.peg.Z',
	where the portions 'fig|' and '.peg.' are literal strings,
	'X', 'Y', and 'Z' are integers, and the subpattern 'X.Y' encodes
	the 'genome_id' that the feature came from. For each feature_ID,
	you will need to be able to extract the 'genome_id' portion of that feature_ID;
	you are allowed to use "regular expressions".
	
        * Find all of the Kmers that occur exactly once in exactly one genome;
	these Kmers are the "Hammers".
	
        * Print a two column tab-separated table of the hammers and the feature_ID
	that the hammer was found in, with column-headers of "hammer" and "fid".
	
	* Print how long the program took to run to STDERR, then exit.
"""
#########################################################
# Pseudocode generated by this solution prompt:
"""
START main:
    START timer
    
    PARSE command-line arguments
        - Kmer length (-K or --kmer_length)
        - Input FASTA file (-i or --input)
        - Output file (-o or --output)

    CALL find_hammers with input file and Kmer length
        STORE returned kmer_counts and fid_to_kmers

    CALL write_hammers with kmer_counts, fid_to_kmers, and output file
    
    PRINT program execution time to STDERR
    END main

START find_hammers(input_file, kmer_length):
    CREATE kmer_counts as a dictionary (default values: dictionary with int values)
    CREATE fid_to_kmers as a dictionary (default values: list)

    FOR each record in the input FASTA file:
        SET feature_id to record identifier
        SET genome_id to extracted genome ID from feature_id (use regex)
        SET sequence to the record’s sequence string

        FOR each substring of length Kmer from the sequence:
            SET kmer to the substring
            INCREMENT count of kmer in kmer_counts[genome_id]
            ADD kmer to the list for feature_id in fid_to_kmers

    RETURN kmer_counts, fid_to_kmers
    END find_hammers

START write_hammers(kmer_counts, fid_to_kmers, output_file):
    OPEN output_file for writing
    WRITE header ("hammer", "fid") to output file

    FOR each fid in fid_to_kmers:
        FOR each kmer in the list of kmers for this fid:
            IF kmer appears exactly once in exactly one genome (check kmer_counts):
                WRITE kmer and fid to the output file

    END write_hammers
"""
#########################################################
# Code generated by this solution prompt:
import sys
import argparse
import re
from Bio import SeqIO
import time
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description="Find hammers from a FASTA file.")
    parser.add_argument('-K', '--kmer_length', type=int, required=True, help='Length of the Kmers to extract')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file for hammers')
    return parser.parse_args()

def extract_genome_id(feature_id):
    # Extract 'X.Y' from 'fig|X.Y.peg.Z'
    match = re.match(r'fig\|(\d+\.\d+)\.peg\.\d+', feature_id)
    if match:
        return match.group(1)
    return None

def find_hammers(input_file, kmer_length):
    kmer_counts = defaultdict(lambda: defaultdict(int))
    fid_to_kmers = defaultdict(list)
    
    # Process FASTA file
    for record in SeqIO.parse(input_file, "fasta"):
        feature_id = record.id
        genome_id = extract_genome_id(feature_id)
        sequence = str(record.seq)
        
        # Extract kmers of the specified length
        for i in range(len(sequence) - kmer_length + 1):
            kmer = sequence[i:i+kmer_length]
            kmer_counts[kmer][genome_id] += 1
            fid_to_kmers[feature_id].append(kmer)

    return kmer_counts, fid_to_kmers

def write_hammers(kmer_counts, fid_to_kmers, output_file):
    with open(output_file, 'w') as out_file:
        out_file.write("hammer\tfid\n")
        for fid, kmers in fid_to_kmers.items():
            for kmer in kmers:
                # A hammer is a kmer that occurs exactly once in exactly one genome
                if len(kmer_counts[kmer]) == 1 and list(kmer_counts[kmer].values())[0] == 1:
                    out_file.write(f"{kmer}\t{fid}\n")

def main():
    start_time = time.time()

    # Parse command-line arguments
    args = parse_arguments()

    # Find all kmers and store their occurrences
    kmer_counts, fid_to_kmers = find_hammers(args.input, args.kmer_length)

    # Write the hammers to the output file
    write_hammers(kmer_counts, fid_to_kmers, args.output)

    # Print execution time
    elapsed_time = time.time() - start_time
    print(f"Program finished in {elapsed_time:.2f} seconds", file=sys.stderr)

if __name__ == "__main__":
    main()