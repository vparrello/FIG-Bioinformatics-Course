########################################################################
# Prompt that generated 'build_representative_set_solution.py'
"""
I will now give you a description of the command-line interface and
program called `build_representative_set.py` that will compute a
"set of representative sequences" (RepGen set) using the "Stingy Addition"
algorithm as defined in the uploaded definitions-file.

The script should accept the following mandatory command-line
arguments, in both long-form and the specified short-form:

* Kmer-length  (integer, short argument-name '-k')
* Sim          (similarity threshold, short argument-name '-s')
* input FASTA-file   (filename, short argument-name '-f')
* output RepSeq-file  (filename, short argument-name '-r')

The measure of similarity to be used is "Number of Kmers in common".

The program should use BioPython to read the FASTA-file.
The first nonwhitespace portion of each FASTA identifier is a "feature-ID",
while the remainder is the "genome name".
Feature-IDs have the format 'fig|X.Y.peg.Z', where 'X', 'Y', and 'Z' are integers,
and the portions 'fig|' and '.peg.' are literal substrings, not variables.
The subpattern 'X.Y' within the feature-ID is a "genome-ID";
you may extract this subpattern using a regular expression.
Return a dictionary that maps feature-IDs to genome-IDs,
a dictionary that maps feature-IDs to genome-names,
and a list of (feature-ID, sequence) pairs.

The main body of the program should construct a subset of the input sequence list
that satisfies the provided definition of a "Representative Set" (RepGen set).

When you are done, please write out the representative set to the RepSeq-file in FASTA format,
where the FASTA identifiers have the form "genome-ID genome-name".
"""
########################################################################
# Pseudocode for 'build_representative_set_solution.py'
"""DEFINE FUNCTION parse_fasta(fasta_file):
    INITIALIZE empty dictionary feature_to_genome
    INITIALIZE empty dictionary feature_to_name
    INITIALIZE empty list sequences
    DEFINE genome_id_pattern as a regular expression for 'fig|X.Y.peg.Z'
    
    FOR each record in fasta_file:
        SET header as the record's description
        SET seq as the record's sequence
        
        SET feature_id as the first portion of header (split by space)
        
        IF feature_id matches genome_id_pattern:
            EXTRACT genome_id from the match
        ELSE:
            RAISE error for invalid feature_id format
        
        SET genome_name as the rest of the header (after feature_id)
        
        ADD genome_id to feature_to_genome[feature_id]
        ADD genome_name to feature_to_name[feature_id]
        APPEND (feature_id, seq) to sequences
    
    RETURN feature_to_genome, feature_to_name, sequences


DEFINE FUNCTION compute_kmers(sequence, k):
    INITIALIZE empty set kmer_set
    FOR each i from 0 to (length of sequence - k):
        ADD subsequence of length k from index i to kmer_set
    RETURN kmer_set


DEFINE FUNCTION stingy_addition(sequences, feature_to_genome, k, sim_threshold):
    INITIALIZE empty list repgen_set
    INITIALIZE empty list repgen_kmer_sets
    
    FOR each (feature_id, seq) in sequences:
        COMPUTE current_kmers from seq using compute_kmers(seq, k)
        
        INITIALIZE is_similar as False
        FOR each existing_kmers in repgen_kmer_sets:
            SET common_kmers as intersection of current_kmers and existing_kmers
            IF size of common_kmers >= sim_threshold:
                SET is_similar as True
                BREAK out of loop
        
        IF is_similar is False:
            APPEND (feature_id, seq) to repgen_set
            APPEND current_kmers to repgen_kmer_sets
    
    RETURN repgen_set


DEFINE FUNCTION write_fasta(repgen_set, feature_to_genome, feature_to_name, output_file):
    OPEN output_file for writing
    
    FOR each (feature_id, seq) in repgen_set:
        SET genome_id as feature_to_genome[feature_id]
        SET genome_name as feature_to_name[feature_id]
        WRITE ">genome_id genome_name" to output_file
        WRITE sequence to output_file
    
    CLOSE output_file


DEFINE FUNCTION main():
    PARSE command-line arguments for kmer length, sim threshold, input fasta file, output fasta file
    
    CALL parse_fasta with input fasta file
    CALL stingy_addition with parsed sequences, feature_to_genome, kmer length, and sim threshold
    
    CALL write_fasta with repgen_set, feature_to_genome, feature_to_name, and output fasta file


IF __name__ == "__main__":
    CALL main()
"""
########################################################################
# Code for 'build_representative_set.py'
import argparse
import re
from Bio import SeqIO
from collections import defaultdict

def parse_fasta(fasta_file):
    """
    Parse the input FASTA file using BioPython and extract feature-IDs, genome-IDs, and genome names.
    Return:
    - A dictionary that maps feature-IDs to genome-IDs.
    - A dictionary that maps feature-IDs to genome-names.
    - A list of (feature-ID, sequence) tuples.
    """
    feature_to_genome = {}
    feature_to_name = {}
    sequences = []
    
    genome_id_pattern = re.compile(r'fig\|(\d+\.\d+)\.peg\.\d+')

    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        seq = str(record.seq)
        
        # Extract the feature-ID (the first non-whitespace portion of the header)
        feature_id = header.split()[0]
        
        # Extract the genome-ID using the regular expression
        genome_id_match = genome_id_pattern.search(feature_id)
        if genome_id_match:
            genome_id = genome_id_match.group(1)
        else:
            raise ValueError(f"Invalid feature-ID format in {feature_id}")
        
        # Extract the genome name (everything after the feature-ID)
        genome_name = header.split(maxsplit=1)[1]
        
        # Store the mappings
        feature_to_genome[feature_id] = genome_id
        feature_to_name[feature_id] = genome_name
        sequences.append((feature_id, seq))

    return feature_to_genome, feature_to_name, sequences

def compute_kmers(sequence, k):
    """
    Given a sequence, return the set of K-mers of length k.
    """
    return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}

def stingy_addition(sequences, feature_to_genome, k, sim_threshold):
    """
    Implement the 'Stingy Addition' algorithm to compute a representative set.
    - sequences: List of (feature-ID, sequence) tuples.
    - feature_to_genome: Dictionary mapping feature-IDs to genome-IDs.
    - k: Kmer length.
    - sim_threshold: Similarity threshold (number of Kmers in common).
    
    Return the representative set of sequences.
    """
    repgen_set = []
    repgen_kmer_sets = []  # Store Kmer sets for the sequences in the repgen set

    for feature_id, seq in sequences:
        # Compute the K-mers of the current sequence
        current_kmers = compute_kmers(seq, k)

        # Check if it is similar to any genome already in the repgen set
        is_similar = False
        for existing_kmers in repgen_kmer_sets:
            common_kmers = current_kmers.intersection(existing_kmers)
            if len(common_kmers) >= sim_threshold:
                is_similar = True
                break

        # If no similar genome is found, add the current sequence to the repgen set
        if not is_similar:
            repgen_set.append((feature_id, seq))
            repgen_kmer_sets.append(current_kmers)

    return repgen_set

def write_fasta(repgen_set, feature_to_genome, feature_to_name, output_file):
    """
    Write the representative set to a FASTA file.
    The format of each FASTA identifier should be 'genome-ID genome-name'.
    """
    with open(output_file, 'w') as f:
        for feature_id, seq in repgen_set:
            genome_id = feature_to_genome[feature_id]
            genome_name = feature_to_name[feature_id]
            f.write(f">{genome_id} {genome_name}\n{seq}\n")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Compute a set of representative sequences using the Stingy Addition algorithm.")
    parser.add_argument('-k', '--kmer', type=int, required=True, help='Kmer length')
    parser.add_argument('-s', '--sim', type=int, required=True, help='Similarity threshold (number of Kmers in common)')
    parser.add_argument('-f', '--fasta', type=str, required=True, help='Input FASTA file')
    parser.add_argument('-r', '--repseq', type=str, required=True, help='Output representative sequence file')
    
    args = parser.parse_args()

    # Parse the input FASTA file
    feature_to_genome, feature_to_name, sequences = parse_fasta(args.fasta)
    
    # Compute the representative set using the Stingy Addition algorithm
    repgen_set = stingy_addition(sequences, feature_to_genome, args.kmer, args.sim)
    
    # Write the representative set to the output file
    write_fasta(repgen_set, feature_to_genome, feature_to_name, args.repseq)

if __name__ == "__main__":
    main()
