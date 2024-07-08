import argparse
from Bio import SeqIO
from collections import defaultdict
import sys

"""
Parse Command-Line Arguments:

Use argparse to handle the -K/--kmer-length and -R/--RepSet arguments.
Read the RepSet FASTA File:

Use BioPython to parse the FASTA file.
Split the sequence descriptions into sequence-ID and sequence-comment.
Store the sequence-ID, sequence-comment, and sequence in a dictionary.
Read Query Sequences from STDIN:

Use BioPython to read the query sequences from standard input.
Extract query-ID, query-comment, and query-sequence.
Compute Kmers:

For each sequence (both RepSet and query), compute all possible kmers of length K.
Compare Query Sequences to RepSet Sequences:

For each query sequence, compute the Jaccard similarity with each RepSet sequence.
Find the RepSet sequence with the highest Jaccard similarity.
Write the TSV Report:

For each query sequence, write the query-ID, query-description, most similar RepSet sequence-ID, RepSet description, number of common kmers, and Jaccard similarity to a TSV file.
"""
def parse_args():
    parser = argparse.ArgumentParser(description="Find the nearest reference sequence by Kmer similarity.")
    parser.add_argument('-K', '--kmer-length', type=int, required=True, help='Length of Kmers')
    parser.add_argument('-R', '--RepSet', type=str, required=True, help='FASTA file of the reference set')
    return parser.parse_args()

def read_fasta(filename):
    sequences = {}
    for record in SeqIO.parse(filename, "fasta"):
        seq_id, seq_comment = record.description.split(None, 1)
        sequences[seq_id] = {
            'comment': seq_comment,
            'sequence': str(record.seq)
        }
    return sequences

def read_query_sequences():
    sequences = {}
    for record in SeqIO.parse(sys.stdin, "fasta"):
        seq_id, seq_comment = record.description.split(None, 1)
        sequences[seq_id] = {
            'comment': seq_comment,
            'sequence': str(record.seq)
        }
    return sequences

def compute_kmers(sequence, k):
    kmers = set()
    for i in range(len(sequence) - k + 1):
        kmers.add(sequence[i:i + k])
    return kmers

def jaccard_similarity(set1, set2):
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection, intersection / union if union > 0 else 0

def main():
    args = parse_args()

    # Read RepSet sequences
    rep_set = read_fasta(args.RepSet)
    rep_kmers = {seq_id: compute_kmers(data['sequence'], args.kmer_length) for seq_id, data in rep_set.items()}

    # Read query sequences from stdin
    query_set = read_query_sequences()

    # Prepare report
    report = ["Query-ID\tQuery-Description\tRepSet-ID\tRepSet-Description\tCommon-Kmers\tJaccard-Similarity"]
    
    for query_id, query_data in query_set.items():
        query_kmers = compute_kmers(query_data['sequence'], args.kmer_length)
        
        best_match = None
        best_similarity = 0
        best_common_kmers = 0
        
        for rep_id, rep_data in rep_set.items():
            common_kmers, similarity = jaccard_similarity(query_kmers, rep_kmers[rep_id])
            if similarity > best_similarity:
                best_match = rep_id
                best_similarity = similarity
                best_common_kmers = common_kmers
        
        if best_match:
            report.append(f"{query_id}\t{query_data['comment']}\t{best_match}\t{rep_set[best_match]['comment']}\t{best_common_kmers}\t{best_similarity:.4f}")
    
    # Output report
    for line in report:
        print(line)

if __name__ == "__main__":
    main()
