import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import os
from collections import defaultdict

def kmer_set(sequence, k):
    return {str(sequence[i:i+k]) for i in range(len(sequence) - k + 1)}

def jaccard_similarity(set1, set2):
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection / union if union != 0 else 0

def main():
    parser = argparse.ArgumentParser(description="Compute K-mer Jaccard similarity")
    parser.add_argument('-K', '--kmer-length', type=int, required=True, help='Length of K-mers')
    parser.add_argument('-R', '--RepSet', type=str, required=True, help='FASTA file of the representative set')
    args = parser.parse_args()

    k = args.kmer_length
    rep_set_file = args.RepSet

    if not os.path.isfile(rep_set_file):
        print(f"Error: Representative set file '{rep_set_file}' does not exist.")
        sys.exit(1)

    # Read the RepSet file
    rep_set = []
    for record in SeqIO.parse(rep_set_file, "fasta"):
        rep_set.append((record.id, record.description, str(record.seq)))

    # Read the query sequences from STDIN
    query_sequences = []
    for record in SeqIO.parse(sys.stdin, "fasta"):
        query_sequences.append((record.id, record.description, str(record.seq)))

    # Process each query sequence
    results = []
    for query_id, query_description, query_seq in query_sequences:
        query_kmers = kmer_set(query_seq, k)
        best_match = None
        best_similarity = -1

        for rep_id, rep_description, rep_seq in rep_set:
            rep_kmers = kmer_set(rep_seq, k)
            similarity = jaccard_similarity(query_kmers, rep_kmers)
            if similarity > best_similarity:
                best_similarity = similarity
                best_match = (rep_id, rep_description)

        results.append((query_id, query_description, best_match[0], best_match[1], best_similarity))

    # Print results as TSV
    for result in results:
        print("\t".join(map(str, result)))

if __name__ == "__main__":
    main()
