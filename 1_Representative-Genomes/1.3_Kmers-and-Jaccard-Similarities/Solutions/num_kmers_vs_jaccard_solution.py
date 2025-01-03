import sys
import argparse
from Bio import SeqIO
import matplotlib.pyplot as plt
from itertools import combinations

def generate_kmers(sequence, k):
    """Generate k-mers from a given sequence."""
    return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}

def jaccard_similarity(set1, set2):
    """Calculate the Jaccard similarity between two sets."""
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection / union if union != 0 else 0

def common_kmers(set1, set2):
    """Find the number of k-mers in common between two sets."""
    return len(set1 & set2)

def main():
    parser = argparse.ArgumentParser(description="Usage: python num_kmers_vs_jaccard.py -k <kmer-length> -f <fasta-file>")
    parser.add_argument("-k", "--kmer_length", type=int, help="Length of the k-mer")
    parser.add_argument("-f", "--fasta_file", type=str, help="Path to the FASTA file")
    args = parser.parse_args()

    k = args.kmer_length
    fasta_file = args.fasta_file

    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))

    kmer_sets = [generate_kmers(seq, k) for seq in sequences]

    common_kmers_list = []
    jaccard_similarity_list = []

    for (seq1, kmers1), (seq2, kmers2) in combinations(zip(sequences, kmer_sets), 2):
        common_kmers_count = common_kmers(kmers1, kmers2)
        if common_kmers_count > 0:
            jaccard_sim = jaccard_similarity(kmers1, kmers2)
            common_kmers_list.append(common_kmers_count)
            jaccard_similarity_list.append(jaccard_sim)

    plt.scatter(jaccard_similarity_list, common_kmers_list)
    plt.xlabel("Jaccard Similarity")
    plt.ylabel("Number of K-mers in Common")
    plt.title(f"Number of K-mers in Common vs Jaccard Similarity (k={k})")
    plt.show()

if __name__ == "__main__":
    main()
