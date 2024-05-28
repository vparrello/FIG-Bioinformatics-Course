import argparse
import sys
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Find most similar RepSet sequence based on K-mer Jaccard similarity.")
    parser.add_argument('-K', '--kmer-length', type=int, required=True, help='Length of K-mers.')
    parser.add_argument('-R', '--RepSet', type=str, required=True, help='FASTA file of representative sequences.')
    return parser.parse_args()

def read_fasta(filename):
    sequences = {}
    for record in SeqIO.parse(filename, "fasta"):
        sequences[record.id] = (record.description, str(record.seq))
    return sequences

def kmer_set(sequence, k):
    return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}

def jaccard_similarity(set1, set2):
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection / union

def find_most_similar(seq_id, seq, rep_set, k):
    seq_kmers = kmer_set(seq, k)
    max_similarity = 0
    best_match = None
    for rep_id, (rep_desc, rep_seq) in rep_set.items():
        rep_kmers = kmer_set(rep_seq, k)
        similarity = jaccard_similarity(seq_kmers, rep_kmers)
        if similarity > max_similarity:
            max_similarity = similarity
            best_match = (rep_id, rep_desc, similarity)
    return best_match

def main():
    args = parse_args()
    k = args.kmer_length
    rep_set_file = args.RepSet

    rep_set = read_fasta(rep_set_file)

    input_sequences = SeqIO.parse(sys.stdin, "fasta")
    
    print("input_id\tbest_match_id\tbest_match_desc\tjaccard_similarity")
    for record in input_sequences:
        input_id = record.id
        input_seq = str(record.seq)
        best_match_id, best_match_desc, similarity = find_most_similar(input_id, input_seq, rep_set, k)
        print(f"{input_id}\t{best_match_id}\t{best_match_desc}\t{similarity:.6f}")

if __name__ == "__main__":
    main()
