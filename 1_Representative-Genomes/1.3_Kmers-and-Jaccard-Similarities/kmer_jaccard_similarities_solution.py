import sys
from itertools import combinations
from Bio import SeqIO

def kmer_set(sequence, k):
    """Return a set of k-mers from the sequence."""
    return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}

def jaccard_similarity(set1, set2):
    """Calculate the Jaccard similarity between two sets."""
    intersection = set1 & set2
    union = set1 | set2
    return len(intersection), len(intersection) / len(union)

def main(k, fasta_file):
    # Parse the FASTA file
    sequences = [(record.id, str(record.seq)) for record in SeqIO.parse(fasta_file, "fasta")]
    
    # Calculate k-mer sets for each sequence
    kmer_sets = {seq_id: kmer_set(seq, k) for seq_id, seq in sequences}
    
    # Compare all pairs of sequences
    for (id1, seq1), (id2, seq2) in combinations(sequences, 2):
        set1, set2 = kmer_sets[id1], kmer_sets[id2]
        num_common, jaccard_sim = jaccard_similarity(set1, set2)
        print(f'{id1}\t{id2}\t{num_common}\t{jaccard_sim:.6f}')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python kmer_jaccard_similarities.py <k> <fasta_file>")
        sys.exit(1)
    
    k = int(sys.argv[1])
    fasta_file = sys.argv[2]
    main(k, fasta_file)