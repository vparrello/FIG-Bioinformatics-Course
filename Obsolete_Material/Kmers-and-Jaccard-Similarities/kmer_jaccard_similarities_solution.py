import sys
from Bio import SeqIO

def generate_kmers(sequence, k):
    return {str(sequence[i:i+k]) for i in range(len(sequence) - k + 1)}

def jaccard_similarity(set1, set2):
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    return len(intersection) / len(union) if union else 0

def main(kmer_length, fasta_filename):
    sequences = list(SeqIO.parse(fasta_filename, "fasta"))
    kmer_length = int(kmer_length)

    # Calculate Jaccard similarity for all pairs
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            kmers1 = generate_kmers(sequences[i].seq, kmer_length)
            kmers2 = generate_kmers(sequences[j].seq, kmer_length)
            jaccard_sim = jaccard_similarity(kmers1, kmers2)
            print(f"{sequences[i].id}\t{sequences[j].id}\t{jaccard_sim}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <kmer_length> <fasta_filename>")
    else:
        kmer_length, fasta_filename = sys.argv[1], sys.argv[2]
        main(kmer_length, fasta_filename)
