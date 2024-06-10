import argparse
import pandas as pd
from itertools import combinations
import matplotlib.pyplot as plt
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Calculate and plot Jaccard similarities for genome data.',
        prog='protein_vs_dna_jaccard.py'
    )
    parser.add_argument('-P', '--protK', type=int, required=True, help='Protein K-mer length')
    parser.add_argument('-D', '--dnaK', type=int, required=True, help='DNA K-mer length')
    parser.add_argument('-d', '--data', type=str, required=True, help='TSV file containing genome data')
    return parser.parse_args()

def read_tsv(file):
    return pd.read_csv(file, sep='\t')

def compute_kmers(sequence, k):
    return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}

def jaccard_similarity(set1, set2):
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    return len(intersection) / len(union) if union else 0

def main():
    args = parse_arguments()
    data = read_tsv(args.data)

    genome_pairs = list(combinations(data.index, 2))
    dna_similarities = []
    protein_similarities = []
    nonzero_protein_similarities = 0
    nonzero_dna_similarities = 0
    num_pairs_processed = 0

    for idx1, idx2 in genome_pairs:
        genome1 = data.iloc[idx1]
        genome2 = data.iloc[idx2]

        prot_kmers1 = compute_kmers(genome1['seed_protein'], args.protK)
        prot_kmers2 = compute_kmers(genome2['seed_protein'], args.protK)
        dna_kmers1 = compute_kmers(genome1['seed_dna'], args.dnaK)
        dna_kmers2 = compute_kmers(genome2['seed_dna'], args.dnaK)

        prot_jaccard = jaccard_similarity(prot_kmers1, prot_kmers2)
        dna_jaccard = jaccard_similarity(dna_kmers1, dna_kmers2)

        dna_similarities.append(dna_jaccard)
        protein_similarities.append(prot_jaccard)

        if prot_jaccard > 0:
            nonzero_protein_similarities += 1
        if dna_jaccard > 0:
            nonzero_dna_similarities += 1

        num_pairs_processed += 1

    plt.scatter(dna_similarities, protein_similarities, s=1)
    plt.xlabel('DNA Jaccard Similarity')
    plt.ylabel('Protein Jaccard Similarity')
    plt.title('Jaccard Similarity Pairs')
    plt.show()

    print(f'Number of pairs processed: {num_pairs_processed}', file=sys.stderr)
    print(f'Nonzero protein similarities: {nonzero_protein_similarities}', file=sys.stderr)
    print(f'Nonzero DNA similarities: {nonzero_dna_similarities}', file=sys.stderr)

if __name__ == '__main__':
    main()
