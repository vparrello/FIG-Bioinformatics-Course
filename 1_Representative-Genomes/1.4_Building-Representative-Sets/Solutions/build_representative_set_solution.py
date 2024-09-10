import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import Counter

def parse_arguments():
    parser = argparse.ArgumentParser(description='Compute a set of representative sequences using the Stingy Addition algorithm.')
    parser.add_argument('-k', '--kmer-length', type=int, required=True, help='Length of the K-mers')
    parser.add_argument('-s', '--similarity-threshold', type=int, required=True, help='Similarity threshold (number of K-mers in common)')
    parser.add_argument('-f', '--input-fasta', type=str, required=True, help='Input FASTA file')
    parser.add_argument('-r', '--output-repseq', type=str, required=True, help='Output representative sequences FASTA file')
    return parser.parse_args()

def extract_kmers(sequence, k):
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def calculate_similarity(seq1_kmers, seq2_kmers):
    return len(set(seq1_kmers) & set(seq2_kmers))

def main():
    args = parse_arguments()
    
    sequences = list(SeqIO.parse(args.input_fasta, "fasta"))
    
    kmer_length = args.kmer_length
    similarity_threshold = args.similarity_threshold
    
    RepGenSet = []
    
    for seq_record in sequences:
        sequence_kmers = extract_kmers(str(seq_record.seq), kmer_length)
        is_representative = True
        
        for rep_seq_record in RepGenSet:
            rep_sequence_kmers = extract_kmers(str(rep_seq_record.seq), kmer_length)
            similarity = calculate_similarity(sequence_kmers, rep_sequence_kmers)
            
            if similarity >= similarity_threshold:
                is_representative = False
                break
        
        if is_representative:
            RepGenSet.append(seq_record)
    
    with open(args.output_repseq, "w") as output_handle:
        SeqIO.write(RepGenSet, output_handle, "fasta")

if __name__ == "__main__":
    main()
