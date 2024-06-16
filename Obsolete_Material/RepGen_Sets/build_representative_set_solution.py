import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_arguments():
    parser = argparse.ArgumentParser(description="Compute a set of representative sequences (RepGen set) using the Stingy Addition algorithm.")
    parser.add_argument('-k', '--kmer-length', type=int, required=True, help="Length of the Kmers")
    parser.add_argument('-m', '--minsim', type=float, required=True, help="Minimum similarity threshold (percent)")
    parser.add_argument('-f', '--input-fasta', type=str, required=True, help="Input FASTA file")
    parser.add_argument('-r', '--output-repseq', type=str, required=True, help="Output Representative Sequences file")
    return parser.parse_args()

def convert_min_sim(minsim_percent):
    return minsim_percent / 100.0

def read_fasta_file(input_fasta):
    return list(SeqIO.parse(input_fasta, "fasta"))

def sort_sequences(sequences):
    return sorted(sequences, key=lambda x: (-len(x.seq), str(x.seq), x.id))

def compute_jaccard_similarity(seq1, seq2, k):
    kmers1 = {str(seq1.seq[i:i+k]) for i in range(len(seq1.seq) - k + 1)}
    kmers2 = {str(seq2.seq[i:i+k]) for i in range(len(seq2.seq) - k + 1)}
    intersection = kmers1 & kmers2
    union = kmers1 | kmers2
    return len(intersection) / len(union)

def stingy_addition(sequences, k, minsim):
    repgen_set = []
    for seq1 in sequences:
        if all(compute_jaccard_similarity(seq1, seq2, k) < minsim for seq2 in repgen_set):
            repgen_set.append(seq1)
    return repgen_set

def write_fasta_file(sequences, output_fasta):
    SeqIO.write(sequences, output_fasta, "fasta")

def main():
    args = parse_arguments()
    k = args.kmer_length
    minsim = convert_min_sim(args.minsim)
    input_fasta = args.input_fasta
    output_repseq = args.output_repseq
    
    sequences = read_fasta_file(input_fasta)
    sorted_sequences = sort_sequences(sequences)
    repgen_set = stingy_addition(sorted_sequences, k, minsim)
    write_fasta_file(repgen_set, output_repseq)

if __name__ == "__main__":
    main()
