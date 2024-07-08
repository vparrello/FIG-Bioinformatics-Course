import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_args():
    parser = argparse.ArgumentParser(description="Compute a Representative Set using the Stingy Addition algorithm.")
    parser.add_argument("-k", "--kmer_length", type=int, required=True, help="Kmer length")
    parser.add_argument("-s", "--similarity_threshold", type=int, required=True, help="Similarity threshold (number of Kmers in common)")
    parser.add_argument("-f", "--input_fasta", type=str, required=True, help="Input FASTA file")
    parser.add_argument("-r", "--output_repseq", type=str, required=True, help="Output RepSeq file")
    return parser.parse_args()

def get_kmers(sequence, k):
    return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}

def is_similar(seq1_kmers, seq2_kmers, threshold):
    common_kmers = seq1_kmers.intersection(seq2_kmers)
    return len(common_kmers) >= threshold

def main():
    args = parse_args()
    k = args.kmer_length
    sim_threshold = args.similarity_threshold
    input_fasta = args.input_fasta
    output_repseq = args.output_repseq
    
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    
    repgen_set = []
    for seq_record in sequences:
        seq_kmers = get_kmers(str(seq_record.seq), k)
        if all(not is_similar(seq_kmers, get_kmers(str(rep.seq), k), sim_threshold) for rep in repgen_set):
            repgen_set.append(seq_record)
    
    SeqIO.write(repgen_set, output_repseq, "fasta")

if __name__ == "__main__":
    main()
