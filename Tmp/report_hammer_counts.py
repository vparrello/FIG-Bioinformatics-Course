import sys
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

def read_hammer_file(file_path):
    """
    Reads a tab-separated value file, skips the header, and loads the first and second columns
    into a dictionary as keys (hammers) and values (feature_ids). Also returns the hammer length.
    """
    hammer_dict = {}
    hammer_length = None
    with open(file_path, 'r') as f:
        # Skip the header line
        next(f)
        for line in f:
            columns = line.strip().split('\t')
            hammer, feature_id = columns[0], columns[1]
            hammer_dict[hammer] = feature_id
            if hammer_length is None:
                hammer_length = len(hammer)
    return hammer_length, hammer_dict

def find_kmers(sequence, k):
    """
    Finds all Kmers of length k in the given DNA sequence.
    """
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]

def process_fasta_and_count_kmers(hammer_length, hammer_dict):
    """
    Processes a FASTA-formatted DNA file from STDIN, finds all Kmers of the hammer length
    in both the forward and reverse complement sequences, and counts occurrences of each Kmer
    that is a hammer.
    """
    kmer_counts = defaultdict(int)

    # Reading FASTA from stdin
    for record in SeqIO.parse(sys.stdin, "fasta"):
        sequence = str(record.seq)
        reverse_complement = str(record.seq.reverse_complement())

        # Find kmers in both forward and reverse complement
        for kmer in find_kmers(sequence, hammer_length) + find_kmers(reverse_complement, hammer_length):
            if kmer in hammer_dict:
                kmer_counts[kmer] += 1

    return kmer_counts

def generate_report(kmer_counts, hammer_dict):
    """
    Generates a tab-separated report of hammers found more than once in the DNA sequences,
    sorted in decreasing order of count. Also includes the associated feature_id.
    """
    report_lines = []
    for hammer, count in sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True):
        if count > 0:
            feature_id = hammer_dict[hammer]
            report_lines.append(f"{hammer}\t{count}\t{feature_id}")
    
    return report_lines

if __name__ == "__main__":
    # Example hammer file path
    hammer_file = sys.argv[1]

    # Step 1: Read the hammer file
    hammer_length, hammer_dict = read_hammer_file(hammer_file)

    # Step 2: Process the FASTA file and count Kmers
    kmer_counts = process_fasta_and_count_kmers(hammer_length, hammer_dict)

    # Step 3: Generate and print the report
    report = generate_report(kmer_counts, hammer_dict)
    print("Hammer\tCount\tFeature_ID")
    print("\n".join(report))
