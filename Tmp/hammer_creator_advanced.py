import sys
import re
from Bio import SeqIO
from collections import defaultdict

def main():
    # Parse command-line arguments
    if len(sys.argv) != 3 or sys.argv[1] != '-K':
        print("Usage: python hammer_creator_advanced.py -K <Kmer-length>", file=sys.stderr)
        sys.exit(1)

    try:
        kmer_length = int(sys.argv[2])
    except ValueError:
        print("Error: Kmer-length must be an integer.", file=sys.stderr)
        sys.exit(1)

    # Initialize data structures
    feature_to_sequence = {}
    feature_to_role = {}
    feature_to_genome = {}
    kmer_dict = defaultdict(lambda: defaultdict(int))

    sequence_count = 0
    kmer_count = 0

    # Process FASTA from STDIN
    for record in SeqIO.parse(sys.stdin, "fasta"):
        sequence_count += 1

        # Extract feature-ID and role
        header = record.description
        match = re.match(r"(\S+)\s+(.*)", header)
        if not match:
            print(f"Error: Unable to parse FASTA header: {header}", file=sys.stderr)
            continue

        feature_id = match.group(1)
        role = match.group(2).strip()

        # Extract genome-ID from feature-ID
        genome_match = re.match(r"fig\|(\d+\.\d+)\.peg\.\d+", feature_id)
        if not genome_match:
            print(f"Error: Invalid feature-ID format: {feature_id}", file=sys.stderr)
            continue

        genome_id = genome_match.group(1)
        
        # Store mappings
        feature_to_sequence[feature_id] = str(record.seq).lower()
        feature_to_role[feature_id] = role
        feature_to_genome[feature_id] = genome_id

    # Build the dictionary of dictionaries for Kmers
    for feature_id, sequence in feature_to_sequence.items():
        for i in range(len(sequence) - kmer_length + 1):
            kmer = sequence[i:i + kmer_length]
            kmer_dict[kmer][feature_id] += 1
            kmer_count += 1

    # Find and output hammers
    print("hammer\tfid\trole")
    hammer_count = 0
    
    for kmer, feature_counts in sorted(kmer_dict.items()):
        if len(feature_counts) == 1:
            feature_id, count = next(iter(feature_counts.items()))
            if count == 1:
                hammer_count += 1
                role = feature_to_role[feature_id]
                print(f"{kmer}\t{feature_id}\t{role}")

    # Print summary statistics to STDERR
    print(f"Sequences read: {sequence_count}", file=sys.stderr)
    print(f"Kmers processed: {kmer_count}", file=sys.stderr)
    print(f"Hammers found: {hammer_count}", file=sys.stderr)

if __name__ == "__main__":
    main()
