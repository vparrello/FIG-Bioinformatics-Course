import argparse
import csv
import re
import sys
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqIO import parse

def parse_hammers_file(hammers_file):
    hammer_to_feature = {}
    feature_to_role = {}
    feature_to_genome = {}

    with open(hammers_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        headers = reader.fieldnames

        for row in reader:
            hammer = row['hammer']
            feature_id = row['fid']
            role = row.get('role', 'Unknown')

            hammer_to_feature[hammer] = feature_id
            feature_to_role[feature_id] = role

            match = re.match(r'fig\|(\d+\.\d+)\.peg\.\d+', feature_id)
            if match:
                genome_id = match.group(1)
                feature_to_genome[feature_id] = genome_id

    kmer_length = len(next(iter(hammer_to_feature))) if hammer_to_feature else 0

    return kmer_length, hammer_to_feature, feature_to_role, feature_to_genome

def parse_genome_names_file(genome_names_file):
    genome_id_to_name = {}

    with open(genome_names_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader)  # Skip header

        for row in reader:
            genome_id, genome_name = row[0], row[1]
            genome_id_to_name[genome_id] = genome_name

    return genome_id_to_name

def process_sequences(kmer_length, hammer_to_feature, feature_to_genome, feature_to_role):
    genome_counts = defaultdict(int)
    genome_role_counts = defaultdict(lambda: defaultdict(int))

    for record in parse(sys.stdin, 'fasta'):
        sequence = str(record.seq).lower()
        reverse_complement = str(Seq(sequence).reverse_complement())

        for strand in (sequence, reverse_complement):
            for i in range(len(strand) - kmer_length + 1):
                kmer = strand[i:i + kmer_length]

                if kmer in hammer_to_feature:
                    feature_id = hammer_to_feature[kmer]
                    genome_id = feature_to_genome[feature_id]
                    role = feature_to_role[feature_id]

                    genome_counts[genome_id] += 1
                    genome_role_counts[genome_id][role] += 1

    return genome_counts, genome_role_counts

def main():
    parser = argparse.ArgumentParser(description="Advanced hammer comparison tool.")
    parser.add_argument('-H', '--hammers', required=True, help="Path to the hammers TSV file.")
    parser.add_argument('-G', '--genome-names', required=True, help="Path to the genome names TSV file.")
    parser.add_argument('-F', '--role-fraction', type=float, default=0.8, help="Minimum fraction of roles required.")

    args = parser.parse_args()

    kmer_length, hammer_to_feature, feature_to_role, feature_to_genome = parse_hammers_file(args.hammers)
    genome_id_to_name = parse_genome_names_file(args.genome_names)

    genome_counts, genome_role_counts = process_sequences(kmer_length, hammer_to_feature, feature_to_genome, feature_to_role)

    total_roles = len(set(feature_to_role.values()))

    output = []
    for genome_id, role_counts in genome_role_counts.items():
        roles_found = len(role_counts)
        score = genome_counts[genome_id]

        if roles_found >= total_roles * args.role_fraction:
            genome_name = genome_id_to_name.get(genome_id, 'Unknown sp.')
            if genome_name == 'Unknown sp.':
                print(f"Warning: Genome ID {genome_id} not found in genome names file.", file=sys.stderr)

            output.append((genome_id, genome_name, score, roles_found))

    output.sort(key=lambda x: x[2], reverse=True)

    print("Genome_ID\tGenome_Name\tScore\tRoles_Found")
    for genome_id, genome_name, score, roles_found in output:
        print(f"{genome_id}\t{genome_name}\t{score}\t{roles_found}")

if __name__ == '__main__':
    main()
