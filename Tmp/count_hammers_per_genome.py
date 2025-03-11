import gzip
import re
import sys
import argparse

def parse_genome_file(genome_file):
    genome_dict = {}
    genome_pattern = re.compile(r"\d+\.\d+")

    with open(genome_file, 'r') as f:
        next(f)  # Skip header line
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            genome_id, genome_name = parts[0], parts[1]
            if genome_pattern.match(genome_id):
                genome_dict[genome_id] = genome_name
    return genome_dict

def process_hammer_file(hammer_file, genome_dict):
    hammer_counts = {}
    fid_pattern = re.compile(r"fig\|(\d+\.\d+)\.peg.\d+")

    open_func = gzip.open if hammer_file.endswith('.gz') else open
    with open_func(hammer_file, 'rt') as f:
        next(f)  # Skip header line
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            fid = parts[1]
            match = fid_pattern.match(fid)
            if match:
                genome_id = match.group(1)
                hammer_counts[genome_id] = hammer_counts.get(genome_id, 0) + 1

    return hammer_counts

def print_results(hammer_counts, genome_dict):
    print("hammers\tgenome_ID\tgenome_name")
    sorted_genomes = sorted(hammer_counts.keys(), key=lambda gid: genome_dict.get(gid, "Unknown"))
    for genome_id in sorted_genomes:
        genome_name = genome_dict.get(genome_id, "Unknown")
        hammer_count = hammer_counts[genome_id]
        print(f"{hammer_count}\t{genome_id}\t{genome_name}")

def main():
    parser = argparse.ArgumentParser(description="Process hammer and genome files.")
    parser.add_argument("-H", "--hammer-file", required=True, help="Hammer file (can be gzipped).")
    parser.add_argument("-G", "--genome-file", required=True, help="Genome file.")

    args = parser.parse_args()
    genome_dict = parse_genome_file(args.genome_file)
    hammer_counts = process_hammer_file(args.hammer_file, genome_dict)
    print_results(hammer_counts, genome_dict)

if __name__ == "__main__":
    main()
