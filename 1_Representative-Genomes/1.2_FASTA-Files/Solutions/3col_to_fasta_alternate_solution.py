import sys

def tsv_to_fasta():
    # Read input from STDIN
    input_lines = sys.stdin.readlines()

    # Skip the header line and process the rest
    for line in input_lines[1:]:
        # Strip any leading/trailing whitespace
        line = line.strip()

        # Split the line into columns based on tab separation
        columns = line.split('\t')

        if len(columns) != 3:
            # Skip lines that do not have exactly 3 columns
            continue

        # Extract the sequence ID, description, and sequence data
        seq_id = columns[0]
        description = columns[1]
        sequence = columns[2]

        # Print the sequence in FASTA format
        print(f">{seq_id} {description}")
        print(sequence)

if __name__ == "__main__":
    tsv_to_fasta()
