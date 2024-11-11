################################################
#        Insert your Pseudo Code Below         #
# between the two rows of "triple-quotes",     #
# which will prevent the pseudocode from       #
# registering as code                          #                             #
"""
1. Initialize Variables:
  * Create variables to store
    - the total sequence count,
    - total sequence length,
    - a list to store sequence records.

2. Read Input (FASTA Format):
  * Read the FASTA file line-by-line.

  * Identify headers (lines starting with >),
    and treat the first word in the header as sequence_id,
    with the rest as description.

  * For subsequent lines, if they do not start with >,
    treat them as part of the sequence.

3. Extract Information:
  * For each sequence entry, calculate its length.
  * Store each sequence ID and its length in a list.

4. Output to STDOUT and STDERR:
  * Print sequence_id and sequence_length as tab-separated values.
  * Print the number of sequences and average sequence length to STDERR.

5. Exit Program:
  * Exit cleanly after printing the required details.
"""
################################################
import sys

def main():
    sequence_count = 0
    total_length = 0
    sequences = []

    current_id = None
    current_sequence = []

    # Reading from STDIN
    for line in sys.stdin:
        line = line.strip()
        
        if line.startswith('>'):  # Header line
            if current_id is not None:
                # Calculate length of current sequence and store result
                sequence_length = len("".join(current_sequence))
                sequences.append((current_id, sequence_length))
                sequence_count += 1
                total_length += sequence_length
            
            # Parse header
            header_parts = line[1:].split(maxsplit=1)
            current_id = header_parts[0]
            current_sequence = []  # Reset sequence storage

        else:  # Sequence line
            current_sequence.append(line)
    
    # Process the last sequence in file
    if current_id is not None:
        sequence_length = len("".join(current_sequence))
        sequences.append((current_id, sequence_length))
        sequence_count += 1
        total_length += sequence_length

    # Print results to STDOUT
    print("sequence_id\tsequence_length")
    for seq_id, seq_len in sequences:
        print(f"{seq_id}\t{seq_len}")

    # Calculate average length and print summary to STDERR
    if sequence_count > 0:
        average_length = total_length / sequence_count
    else:
        average_length = 0

    print(f"{sequence_count} sequences read.", file=sys.stderr)
    print(f"Average sequence length: {average_length:.2f}", file=sys.stderr)

if __name__ == "__main__":
    main()
