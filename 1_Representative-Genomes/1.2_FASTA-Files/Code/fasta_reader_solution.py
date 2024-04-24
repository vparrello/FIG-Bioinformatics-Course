#           Insert your Psuedo Code Below      #
# Use ### before your text to ensure it does   #
# not register as code. The three quotes below #
# will protect you                             #
"""
To write a Python program that reads a named FASTA file
and performs the requested operations, we'll follow these steps
in the program:

1. Open and Read the FASTA File: Parse the FASTA file line by line.
  A FASTA file starts each new sequence with a line that begins
  with ">", followed by the sequence ID. The lines that follow,
  up until the next ">" or the end of the file, are the sequence
  itself.

2. Process Each Sequence: We'll store the ID and the sequence
  length, printing them to standard output (STDOUT) in a
  tab-separated format.

3. Calculate the Average Length: Once all sequences are processed,
  we'll calculate the average length of the sequences and print this
  to standard error (STDERR).

4. Error Handling: Add basic error handling for file operations. 
"""
################################################

import sys

def read_fasta_file(filename):
    """
    Reads a FASTA file and returns a list of tuples containing (id, sequence).
    The sequence ID is defined as the portion of the header line immediately
    following '>' and up to the first whitespace character. Any text following
    the first whitespace after the ID is considered a description and ignored
    for the purposes of this script.
    """
    with open(filename, 'r') as file:
        sequences = []
        seq_id = ''
        seq = ''
        
        for line in file:
            if line.startswith('>'):
                if seq_id:  # Save previous sequence before starting a new one
                    sequences.append((seq_id, seq))
                    seq = ''
                # Extract sequence ID up to the first whitespace, ignore description
                seq_id = line.strip().split()[0][1:]
            else:
                seq += line.strip()
        
        if seq_id:  # Save the last sequence
            sequences.append((seq_id, seq))
    
    return sequences

def main(filename):
    sequences = read_fasta_file(filename)
    total_length = 0
    
    for seq_id, seq in sequences:
        print(f"{seq_id}\t{len(seq)}")
        total_length += len(seq)
    
    average_length = total_length / len(sequences) if sequences else 0
    num_sequences = len(sequences)
    
    # Print to STDERR
    print(f"Number of sequences read: {num_sequences}", file=sys.stderr)
    print(f"Average length: {average_length}", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <FASTA filename>", file=sys.stderr)
        sys.exit(1)
    
    fasta_filename = sys.argv[1]
    try:
        main(fasta_filename)
    except FileNotFoundError:
        print(f"Error: File '{fasta_filename}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)
