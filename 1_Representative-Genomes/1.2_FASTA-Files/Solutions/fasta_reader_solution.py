################################################
#        Insert your Pseudo Code Below         #
# between the two rows of "triple-quotes",     #
# which will prevent the pseudocode from       #
# registering as code                          #                             #
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

def process_fasta_file():
    sequence_count = 0
    total_length = 0
    current_id = None
    current_sequence = []

    for line in sys.stdin:
        line = line.strip()
        if line.startswith('>'):
            if current_id is not None:
                sequence_length = len(''.join(current_sequence))
                print(f"{current_id}\t{sequence_length}")
                sequence_count += 1
                total_length += sequence_length
            
            current_id = line[1:].split()[0]
            current_sequence = []
        else:
            current_sequence.append(line)
    
    if current_id is not None:
        sequence_length = len(''.join(current_sequence))
        print(f"{current_id}\t{sequence_length}")
        sequence_count += 1
        total_length += sequence_length
    
    average_length = total_length / sequence_count if sequence_count > 0 else 0
    print(f"Number of sequences: {sequence_count}", file=sys.stderr)
    print(f"Average sequence length: {average_length:.2f}", file=sys.stderr)

if __name__ == "__main__":
    process_fasta_file()
