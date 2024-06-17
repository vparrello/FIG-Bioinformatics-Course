import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def main():
    first_line = True
    records = []  # List to store SeqRecord objects
    
    for line in sys.stdin:
        if first_line:
            first_line = False
            continue
        
        parts = line.strip().split('\t')
        if len(parts) < 3:
            continue  # Skip lines that don't have enough columns
        
        seq_id = parts[0]
        description = parts[1]
        sequence_data = parts[2]

        # Create a SeqRecord object
        seq_record = SeqRecord(Seq(sequence_data),
                               id=seq_id,
                               description=description)
        records.append(seq_record)

    # Write all records to standard output in FASTA format
    SeqIO.write(records, sys.stdout, "fasta")

if __name__ == "__main__":
    main()
