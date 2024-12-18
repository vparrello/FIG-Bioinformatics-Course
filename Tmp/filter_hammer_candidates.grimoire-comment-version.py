"""
Python Program: filter_hammer_candidates.py

Requirements:

1. Accept Arguments:
   - Accept a mandatory argument `-C` or `--contigs-directory` specifying the path to a directory containing FASTA files.

2. Read Input:
   - Read a tab-separated file from standard input.
     - Save the first line (header).
     - Store the rest in a dictionary (`hammer_dict`) with:
       - Keys: The first field.
       - Values: The entire line.
     - Return:
       - Number of entries read.
       - Hammer length.
       - The header.
       - The dictionary.

3. Process Contigs Directory:
   - For each file in the directory:
     - Print `Processing '<filename>'` to standard error.
     - Use BioPython to parse the FASTA sequences in the file.

4. Find Kmers:
   - For each contig:
     - Extract all kmers of hammer length from:
       - The sequence.
       - Its reverse complement.

5. Match Hammers:
   - For each kmer:
     - Check if it exists in `hammer_dict`.
     - If found:
       - Increment its count in a new dictionary.
       - If the count exceeds `1`, remove it from `hammer_dict`.

6. Output Results:
   - Print:
     - The saved header.
     - All remaining entries in `hammer_dict`, sorted by their keys.

7. Report Statistics:
   - Print the number of:
     - Candidates read.
     - Candidates eliminated.
     - Hammers accepted.
"""
