import csv
import sys


if len(sys.argv) != 2
    print("Usage: python tsv_headers.py <filename>")
    sys.exit(1)

filename = sys.argv[1]

try:
    with open(filename, newline='') as file
        reader = csv.reader(file, delimiter='\t')
        headers = next(reader)
        print("Field names in the TSV file are:")
        for header in headers
            print(header)
except Exception as e:
    print(f"An error occurred while reading the file: {e}")
    sys.exit(1)

