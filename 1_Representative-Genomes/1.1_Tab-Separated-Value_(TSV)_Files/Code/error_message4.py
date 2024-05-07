import csv
import sys


if len(sys.argv) != 2:
    print("Usage: python tsv_headers.py <filename>")
    sys.exit(1)

filename = sys.argv[1]
counter = 1
try:
    with open(filename, newline='') as file:
        # Create a reader object that will read the file as a TSV
        reader = csv.reader(file, delimiter='\t')
        # Extract the headers (first row of the TSV file)
        headers = next(reader)
        # Print the headers
        print("Field names in the TSV file are:")
        for header in headers:
            print(header)
            print(f'Next up is: {headers[counter]}')
            counter = counter + 1
except Exception as e:
    print(f"An error occurred while reading the file: {e}")
    sys.exit(1)
