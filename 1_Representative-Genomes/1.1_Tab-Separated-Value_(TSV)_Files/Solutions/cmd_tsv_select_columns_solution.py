import sys
import csv

def main():
    # Check if there are command-line arguments for column headers
    if len(sys.argv) < 2:
        sys.stderr.write("No column names provided. Please specify the columns you want to extract.\n")
        sys.exit(1)

    # Extract column names from command-line arguments
    column_names = sys.argv[1:]

    try:
        # Setup a CSV reader to read from standard input assuming tab-delimited format
        reader = csv.reader(sys.stdin, delimiter='\t')
        
        # Read the header line
        headers = next(reader, None)
        if headers is None:
            sys.stderr.write("The input appears to be empty.\n")
            sys.exit(1)

        # Find indices of the requested columns
        column_indices = []
        for name in column_names:
            if name in headers:
                column_indices.append(headers.index(name))
            else:
                sys.stderr.write(f"Column '{name}' not found in the header.\n")
                sys.exit(1)
        
        # Setup a CSV writer to write to standard output in tab-delimited format
        writer = csv.writer(sys.stdout, delimiter='\t')

        # Write the selected columns' headers
        writer.writerow([headers[i] for i in column_indices])
        
        # Process the remaining rows
        for row in reader:
            writer.writerow([row[i] for i in column_indices])
        
        # Inform the user of completion
        sys.stderr.write("Columns extracted and written successfully.\n")

    except Exception as e:
        sys.stderr.write(f"An error occurred: {str(e)}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
