################################################
#   Prompt that generated this code:
"""
I have uploaded a program that I'd like you to modify.
Please make the following revisions to the uploaded program:

* The program should accept its input datafile-name
via the command-line argument '-i' (short for "Input"),
e.g. '-i datafilename'.

* The program should extract the header-names
from the first line of the input datafile.

* The program can accept an optional '-m' argument
defining a "skip-factor" or "step-factor" by which
it will skip over  unselected data-columns,
e.g., '-m 10' means "print every 10th column",
so that columns [1, 11, 21, 31...] will be selected.

* The program can accept an optional argument '-n'
to specify the "total number of selected data-columns
that will be printed to STDOUT", e.g. -n 4 means
"Print a TOTAL of 4 selected columns to STDOUT".
Please note that this argument specifies the TOTAL number of columns
to be printed to STDOUT --- it does _NOT_ specify the maximum 
column-number to be printed!
If this optional argument is not specified,
then the program should print all of the selected columns
in the input datafile.
"""
################################################
#   Pseudocode for this program:
"""
function main():
    # Set up argument parser
    parser = create_argument_parser(description="Process a TSV file and print selected columns")
    
    # Define required and optional arguments
    parser.add_argument('-i', '--input', required=True, help="Input data file name")
    parser.add_argument('-m', '--skip_factor', type=int, default=1, help="Step factor to select columns")
    parser.add_argument('-n', '--num_columns', type=int, help="Total number of selected columns to print")
    
    # Parse command-line arguments
    args = parser.parse_args()

    try:
        # Open the input file in read mode
        with open(args.input, 'r') as tsvfile:
            # Initialize CSV reader with tab ('\\t') delimiter
            reader = csv.reader(tsvfile, delimiter='\\t')

            # Extract the header (first row) from the TSV file
            header = next(reader)

            # Select every mth column from the header based on the skip factor
            step = args.skip_factor
            selected_columns = header[::step]

            # If the number of columns to print is specified, truncate the selection
            if args.num_columns:
                selected_columns = selected_columns[:args.num_columns]

            # Print the selected columns to STDOUT, joined by tabs ('\\t')
            print(join(selected_columns, '\\t'))

    except Exception as e:
        # In case of any error, print the error message to STDERR
        print("Error: " + str(e), to=STDERR)

# Execute the main function
if __name__ == "__main__":
    main()
"""
################################################
#   Grimoire-generated code:
import csv
import argparse
import sys

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process a TSV file and print selected columns.")
    parser.add_argument('-i', '--input', required=True, help="Input data file name")
    parser.add_argument('-m', '--skip_factor', type=int, default=1, help="Step factor to select columns")
    parser.add_argument('-n', '--num_columns', type=int, help="Total number of selected columns to print")
    
    # Parse the arguments
    args = parser.parse_args()
    
    try:
        # Open the input file
        with open(args.input, 'r') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            
            # Extract the header (first row)
            header = next(reader)
            
            # Select every mth column from the header
            step = args.skip_factor
            selected_columns = header[::step]
            
            # If the number of columns to print is specified, truncate the selection
            if args.num_columns:
                selected_columns = selected_columns[:args.num_columns]
            
            # Print the selected columns to STDOUT
            print('\t'.join(selected_columns))
    
    except Exception as e:
        # If any error occurs, print to STDERR
        print(f"Error: {str(e)}", file=sys.stderr)

if __name__ == "__main__":
    main()
