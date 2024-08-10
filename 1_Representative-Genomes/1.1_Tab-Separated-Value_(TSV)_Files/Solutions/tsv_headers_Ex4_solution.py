################################################
#   Prompt that generated this code:
"""
Please make the following revisions to the upload program:
* The program should accept its input datafile-name
via the command-line argument '-i' (short for "Input"),
e.g. '-i datafilename'.

* The program should extract the header-names
from the first line of the input datafile.

* The program should accept an optional argument '-n'
to specify the "total number of selected data-columns"
that will be printed to STDOUT, e.g. -n 4 means
"Print a total of 4 selected columns".
If this optional argument is not specified,
then the program should print all of the selected columns
in the input datafile.

* The program can take an optional '-m' argument
defining a "skip-factor" or "step-factor" by which
it will skip over the unselected data-columns,
e.g., '-m 10' means "print every 10th column",
so that columns [1, 11, 21, 31...] will be selected,
up to -n total columns.
"""
################################################
#   Pseudocode for this program:
"""
Function main:
    Create a parser object from ArgumentParser with the description "Extract headers from a TSV file and print selected columns."
    Add an argument '-i' or '--input' to the parser, which is required and represents the input data file name.
    Add an optional argument '-n' or '--num_columns' to the parser to specify the total number of data columns to print, with a default of printing all columns.
    Add an optional argument '-m' or '--skip_factor' to the parser to define how many columns to skip, with a default value of 1.

    Parse the arguments using the parser and store them in the variable args.

    Try to do the following:
        Open the file specified by args.input and read all its lines.
        
        If there are no lines in the file, print "No data found." to standard error and return from the function.

        Extract the first line from the lines as the header.
        Split the header line by tabs to get the field names.

        Select columns from the field names using the skip factor specified by args.skip_factor.
        If args.num_columns is specified, limit the selected columns to that number.

        For each selected field name, print it to standard output.

    If an exception occurs, catch it and print an error message to standard error.

If this script is being run as the main program, call the main function.
"""
################################################
#   Grimoire-generated code: 
import sys
import argparse

def main():
    try:
        # Parse command-line arguments
        parser = argparse.ArgumentParser(description="Extract header names from a TSV file")
        parser.add_argument("-i", "--input", required=True, help="Input datafile name")
        parser.add_argument("-n", "--number", type=int, help="Total number of selected data columns")
        parser.add_argument("-m", "--skip", type=int, default=1, help="Skip-factor or step-factor for selecting columns")
        
        args = parser.parse_args()
        
        # Read from the input file
        with open(args.input, 'r') as file:
            tsv_input = file.read()
        
        # Split the input into lines
        lines = tsv_input.strip().split('\n')
        
        # Check if there are any lines in the input
        if not lines:
            print("No input provided", file=sys.stderr)
            sys.exit(1)
        
        # Extract the header line
        header_line = lines[0]
        
        # Split the header line by tab characters to get field names
        field_names = header_line.split('\t')
        
        # Determine the columns to print
        selected_columns = field_names[::args.skip]
        
        # Limit the number of columns if the -n argument is specified
        if args.number:
            selected_columns = selected_columns[:args.number]
        
        # Print selected field names to STDOUT
        for field_name in selected_columns:
            print(field_name)
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
