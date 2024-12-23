######################################################################
#...Prompt that generated this program:
"""
Please write a Python program that will:
* read a tab-separated-value data-file with header-line from STDIN
* print the field-names in that TSV-file's header-line columns to STDOUT, and then exit
* any error-messages should be printed to STDERR.
"""
######################################################################
#...Pseudocode for this program:
"""
function main:
    try:
        # Read the TSV input from STDIN
        reader = initialize CSV reader using STDIN with tab ('\t') delimiter

        # Retrieve the first row (header line) from the TSV input
        header = get next row from reader

        # Print the header fields to STDOUT, joined by tabs
        print join header fields with tab ('\t') and output to STDOUT

    catch CSV error as e:
        # If a CSV processing error occurs, print the error message to STDERR
        print "Error processing CSV file: " + error message e to STDERR

    catch general error as e:
        # Handle any other exceptions by printing to STDERR
        print "An error occurred: " + error message e to STDERR

if __name__ == "main":
    call main
"""
######################################################################
#...Program generated by Grimoire:
import sys
import csv

def main():
    try:
        # Reading from STDIN
        reader = csv.reader(sys.stdin, delimiter='\t')
        
        # Get the header (field names)
        header = next(reader)
        
        # Print the field names (header) to STDOUT
        print("\t".join(header))
    
    except csv.Error as e:
        # Print error message to STDERR
        print(f"Error processing CSV file: {e}", file=sys.stderr)
    except Exception as e:
        # Handle other errors and print them to STDERR
        print(f"An error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
