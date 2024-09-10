import sys
import argparse

# Use Argparse module to grab the information needed from the command line based on the arguments given.
parser = argparse.ArgumentParser()
# This is the input arguement. This is the file that will be parsed for its headers
parser.add_argument("-i", "--input", help="The input file that gets read from the program. It must be in a Tab Seperated format (TSV)")
# This is the number of columns that will get printed into standard out
parser.add_argument("-n", "--number", help="Number of data columns to be printed into the standard out output. If not specifies, the program should print all of the columns", default=200)
# This is the skip counter for which multiple of the columns you want to see (Examples: Every 4th column, every other column, and so on)
parser.add_argument("-m", "--multiple", help="The number of columns that can get skipped. For example: if 4 is specified, it will print every fourth column.", default=1)

# This creates variables for each of the parsed inputs above and creates values in those inputs. The solution example has these values inserted into the variables. Use these inputs to help you visualized what the code does
#   args.input = 'Data/data.tbl'
#   args.number = 8
#   args.multiple = 4

args = parser.parse_args()
try:
    input_file = args.input
except:
    print("No file given. Be sure to specify an input file using -i or --input")
    quit()

def main():

    try:
        counter = 0
        # Open the input file
        with open(input_file) as input:
            
            # Read from standard input and split into lines. The data is now a list
            input = input.readlines()
            
            # Seperate the first line out from the rest of the file to only grab the header line. It is now a string.
            header = input[0]
            
            # Split the header line into its columns. It is a list again.
            field_names = header.split('\t')
            
            # Iterate through the field names list
            for i in range(len(field_names)):
                # Check to make sure we are still under the column count. This is a seperate counter because not ever column in i will be printed.
                if counter <= int(args.number):
                    # Skip every column that does not count for the skip counter
                    if i % int(args.multiple) == 0:
                        counter += 1
                        print(field_names[i])

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
