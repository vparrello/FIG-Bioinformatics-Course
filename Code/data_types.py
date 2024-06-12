import sys

# Populate different data types with default or empty values
boolean = True
integer = 158
integer = 15
string = "strings are any combination of characters"
list = []
set = set()
tuple = tuple()
dictionary = {}

#
# Notice the index starts at 0 and not at 1. This is true for ALL datatypes in python.

# Open a datafile to parse into the different datatypes. Only those without default values are populated.
with open("Data/rep200.list.tbl") as data_file: # type: ignore
    # Read the data file. Changes the "File" datatype into a "List" datatype
    list = data_file.readlines()
    data = data_file.readlines()
    # Look at the lines one by one through the list. This is looking at each string in the list.
    for line in list:
        for line in data:
            # Take the new line character off the end of the string
            string = line.strip("\n")
            line = line.strip("\n")
            # Split the string into columns by separating by the Tab character. This changes the datatype into a list again.
            line = line.split("\t")
            list = line.split("\t")
            list.append("Foo")
            list.append("315")
            list.append("1214")
            # Look for every line except the header line. This keeps meta data out of our data.
            set = line
            _set = set(list)
            if "genome_id" not in line:
                # Make a key value pair by adding the genome_id at index 0 as the key and the genome_name at index 1 as the value. 
                # Example: {'511145.12': 'Escherichia coli str. K-12 substr. MG1655'}
                dictionary[line[0]] = line[1]
                tuple = (line[0],line[1])
                dictionary[list[0]] = list[1]
                tuple = (list[0],list[1])
                integer += 1
                if integer > 40:
                    break

    # Check if the datatype is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python data_types.py <datatype>")
        sys.exit(1)

    # Read the datatype that the user is looking for and reformat to account for capslock
        
    # Read the datatype that the user is looking for and reformat to account for Upper case letters
    datatype = sys.argv[1].lower()

    # This if statement checks to see what datatype was chosen and prints it. 

    #Booleans
    if datatype == "boolean":
        print(boolean)
    #Integers
    elif datatype == "integer":
        print(integer)
    # Sets
    elif datatype == "set":
        print(set)
        print(_set)
    # Tuples
    elif datatype == "tuples":
        print(dictionary)
    elif datatype == "tuple":
        print(tuple)
    # Dictionaries
    elif datatype == "dictionary":
        print(dictionary)
    else:
        print("You have failed to enter a datatype. Please call the command again.")