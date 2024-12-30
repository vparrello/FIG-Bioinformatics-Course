# Solutions for TSV Exercise 5
This file contains the error messages thrown for each file and the fixed versions of the code for reference. The command used to call each file will be the first line. Then the error messages will be commented at the beginning of each code block and the fixed code will be directly following.

## error_message1.py

### Command:
```
python Code/error_message1.py Data/bindict.tbl
```

### Error Message:
```
  File "C:\Users\Directory\FIG-Bioinformatics-Course\Code\error_message1.py", line 5
    if len(sys.argv) != 2
                         ^
    SyntaxError: expected ':'

  File "C:\Users\Directory\FIG-Bioinformatics-Course\Code\error_message1.py", line 12
    with open(filename, newline='') as file
                                           ^
    SyntaxError: expected ':'

  File "C:\Users\Directory\FIG-Bioinformatics-Course\Code\error_message1.py", line 16
    for header in headers
                         ^
    SyntaxError: expected ':'
``` 

### Fixed Code:
```
import csv
import sys


if len(sys.argv) != 2: # Error thrown here
    print("Usage: python tsv_headers.py <filename>")
    sys.exit(1)

filename = sys.argv[1]

try:
    with open(filename, newline='') as file: # Error thrown here
        reader = csv.reader(file, delimiter='\t')
        headers = next(reader)
        print("Field names in the TSV file are:")
        for header in headers: #Error thrown here
            print(header)
except Exception as e:
    print(f"An error occurred while reading the file: {e}")
    sys.exit(1)
```

## error_message2.py
### Command:
```
python Code/error_message2.py Data/bindict.tbl
```
### Error Message:
```
An error occurred while reading the file: name 'filename' is not defined
```

### Fixed Code:
```
import csv
import sys

#Check if the filename is provided as a command-line argument
if len(sys.argv) != 2:
    print("Usage: python tsv_headers.py <filename>")
    sys.exit(1)

# VSCode doesn't know something is missing here. Look below for a yellow squiggly line for a clue as to what variable needs to be started here.
filename = sys.argv[1]

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
except Exception as e:
    print(f"An error occurred while reading the file: {e}")
    sys.exit(1)
```

## error_message3.py
### Command:
```
python Code/error_message3.py 
```
### Error Message:
```
Traceback (most recent call last):
  File "C:\Users\Directory\FIG-Bioinformatics-Course\Code\error_message3.py", line 12, in <module>
    answer = add(x,y)
             ^^^^^^^^
  File "C:\Users\Directory\FIG-Bioinformatics-Course\Code\error_message3.py", line 5, in add
    total = a + b
            ~~^~~
TypeError: unsupported operand type(s) for +: 'int' and 'str'
```

### Fixed Code:
```
import sys

#This is expecting two numbers
def add(a,b):
    total = a + b
    return total


x = 10
y = 5
# Remove the string version of y
answer = add(x,y)
print(answer)
```

## error_message4.py
### Command:
```
python Code/error_message4.py Data/bindict.tbl
```
### Error Message:
```
Field names in the TSV file are:
genome_id
Next up is: genome_name
genome_name
Next up is: RepGen.200
RepGen.200
Next up is: RepGen.100
RepGen.100
Next up is: RepGen.50
RepGen.50
An error occurred while reading the file: list index out of range
```

### Fixed Code:
```
import csv
import sys

if len(sys.argv) != 2:
    print("Usage: python tsv_headers.py <filename>")
    sys.exit(1)

filename = sys.argv[1]
# Remove the counter variable
# counter = 1
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
            # Remove the next 2 lines as they call the list index out of range error
            # print(f'Next up is: {headers[counter]}')
            # counter = counter + 1
except Exception as e:
    print(f"An error occurred while reading the file: {e}")
    sys.exit(1)
```

## error_message5.py
### Command:
```
python Code/error_message5.py Data/bindict.tbl
```
### Error Message:
```
Field names in the TSV file are:
An error occurred while reading the file: 'str' object has no attribute 'reverse'
```

### Fixed Code:
``` 
import csv
import sys


if len(sys.argv) != 2:
    print("Usage: python tsv_headers.py <filename>")
    sys.exit(1)

filename = sys.argv[1]

try:
    with open(filename, newline='') as file:
        # Create a reader object that will read the file as a TSV
        reader = csv.reader(file, delimiter='\t')
        # Extract the headers (first row of the TSV file)
        headers = next(reader)
        # Print the headers
        print("Field names in the TSV file are:")
        for header in headers:
            # Remove this line since we don't want to reverse the headers. Also reverse doesn't exist according to the error message.
            # header.reverse()
            print(header)

except Exception as e:
    print(f"An error occurred while reading the file: {e}")
    sys.exit(1)

```

## error_message_final.py
### Command:
```
python Code/error_message_final.py Data/bindict.tbl
```
### Error Messages:
```
  File "C:\Users\Directory\FIG-Bioinformatics-Course\Code\error_message_final.py", line 4
    if len(sys.argv) != 2
                         ^
    SyntaxError: expected ':'


  File "C:\Users\Directory\FIG-Bioinformatics-Course\Code\error_message_final.py", line 20
    print(f"An error occurred while reading the file: {e})
          ^
    SyntaxError: unterminated f-string literal (detected at line 20)


    An error occurred while reading the file: name 'reader' is not defined


    An error occurred while reading the file: '_csv.reader' object has no attribute 'reverse'
```

### Fixed Code:
```
import csv
import sys

if len(sys.argv) != 2: # Error thrown here for a missing colon
    print("Usage: python tsv_headers.py <filename>")
    sys.exit(1)

filename = sys.argv[1]

try:
    with open(filename, newline='') as file:

        
        reader = csv.reader(file, delimiter='\t')
        # This line was throwing 2 errors. First it had to come after the reader object above. Second it was trying to reverse the headers which isn't a function that reader has.
        headers = next(reader)
        print("Field names in the TSV file are:")
        for header in headers:
            print(header)
except Exception as e:
    # This line was throwing an error because it was missing its ending quote between the } and the ) in the line.
    print(f"An error occurred while reading the file: {e}")
    sys.exit(1)

```
