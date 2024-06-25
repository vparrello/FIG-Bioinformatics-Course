import sys

def main():
    try:
        # Read from standard input
        input_data = sys.stdin.read()
        
        # Split the input data into lines
        lines = input_data.splitlines()
        
        if not lines:
            print("No data found.", file=sys.stderr)
            return
        
        # The first line is the header
        header = lines[0]
        
        # Split the header line by tab to get the field names
        field_names = header.split('\t')
        
        # Print the field names to standard output
        for field_name in field_names:
            print(field_name)
    
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
