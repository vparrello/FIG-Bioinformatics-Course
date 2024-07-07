import sys

def interactive_exercise(data_type=None):
    valid_data_types = ["boolean", "integer", "string", "list", "set", "tuple", "dictionary"]
    python_data_types = ["int", "float", "complex", "bool", "str", "list", "tuple", "range", "dict", "set", "frozenset", "bytes", "bytearray", "memoryview"]

    if data_type is None or data_type.lower() == "types":
        print("Datatypes described in this program:")
        for dt in valid_data_types:
            print(f"- {dt.capitalize()}")
        return

    data_type = data_type.lower()

    if data_type not in valid_data_types:
        if data_type in python_data_types:
            print("This datatype is valid, but will not be discussed within this exercise.")
        else:
            print("This datatype is not a valid python datatype")
        return

    def explain_and_example(data_type, description, example):
        print(f"\n{data_type.capitalize()}")
        print(description)
        input("Press Enter to see an example...")
        print(example)

    if data_type == "boolean":
        explain_and_example(
            "Boolean",
            "Booleans represent two values: True or False. Useful for conditional checks.",
            """
# Example: Checking if a DNA sequence is valid
def is_valid_dna(sequence):
    valid_nucleotides = {'A', 'T', 'C', 'G'}
    return all(nucleotide in valid_nucleotides for nucleotide in sequence)

dna_sequence = "ATCGTTAGC"
is_valid = is_valid_dna(dna_sequence)
print(f"Is the DNA sequence valid? {is_valid}")  # Output: True
""")

    elif data_type == "integer":
        explain_and_example(
            "Integer",
            "Integers represent whole numbers. Useful for counting or indexing.",
            """
# Example: Counting the number of nucleotides in a DNA sequence
dna_sequence = "ATCGTTAGC"
nucleotide_count = len(dna_sequence)
print(f"The number of nucleotides: {nucleotide_count}")  # Output: 9
""")

    elif data_type == "string":
        explain_and_example(
            "String",
            "Strings represent sequences of characters. Useful for storing textual data like DNA sequences.",
            """
# Example: Reverse complement of a DNA sequence
def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[nucleotide] for nucleotide in reversed(sequence))

dna_sequence = "ATCG"
reverse_complement_sequence = reverse_complement(dna_sequence)
print(f"Reverse complement: {reverse_complement_sequence}")  # Output: CGAT
""")

    elif data_type == "list":
        explain_and_example(
            "List",
            "Lists represent ordered collections of items. Useful for storing sequences of data.",
            """
# Example: Storing multiple DNA sequences
dna_sequences = ["ATCG", "CGTA", "TTAGC"]
for sequence in dna_sequences:
    print(f"DNA sequence: {sequence}")
""")

    elif data_type == "set":
        explain_and_example(
            "Set",
            "Sets represent unordered collections of unique items. Useful for ensuring no duplicates.",
            """
# Example: Finding unique nucleotides in a DNA sequence
dna_sequence = "ATCGATTG"
unique_nucleotides = set(dna_sequence)
print(f"Unique nucleotides: {unique_nucleotides}")  # Output: {'A', 'C', 'G', 'T'}
""")

    elif data_type == "tuple":
        explain_and_example(
            "Tuple",
            "Tuples represent ordered collections of items that are immutable. Useful for fixed collections of data.",
            """
# Example: Storing a DNA sequence and its length
dna_info = ("ATCG", 4)
sequence, length = dna_info
print(f"DNA sequence: {sequence}, Length: {length}")  # Output: DNA sequence: ATCG, Length: 4
""")

    elif data_type == "dictionary":
        explain_and_example(
            "Dictionary",
            "Dictionaries represent collections of key-value pairs. Useful for mapping relationships.",
            """
# Example: Mapping nucleotide counts in a DNA sequence
def count_nucleotides(sequence):
    counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for nucleotide in sequence:
        counts[nucleotide] += 1
    return counts

dna_sequence = "ATCGATTG"
nucleotide_counts = count_nucleotides(dna_sequence)
print(f"Nucleotide counts: {nucleotide_counts}")  # Output: {'A': 3, 'T': 3, 'C': 1, 'G': 1}
""")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        data_type = sys.argv[1]
        interactive_exercise(data_type)
    else:
        interactive_exercise()
