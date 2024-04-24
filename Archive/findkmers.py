'''
@author Victoria Parrello
Week 1 - Representative Genomes
Activity 1: Find the secret sequences
'''
dna_sequences = "acgt"
header = ""
counter = 0

with open("1 - Representative Genomes/Week 1/Data/Sample1.fasta") as datafile:
    data = #TODO put something about *readlines* in here
    for line in data:
        if "TODO insert header pattern here" in line:
            #TODO if the line is a header and not a sequence, catch the value in the header variable.
        #TODO Find the sequence in the line and catch it
        #TODO can you *strip* the sequence from all of the dna pieces?
        #TODO print your header and sequences. Bonus points if you can number them