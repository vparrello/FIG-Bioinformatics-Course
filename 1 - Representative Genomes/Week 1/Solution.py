'''
@author Victoria Parrello
Week 1 - Representative Genomes
Activity 1: Find the secret sequences
'''

dna_sequences = "acgt"
header = ""
counter = 0

with open("1 - Representative Genomes/Week 1/Data/Sample1.fasta") as datafile:
    data = datafile.readlines()
    for line in data:
        if ">REP" in line:
            header = line.strip("\n")
        elif "introducing17mers" in line:
            counter += 1
            kmersequence = line.strip("\n")
            kmersequence = kmersequence.strip("acgt")
            print(f'{counter}. {header} has the sequence {kmersequence}.')
            
#Expected Output
# 1. >REP_104628.44_1_covg_50.0 has the sequence introducing17mers.
# 2. >REP_104628.44_1_covg_50.0 has the sequence introducing17mers.
# 3. >REP_104628.44_4_covg_50.0 has the sequence introducing17mers.
# 4. >REP_106588.37_128_covg_50.0 has the sequence introducing17mers.
# 5. >REP_871585.3_286_covg_50.0 has the sequence introducing17mers.
# 6. >REP_797304.7_378_covg_50.0 has the sequence introducing17mers.
# 7. >REP_797304.7_380_covg_50.0 has the sequence introducing17mers.
# 8. >REP_1179773.3_387_covg_50.0 has the sequence introducing17mers.
# 9. >REP_1179773.3_387_covg_50.0 has the sequence introducing17mers.
# 10. >REP_1179773.3_387_covg_50.0 has the sequence introducing17mers.
# 11. >REP_1179773.3_387_covg_50.0 has the sequence introducing17mers.
# 12. >REP_568816.4_391_covg_50.0 has the sequence introducing17mers.
# 13. >REP_568816.4_440_covg_50.0 has the sequence introducing17mers.
# 14. >REP_401562.4_499_covg_50.0 has the sequence introducing17mers.
# 15. >REP_401562.4_542_covg_50.0 has the sequence introducing17mers.
# 16. >REP_391037.6_587_covg_50.0 has the sequence introducing17mers.
# 17. >REP_36809.5_620_covg_50.0 has the sequence introducing17mers.

# 4mer - 4
# this8mer - 8
# Whatabouta15mer - 15
# canigeta16merhere - 16
# introducing17mers - 17
# whatisthepatternofthe26mer - 26