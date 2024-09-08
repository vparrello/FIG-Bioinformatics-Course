# Hammer Precision Exercise 1

Objective: Understand the concept of hammer precision and learn how to evaluate the effectiveness of genetic barcodes (hammers) in identifying specific genomes.

Hammers are unique k-mers that occur exactly once in a genome. The precision of a hammer refers to how specific it's sequence is to the genomes inside of the representative genome's space. A highly precise hammer will only match a single genome or two, while a less precise hammer may match all of the genome peers inside that genome space. This exercise will help you understand how to evaluate and improve the precision of your hammers by comparing them to the other hammers inside of the hammer set. 

## Materials:
[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

FIG-Bioinformatics-Course/
├── 2_Hammers
│   └── 2.2_Hammer-Tuning/
│       └── Hammer-Tuning-Exercise-1_Precision.md (you are here)
└── Data/
    ├── rep10.hammers.tbl (from previous exercise)
    ├── rep10.seqs.tbl
    └── test_genomes.fna (new file with multiple test genomes)
    Code/
    └── find_duplicate_hammers.py
    └── hammer_precision.py

## Exercises:

1. Hammers have multiple functions. One of the most important things is that they allow us to identify genomes. If a genome shares a hammer with its representative genome, we know that it is a neighbor of that genome. But sometimes a genome may share a hammer with another representative genome. In this case, we need to know if this occurs since it can skew our data to multiple genomes inside the hammer set. Ask Grimoire how it would determine if duplicate hammers exist in the `rep10.hammers.tbl` file.

2. We need to evaluate the precision of our hammers. Ask Grimoire to write a Python script named `hammer_precision.py` that does the following:
   - Reads the `rep10.hammers.tbl` file
   - For each hammer, create a counter of how many times it appears in the `rep10.hammers.tbl` file
   - Iterates through each hammer and deletes any hammers that have a counter of 1.
   - Outputs a table with columns: hammer, genome_id, counter

3. Run your `hammer_precision.py` script and analyze the results. Are there any genomes that show up more often than others? Ask Grimoire to explain to you why a genome might have more duplicate hammers than others and what it means for the accuracy of the genome inference.

4. Having more precise hammers means that you can more accurately determine where in genome space your sample is from. However, as you take hammers out of your set, you are decreasing the possibilities of finding a matching hammer for your sample. Making our hammer set more like a chisel than a sledgehammer. Ask Grimoire to give you a use case of when you would want to have more precise hammers vs when you would want to have more hammers in your set.

5. We need to create a new set of hammers that are more precise. Ask Grimoire to write a modify the Python script named `hammer_precision.py` that does the following:
   - Reads the `rep10.hammers.tbl` file
   - For each hammer, create a counter of how many times it appears in the `rep10.hammers.tbl` file
   - Iterates through each hammer and outputs any hammer that has a counter of 1 to a new file called `rep10.precise.hammers.tbl`
   - Outputs a table with columns: hammer, genome_id, counter
   - Takes a command line argument for how many duplicates of a hammer are allowed before it is removed from the set.

6. Run your `hammer_precision.py` script on the 'rep10.hammers.tbl' file with a command line argument of 5. Look at the difference in size between the original 'rep10.hammers.tbl' file and the new 'rep10.precise.hammers.tbl' file. Giving some leeway to the number of duplicates allowed, can sometimes help to improve the performance of your hammers for different use cases. Ask Grimoire to explain to you how this might work.

7. Bonus: Test out different values for the number of duplicates allowed in the `rep10.precise.hammers.tbl` file and use the different hammer sets to identify the genomes in the `test_genomes.fna` file. See if there is a relationship between the number of duplicates allowed and the performance of the hammers.

In the next exercise, we will explore additional factors that contribute to hammer performance and how to incorporate them into our hammer selection process.

