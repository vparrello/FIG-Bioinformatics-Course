# Hammers as Genetic Barcodes Exercise 2 - New Genome and his Representative

Objective: Find a "Genetic Barcode" for a genome in order to classify one inside a sample.

Normally when working with bioinformatics data, we are given a genome that has been annotated. This means that it has been compared to the rest of genome space and has been sorted into the most likely genomes for that sample. If a sample has not been annotated, then we need to use our hammers to try and find them.

In this exercise you will be given a sample that has not been annotated. We will then walk you through the steps of finding the hammers from a single representative genome to see if any of those hammers are present in the mystery genome. If they are, that tells us that they also have the same barcodes in their genome and therefore are more likely to be the same genome.

## Materials

TODO insert a sample data here with planted genomes in it. -  MysteryGenome.fasta
TODO have a sample representative genome sample for creating hammers. - Rep10.list.tbl, Rep10.seqs.tbl, genomeid.fasta
TODO program that will grab 20-mers without killing the memory in most machines

## Exercises

1. First we need to make a set of hammers that can be used as "barcodes" for at least 1 representative genome. We will do this by taking out all the unique 20-mers from the PheS inside of Rep10.seqs.tbl. Ask Grimoire to write a program that will take a string and create all the unique substrings of that sequence and then print a two column table with the hammers and their representative genome as the two columns. Paste that program into 'hammer_creator.py'.
    *Hint: Can you remember which datatype can only have unique strings?

2. Next we need to connect that program to the tsv_headers.py program because all of our PheS sequences are in Rep10.seqs.tbl. Ask Grimoire to adjust the tsv_headers.py program to take in a genome name and the file name for the Rep10.seqs.tbl file as command line arguments. It then needs to grab the PheS sequence for that genome and send it to the hammer_creator.py program. 

3. A lot of the times, people will use a Command Line Script (or a bash script) to call more than one program. Other times people will use Object Oriented Programming to connect two programs together through the code itself. Ask Grimoire for explanations of both frameworks and decide which one you are more comfortable with. Then tell Grimoire to adjust tsv_headers.py and hammer_creator.py to follow that framework.

4. Once you have those programs connected, it is time to create hammers for a specific genome. Using the framework you chose in step 3, create hammers for the first genome in the list. 

5. We now need to compare these hammers to the MysteryGenome. Ask Grimoire to make another program. One that will read a sequence, make a 20mer, compare it with the set, and then add to a counter every time it has a match. Paste this code into hammer_compare.py.

6. Add hammer_compare.py to your framework from Step 3. If you need help, ask Grimoire to adjust the format to include all the programs. 

7. Run your script on the MysteryGenome.fasta. You should have a count print out of how many hammers were hit. 


Congratulations! You have created the genetic barcodes from the PheS SOUR for the TODO REPGEN NUMBER HERE genome. We need to do the same for the rest of the representative genome set and see how the counts differ. Please complete the hammers for all of the representative genomes so that you are prepared for the next exercise.


### Self Check

Below are the list of genomes in the representative genome set. Use the counts for the mystery genome to check your work.

TODO put answers here