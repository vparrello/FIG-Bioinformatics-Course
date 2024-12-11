# Hammer Exercise 4 - Hammer Quality-Control

**Objective:** Filter out candidate hammers that occur elsewhere
within the complete sequences of the RepGen genomes in addition to the SOURs.

So far we have only looked for Hammer candidates
within the PheS SOUR sequences from our Universe of Genomes;
however, PheS is only one gene out of several thousand genes
within a typical bacterial genome, and it is only about 
1000 basepairs long, whereas a typical bacterial genome is 
several million basepairs long.
Therefore, it is possible that a Kmer that was found exactly once
in exactly one representative genome when only looking at the sequences
for a set of PheS SOURs within a RepGen set
might stilll occur more than once when looking at the full-length 
contig sequences for the members of that RepGen set. 
Hence, once a set of hammer-candidates has been found,
we must doublecheck them to filter out those candidates
that occur elsewhere in any of the representative genomes
in ddition to occuring within the selected set of SOUR genes.

## Generating the filter program

By definition, a `hammer` must occur exactly once in exactly one representative genome,
and that single occurence must be within a SOUR.
When we constructed our set of hammer candidates,
we only looked at the sequences for the SOURs,
not the entire set of contig sequences for the RepGen set,
so now we need to check the hammer-candidates against the contigs
to see if any of them are found elsewhere within one of the genomes.
If after checking all of the contigs sequences for the entire RepGen set,
a given hammer-candidate is still found only once, then it is a valid hammer;
however, if additional instances of a candidate were found elsewhere
in the set of contigs,
then that candidate does not satisfy the definition of a `hammer`,
and it must be eliminated.
We may think of this doublechecking step as a "filtering" process.

Since the number of hammer-candidates is small compared
to the number of Kmers in the set of contigs,
the strategy we will use is to read the candidates into a dictionary,
and then foreach genome in the RepGen set, scan the Kmers in its contigs
to see if any of them are in the dictionary.
If the number of times we have seen a Kmer that is within the dictionary of candidates
ever exceeds 1,
then it cannot be a hammer, and we must delete it from the dictionary.
(It may be helpful at this point to go back and review the interactive exercise on adding and deleting dictionary entries
in `TSV-Exercise-2_Python-Datatypes`.)

We will organize the set of representative-genome contigs
by storing the contigs for each genome as a FASTA file with name-format `genome_id.fna`, 
and store all of contig-files within a directory. 
We read the hammer-candidates from `STDIN`, and load them into a dictionary.
Next, we read the contents of the genome directory,
and foreach genome compare its Kmers against the dictionary of candidates,
removing any dictionary entries that are found more than once within the contigs.
Finally, after checking all genomes in the directory,
we write the remaining contents of the candidate dictionary to `STDOUT`.
Here is a prompt that will generate the algorithm we have just described:

```
Please write a python program named 'filter_hammer_candidates.py' that:

* Accepts a mandatory directory-name argument for a directory
containing a set of FASTA contig files, short name '-C',
long name '--contigs-directory'.

* Reads a hammer-candidate file as a tab-seperated-value file from STDIN,
saves the header-line, and loads the remainder of the file into
a dictionary called the "hammer dictionary",
with the first field as the key, and the entire line as the value.
The keys are named "hammer candidates"; note that all hammers have the same length.
The hammer read-function should return the number of candidates read,
the hammer-length, the header-line, and the hammer-dictionary.

* Reads the contents of the contigs directory.

* Foreach file in the contigs directory,
print the progress-message f"Processing '{filename}'" to STDERR,
where 'filename' is the name of the file, 
and then use BioPython to read the file as a set of 
FASTA-formatted DNA contigs.

* Foreach contig, convert the contig-sequence to lower-case,
and find all Kmers with length equal to the hammer-length
on both the forward and reverse-complement sequences.

* For each Kmer extracted from the contig and its reverse-complement:
  - Check the hammer-dictionary to see if that Kmer is a hammer candidate,
  
  - If the Kmer is a hammer candidate, increment the count for that Kmer
    in the candidate-count dictionary.
    
  - If the count for a candidate exceeds '1', delete that candidate's entry
    from the hammer dictionary, and increment a count of the number
    of candidates eliminated.

* Print the saved header-line to STDOUT, sort the remaining keys
in the hammer-dictionary (now called "Accepted hammers"),
and then foreach accepted candidate in the hammer-dictionary,
print its value to STDOUT.

* Finally, print a blank line to STDERR,
then the number of candidates read,
the number of candidates eliminated,
and the number of hammers accepted, then exit.
```