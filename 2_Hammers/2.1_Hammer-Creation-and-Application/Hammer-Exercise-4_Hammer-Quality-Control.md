# Hammer Exercise 4 - Hammer Quality-Control

**Objective:** Filter out candidate hammers that also occur elsewhere within the complete sequences of a RepGen genome,
in addition to the SOUR that the candidate came from.

So far we have only looked for Hammer candidates
within the PheS SOUR sequences from our Universe of Genomes.
However, PheS is only one gene, and it is only about 
1000 basepairs long; by contrast, a typical bacterial genome
contains several thousand genes, and is several million basepairs long.
Therefore, it is possible that a Kmer that was found exactly once
in exactly one representative genome when one is only looking
at the sequences for the PheS SOURs within a RepGen set
might still occur more than once when one looks at the
complete contig sequences for the members of that same RepGen set
— and by definition, a `hammer` can only occur within
a single RepGen, and it can only occur once within that RepGen.
Hence, once a set of hammer-candidates has been found,
we must doublecheck them to filter out any candidates
that have any additional occurances within any of the representative genomes.

## Generating the filter program

By definition, a `hammer` must occur exactly once
in exactly one representative genome,
and that single occurence must be within a SOUR.
When we constructed our set of hammer candidates,
we only looked at the sequences for the selected SOURs,
not the complete set of contig sequences for the RepGen set,
so we need to check the hammer-candidates against the contigs
to see if any of them are found elsewhere within one of the RepGens before we accept that candidate as a `hammer`.
If, after checking all of the contigs sequences for the entire RepGen set,
a given hammer-candidate is still found only once, then it is a valid hammer;
however, if additional instances of that candidate are found elsewhere within the set of contigs,
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
Once all of the contigs have been scanned,
the remaining contents of the dictionary will be valid hammers, and we can write them out to a file.

(It may be helpful at this point to go back and review the interactive exercise on adding and deleting dictionary entries
in `TSV-Exercise-2_Python-Datatypes`.)

We will organize the set of RepGen contigs
by storing the each genome's contigs in a FASTA file
with name-format `genome_id.fna`, 
and store all of these contig-files within a directory. 
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
containing a set of FASTA contig files, short name '-D',
long name '--contigs-directory'.

* Reads a hammer-candidate file as a tab-separated-value file from STDIN,
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
and the number of hammers accepted, and exit.
```

Ask Grimoire to generate code from this prompt, and then use VScode to save it as usual.

Before we can run the filter, we will need to fetch the contigs for your Repgen from BV-BRC. While it is possible to do this task
directly from the command-line, doing so would require more
knowledge of the command-line syntax than you may corrently have.
Fortunately, we can sidestep this problem by asking Grimoire
to code the task for us! :-)
The following prompt will do the job:

```
Please write a python program named 'get_contigs_for_genomes.py' that:

* Accepts the following two named arguments
  - Filename of a tab-separated-value list containing genome-IDs
    short name '-G', long name '--genome-list'

  - Name of an output directory, short name '-D',
    long-name '--output-data'

* Reads the genome-ID file, skips the header-line,
and loads the remainder of the first column into a list.

* Creates the output directory if it does not already exist,
  else warns the user and exits.

* Foreach genome-ID in the genome-list, prints a progress-message to STDERR,
and then executes the command described by the following format:

   f"p3-genome-fasta --contig {genome_ID} > {Output_Directory}/{genome_ID}.fna"

If the command fails, print an error-message to STDERR, and then exit.
```

Please ask Grimoire to also generate code from this prompt,
and then again use VScode to save it as usual.
We are now in a position to fetch your genome contigs.
You will need to run the following command within
the BV-BRC app's window.
Also, please remember that the BV-BRC app
opens its window in your home-directory,
so you will first need to change directory
to the course directory.

```
python Code/get_contigs_for_genomes.py -G Data/myrep10.genomes.tab -D Data/Myrep10_Genomes
```

(Again, remember that all of the above should be entered on a single line, even if it looks like it's wrapped across multiple lines on your screen.)

The genome-fetching process will take some time,
but every 10 or 15 seconds or so you should see a progress message
indicating which genome is currently being fetched.

Once you have your set of genome contigs,
it's time to run the filtering program:

```
python Code/filter_hammer_candidates.py -D Data/Myrep10_Genomes < Data/myrep10.PheS.hammers.tbl > Data/myrep10.PheS.hammers.filtered.tbl
```
(Again, beware line-wrapping!)

The filtering program will once again take some time,
but it should print progress-messages every few seconds
to let you know that it's running.

Congratulations! You now have a guaranteed valid set of hammers!


### Bonus exercise

Repeat the contig-fetching process and filtering process
of the `myrep50` RepGen set.

## Self-Check

(TDB)
