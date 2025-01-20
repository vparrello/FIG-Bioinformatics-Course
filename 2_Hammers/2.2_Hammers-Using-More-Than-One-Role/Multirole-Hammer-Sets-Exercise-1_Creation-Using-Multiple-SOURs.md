# Multirole Hammer Sets Exercise 1 - Building Hammers Using More Than One SOUR

**Obejective:** Build a set of hammers starting from a set of sequences that represent more than one SOUR,
to support more reliable RepGen assignments to a sample,
by requiring a "consensus" among the assignemnts made 
by each SOUR.

So far we have built hammers for a RepGenSet
using only the PheS SOUR sequences from the RepGens.
However, we found that in addition to the reliably identified
highest-scoring RepGens found within the sample,
there were also a number of low-scoring "noise hits".

In this set of excercise, we will see that the "noise hits"
can be largely or even completely eliminated by requiring
a "consensus" on which genomes are present within a sample
based on hammer evidence from more than one SOUR. 
You may think of this consensus-building process as being analogous
to forming a "jury" composed of several different SOURs,
and considering their scores to represent "votes"
on which RepGens are closest to the genome or genomes within a sample.
For example, suppose that we have built hammers using sequences
implementing 5 different Singly-Occuring Universal Roles;
we can require that in order for a RepGen to be asserted
to be "close" to some genome within the sample, it must have support
from (for exampl) at least 4 out of the 5 SOURs.
Such a requirement of support from multiple lines
of evidence is fundamental to the scientific process;
for example, we do not consider an experiment to be valid
unless it can be reproduced.

## Materials

```
FIG-Bioinformatics-Course/
├── 2_Hammers
│   └── 2.2_Hammers-Using-More-Than-One-Role/
│       └── Multirole-Hammer-Sets-Exercise-1_Using-Multiple-SOURs.md (you are here)
├── Code/
│   ├── fasta_reader.py
│   └── filter_hammer_candidates.py
└── Data/
    ├── five_roles.tbl
    ├── myrep10.genomes.tbl
    ├── myrep50.genomes.tbl
    ├── Myrep10_Genomes/
    │   └── (141 genome contig-FASTAs)
    └── Myrep50_Genomes/
        └── (921 genome contig-FASTAs)
```

## Programs to be Generated

* get_seqs_for_roles_and_genomes.py

* hammer_creator_multirole.py


## Exercises

### 1. Getting the SOUR sequences

In Hammer-Exercise-1, you fetched the sequences for the PheS SOUR
using a command-line pipeline. While one could in principle
just repeat this procedure 4 more times to get the sequences
for the other 4 SOURs used in this excercise, such a process
will be tedious and error-prone; hence, it is better to ask Grimoire
to automate fetching the sequences for the specified set of roles
and the RepGenSet you will be building hammers for.
As usual, you should upload the file `Definitions.html` to Grimoire
as in previous exercises, and ask it to learn the definitions in the file.
Then, enter the following prompt to generate the automated sequence-fetcher:
```
Please write a python program named 'get_seqs_for_roles_and_genomes.py' that:

* Accepts the following three mandatory named arguments:
  - Type-flag, short name '-T', long name '--type',
    allowed values 'dna' or 'protein'.
    
  - Filename of a tab-separated-value list containing genome-IDs
    short name '-G', long name '--genome-list'.

  - Filename of a tab-separated-value list of roles,
    short name '-R', long-name '--role-list'.

* Reads the role-name file, skips the header-line,
and loads the remainder of the second column into a list.

* Foreach role in the role-list, print a progress-message to STDERR,
and then execute the following two commands, piping the output of 'cmd1'
to the input of 'cmd2', trapping any errors thrown by either command,
and printing the output of 'cmd2' to STDOUT;
note that the role can contain whitespace, and so this argument must be quoted:

   cmd1 = f"p3-get-genome-features --selective --input {genomes_filename} --col genome_id --eq product,'{role}' --attr patric_id"

   cmd2 = f"p3-get-feature-sequence --col feature.patric_id --{type}"

If either command fails, print the trapped error-message to STDERR, and then exit.
```
Ask Grimoire to generate pseudocode for the script if it has not done so,
then use VScode to paste the pseudocode and code into the template
`Code/get_seqs_for_roles_and_genomes.py` and save it.

Because this script calls two "P3 commands",
it will need to be run from within the BV-BRC app's command-window.
To run the command for your 'myrep10' genomes, 
launch the BV-BRC app using the method appropriate
for your computer's operating-system as discussed in
earlier exercises such as `RepGen-Exercise-1`,
enter the `cd` command in the BV-BRC window
to go to the course-material's folder,
(for example, if you installed the course-material
in your `Documents` directory you would enter
`cd ~/Documents/FIG-Bioinformatics-Course`),
and enter the following command in the BV-BRC window:

```
    python Code/get_seqs_for_roles_and_genomes.py -T dna -G Data/myrep10.genome.tbl -R Data/five_roles.tbl > Data/myrep10.five_roles.dna-sequences.fna
```

This program will take some time to run, but it should start printing progress-messages
to the screen to let you know it's working.

### 2. Building the Hammers

Building hammers for a set of roles rather than a single role
will require a slight modification to the prompt you used in Hammer-Exercise-1;
the following revised prompt will generate the desired code:

```
Please write a python program named 'hammer_creator_multirole.py' that will:

* Accept as '-K' a mandatory integer Kmer-length command-line argument;
        
* Read a FASTA-formatted file from 'STDIN' using BioPython,
and foreach FASTA entry extract the following:

  - The portion of the FASTA header up to the first whitespace-character is the
  "feature-ID" ('fid') for the FASTA entry;
  the remainder of the FASTA header starting with the first non-whitespace character
  after the whitespace is that feature's "role".
  There may be more than one whitespace-character separating the feature-ID
  from its role.
  Please extract the feature-ID and the role frome the FASTA header,
  then convert the sequence to lower-case,
  and build dictionaries mapping feature-IDs to sequences,
  and feature-IDs to their roles.

  - Each feature-ID has the format 'fig|x.y.peg.z', where 'x', 'y', and 'z' are integers,
  and the 'fig|' and '.peg' portions are literal strings, not variables.
  The substring 'x.y' is the 'genome_id' for the feature;
  please use a regular expression to extract the genome-ID from the feature-ID,
  and build a dictionary mapping feature-IDs to genome-IDs. 

  - Return the dictionary mapping feature-IDs to genomes,
  the dictionary mapping feature-IDs to roles,
  and the dictionary mapping feature-IDs to sequences.

* Print a tab-separated header-line to STDOUT, whose column-names are
  "hammer", "fid", and "role".

* Create an empty "dictionary of dictionaries" whose outer-keys will be Kmers,
  inner-keys will be feature-IDs, and inner-values will be integer counts.

* Foreach sequence, find all of its Kmers, and foreach Kmer
  increment the count for that sequence's feature-ID
  in the dictionary-of-dictionaries.

* A "hammer" is defined as a Kmer that occur exactly once in exactly one genome-ID.
  Each feature-ID can only occur in one genome.
  So if a Kmer occurs with a count of '1' in exactly one feature, 
  then it must be a "hammer".
  Find all of the "hammers" in the dictionary-of-dictionaries,
  and print to STDOUT a three-column table table of "hammer", "feature-ID",
  and "role" for that feature-ID, sorted by hammer.

* Finally, please print to 'STDERR' the number of sequences that were read,
  the number of Kmers that were processed, and the number of Kmers that were hammers,
  then exit.
```

Again use VScode to paste the pseudocode and code into the template
`Code/hammer_creator_multirole.py` and save it.

To run the command for your 'myrep10' genomes, enter the following:

```
    python3 Code/hammer_creator_multirole.py -K 20 < Data/myrep10.five_roles.dna-sequences.fna > Data/myrep10.five_roles.hammers.tbl
```


## Self-Check

### Exercise 1.

To do a quick check that the correct set of DNA sequences
have been fetched from BV-BRC, run them through the program
`fasta_reader.py` that you created in `FASTA-Exercise-3`:

```
python3 Code/fasta_reader.py < Data/myrep10.five_roles.dna-sequences.solution.fna > /dev/null
```

**NOTE:** In the above command, you are redirecting the data sent to `STDOUT` to a special "device" called `/dev/null`.
`/dev/null` is a program that silently "eats" any input
that sent to it;
if you had left off the `> /dev/null` portion of the command,
you would have seen the names and lengths of all 705 sequences in the FASTA file
before seeing the summary-statistics sent to `STDERR`.
`/dev/null` is useful when you do *not* want to see or save
output from a program for some reason.

The result printed to `STDERR` should be:

```
Number of sequences: 705
Average sequence length: 1823.90
```

Note that "705 sequences" is consistent with fetching
5 roles each for 141 genomes.

### Exercise 2.

The number of lines in the file `Data/myrep10.five_roles.hammers.tbl` should be 1259712;
you can check this as follows:

```
wc -l Data/myrep10.five_roles.hammers.tbl
```

Recal that `wc` (short for "word count") is the command
that counts the number of characters, words, and lines in a file,
and that including the `-l` option tells the program
to only returns the number of lines, not the number of characters and words.

For completeness, we have included solutions-files
for the fetched sequences, the sets of hammers,
and the sets of filtered hammers
in the `Solutions/` subdirectory for this module,
however they are in `gzip-compressed form`
because they were too large to upload to GitHub,
so in order to look at them you will need to "uncompress" them.
gzip-compressed files can be uncompressed from the command-line
by using the `gunzip` command;
however most operating-systems support uncompressing
a commpressed file simply by "double-clicking" on them
in your file-browser.

We have also included soultion "hammer-hits" files
for each of the test-cases in this exercise.
