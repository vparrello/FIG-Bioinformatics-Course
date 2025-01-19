# Hammer Exercise 4 - Building Hammers Using More Than One SOUR

**Obejective:** Build a set of hammers starting from a set of sequences that represent more than one SOURs,
and make more reliable genome assignments by requiring a "consensus" for the assignemnts based on each SOUR.

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
from at least 4 out of the 5 SOURs. Such a requirement of support
from multiple lines of evidence is fundamental to the scientific process;
for example, we do not consider an experiment to be valid unless
it can be reproduced.

## Materials

```
FIG-Bioinformatics-Course/
├── 2_Hammers
│   └── 2.1_Hammers-Creation-and-Application/
│       └── Hammers-Exercise-1_Using-Multiple-SOURs.md (you are here)
├── Code/
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

* hammer-compare_multirole.py


## Getting the SOUR sequences

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

To run the command for your 'myrep10' genomes, enter the following:

```
    python Code/get_seqs_for_roles_and_genomes.py -T dna -G Data/myrep10.genome.tbl -R Data/five_roles.tbl > Data/myrep10.five_roles.dna-sequences.fna
```

This program will take some time to run, but it should start printing progress-messages
to the screen to let you know it's working.

## Building the Hammers

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

## Using the Hammers

The hammer-application program also needs revision to support multiple roles;
the following prompt will do it:

```
Please write a program named `hammer_compare_multirole.py` that will:

* Accept a mandatory TSV "hammers" filename argument,
short name `-H`, long name `--hammers`.

* Accept a mandatory TSV genome-names filename argument,
short name `-G`, long name `--genome-names`.

* Accept an optional "minimum fraction of roles required" argument,
short name `-F`, long name `--role-fraction`, default-value 0.8.

* Open the hammer-file and perform the following:
    
    - Read the header-line and extract the column-names
    
    - Foreach subsequent line, read the "hammer" and "fid"
    columns as `hammer` and `feature_id` variables,
    respectively, and built a dictionary mapping
    hammers to feature_ids.

    - If a "role" column exists, read it as the `role` variable,
    else the role should default to "Unknown",
    and build a dictionary mapping feature_ids to roles.

    - A `feature_id` has the format 'fig|x.y.peg.z',
    where 'x', 'y', and 'z' are integers, and 'fig|' and '.peg.'
    are literal substrings; the portion 'x.y' is the `genome_id`
    for that 'feature_id'. Extract the `genome_id` from the `feature_id`
    using a regex, and build a dictionary mapping feature_id to its genome_id.

    - Determine the Kmer-length `K` of the hammers.
    (NOTE: all the hammers in the file will have the same length.)

    - Return the hammer-length `K`, the list of roles,
    the hammer-to-feature_id dictionary,
    the feature_id-to-genome_id dictionary,
    and the feature_id-to-role dictionary

* Open the 'genome-names' file, skip the header-line,
and for subsequent lines read the first and second columns
into a dictionary as a `genome_id` and `genome_name`, respectively;
then return the genome_id-to-genome_names dictionary.

* Create an empty directory whose keys will be genome-IDs,
and values will be an integer count of the number of hammers
found for that genome.

* Create an empty "dictionary of dictionaries" whose outer-keys will be genome_ids,
inner-keys will be roles, and inner-values will be the integer count
for the number of hammers found for that role in that genome.

* Use BioPython to read the sequences of the genome from `STDIN`;
the sequence should then be converted to lower-case.

* For each sequence, extract all possible Kmers, and if a Kmer
is a hammer, increment the hammer-count for its associated `genome_id`
in the genome-counts dictionary, and also the hammer-count for that role
in the genome-to-role-to-counts directory-of-directories;
then repeat this operation on the reverse-complement of that sequence,
since a gene can face in either direction.

* Once all the sequences have been processed,
foreach genome_id in the genome-to-roles-to-counts directory
of directories, if the number of roles for that genome
is greater than or equal to the total number of roles
times the minimum fraction of roles,
print to STDOUT a TSV file of the genome_ids found
and their associated genome_name and score,
sorted by decreasing genome-to-counts score.
If the total number of roles exceeds 1,
also add a column for the number of roles for the genome.
Please handle missing genome-names gracefully;
if a genome_id does not have an associated genome-name,
display the genome_name as 'Unknown sp.' in the output TSV file,
and send a warning to STDERR that the name of genome_id was not in 
the genome-names file.
```

To run the program on 'MysteryGenome1.fna', enter the following:

```
python3 Code/hammer_compare_multirole.py -H Data/myrep10.five_roles.hammers.tbl -G Data/myrep10.genome.tbl < DataMysteryGenome1.fna > Data/MysteryGenom1.myrep10.five_roles.hammers.hits
```

The resulting output table is:

| Genome_ID | Genome_Name | Score | Roles_Found |
| --- | --- | --- | --- |
| 511145.12 | Escherichia coli str. K-12 substr. MG1655 | 8257 | 5 |

Recall that in `Hammer-Exercise_1`, the program `hammer_compare.py` found several low-scoring "noise genomes" in addition to the correct RepGen `511145.12`, and that in order to remove the "noise hits" we had to "filter" the hammers to remove candidates that had additional occurances outside the SOUR sequences. However, using multirole hammers and requiring at least 4 out of 5 roles to be found before we will accept a RepGen hit, the program returns only correct RepGen and no "noise genomes" even without filtering!
So requiring multiple confirming-evidences for a repGen does appear to eliminate the "noise problem", even without "filtering". We can confirm the preceeding claim by reruning the program with the minimum role-fraction set to 0
so that all hits will be accepted;  the resulting top 5 hits are:

| Genome_ID | Genome_Name | Score | Roles_Found |
| --- | --- | --- | --- |
| 511145.12 | Escherichia coli str. K-12 substr. MG1655 | 8257 | 5 |
| 1637999.4 | Verrucomicrobia bacterium IMCC26134 | 4 | 1 |
| 2024963.9 | Caldimicrobium sp. strain SpSt-122 | 3 | 1 |
| 1062.102 | Rhodobacter sp. strain HK-STAS-PROT-59 | 3 | 2 |
| 309798.4 | Coprothermobacter proteolyticus DSM 5265 | 3 | 1 |

Note that the low-scoring "noise genomes" only had hits for 1 or 2 roles, so requiring a "4 out of 5 jury" eliminates them.


