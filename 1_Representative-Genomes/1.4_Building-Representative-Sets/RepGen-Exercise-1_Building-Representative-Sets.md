# RepGen Exercise 1 - Building a Representative Set

Objective: Combining what we have learned to create a program that builds a set of Representative Sequences from an input FASTA file.

## Materials

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

```
FIG-Bioinformatics-Course/
├── 1_Representative-Genomes/
│   └── 1.4_Building-Representative-Sets/
│       └── RepGen-Exercise-1_Building-Representative-Sets.md (you are here)
└── Data/
    └── Universe.fasta
```
## Set up:

At the beginning of this course, you downloaded the BV-BRC application. It should have created a desktop shortcut that looks like a terminal window. The following exercises will require that you use this application. 

Under `macOS` you can open the BV-BRC app by either double-clicking on its icon
in your `/Applications` folder, or by using the following command-line method:
```
open /Applications/BV-BRC.app
```

Under `Windows`, you can double-click on the shortcut icon that the installation-wizard should have created on your desktop.
=======
Windows Example:
```
start "C:\Users\Default\AppData\Local\bv-brc.exe"

```

## Overview:

As a result of the rapidly decreasing cost of sequencing genetic data, we are now at the stage where on the order of a million  bacterial genomes have been sequenced and deposited in the public databases along with several billion associated gene-sequences, leading to an abundances of riches that have become difficult to manage.
There has thus developed a need to construct manageably small sets of "representatives" that capture the diversity of the full set of genomes and genes.

We define a "Set of Representatives" as follows:

* Let U be the Universe of entities being considered (in our case, genomes or sequences).

* Let there be some measure of similarity between two entities.

* Then a "Set of Representatives" or "RepSet" is a subset of U such that:

    - no two members of the RepSet are more similar than some threshold similarity `Sim`, and

    - every member of U is similar to at least one member of the RepSet.

When the set of entities being represented is a set of genomes, we usually refer to it as a "RepGen Set" for short, and a member of the "RepGen Set" is called a "RepGen". In this exercise, we will look at the simpler problem of constructing a set of "Representative Sequences".

Many measures of similarity can be used, but the measure we will use in this exercise is "Number of Kmers in Common". 

* NOTE: In previous exercises, we have worked with data-files having names like `rep10` or `rep200`; each of these files were obtained by constructing a set of representative genomes, and the number in the filename was the threshold number of Kmers in common. For example, no two sequences in `rep10` have 10 or more short Kmers in common; no two members of `rep200` have 200 or more short Kmers in common, etc.
In our "Production" RepGenSets such as `rep10` or `rep200`, we have selected the protein-sequence that implements the role "Phenylalanine tRNA Synthetase alpha subunit" (abbreviated as "PheS") as the sequence to base our RepGenSets on, and "number of protein 8-mers in common between two sequences" as our measure of similarity.
"PheS" is an example of a "Singly-Occurring Universal Role" (SOUR), which is a role implemented by exactly one gene (it is "Singly Occuring"), and that can be expected to be found in every genome (it is "Universal").

There are several algorithms for building RepSets. In this course, we will be using the "Stingy Addition" algorithm. The pseudocode for "Stingy Addition" is as follows:
```
  RepGenSet = ()   #... i.e. the "empty list"
  ToCheck   = U    #... The "Universe of Genomes"
  Foreach G1 in ToCheck:
    If no G in RepGenSet is more similar to G1 than MinSim:
      Add G1 to RepGenSet
            	
  Return RepGenSet
```

## Exercises:

1. Attach the `Definitions.html` file as in previous exercises.

2. The following prompt provides the program specification for the representative-set algorithm:

```
I will now give you a description of a command-line interface and
program called `build_representative_set.py` that will compute a
"set of representative sequences" (RepGen set) using the "Stingy Addition"
algorithm as defined in the uploaded definitions-file.

The script should accept the following mandatory command-line
arguments, in both long-form and the specified short-form:

* Kmer-length  (integer, short argument-name '-k')
* Sim          (similarity threshold, short argument-name '-s')
* input FASTA-file   (filename, short argument-name '-f')
* output RepSeq-file  (filename, short argument-name '-r')

The measure of similarity to be used is "Number of Kmers in common".

The program should use BioPython to read the FASTA-file.
The first nonwhitespace portion of each FASTA identifier is a "feature-ID",
while the remainder is the "genome name".
Feature-IDs have the format 'fig|X.Y.peg.Z', where 'X', 'Y', and 'Z' are integers,
and the portions 'fig|' and '.peg.' are literal substrings, not variables.
The subpattern 'X.Y' within the feature-ID is a "genome-ID";
you may extract this subpattern using a regular expression.
Return a dictionary that maps feature-IDs to genome-IDs,
a dictionary that maps feature-IDs to genome-names,
and a list of (feature-ID, sequence) pairs.

The main body of the program should construct a subset of the input sequence list
that satisfies the provided definition of a "Representative Set" (RepGen set).

When you are done, please write out the representative set to the RepSeq-file in FASTA format,
where the FASTA identifiers have the form "genome-ID genome-name".
```

Save the program that Grimoire generates as `build_representative_set.py`.

3. Run the program as follows:

```
python Code/build_representative_set.py -k 8 -s 10 -f Data/Universe.fasta -r Data/myrep10.faa
```

4. Please repeat the previous step with a similarity-threshold of 50,
as we will need this result in a later exercise:

```
python Code/build_representative_set.py -k 8 -s 50 -f Data/Universe.fasta -r Data/myrep50.faa
```

## Self-Check

### myrep10

Check your `myrep10` result by running:
```
python Code/fasta_reader.py < Data/myrep10.faa > Data/myrep10.genomes-and-lengths.txt
```
The result should print the following to STDERR (which in this case will just be the screen):
```
Number of sequences: 153
Average sequence length: 356.73
```
You can make a detailed comparison of your results with the solution results using the `diff` command, which compares two files:
```
diff Data/myrep10.genomes-and-lengths.txt 1_Representative-Genomes/1.4_Building-Representative-Sets/Solutions/myrep10.genomes-and-lengths.solution.txt
```
since there should be no difference between these two files, the `diff` command should emit nothing. (NOTE: once again, the above will probably appear to be "wrapped" onto multiple lines on your screen, but it should all be entered as a single line.)

### myrep50

Similarly, you can check your `myrep50` result by running:
```
python Code/fasta_reader.py < Data/myrep50.faa > Data/myrep50.genomes-and-lengths.txt
```
The result should print the following to STDERR (which in this case will just be the screen):
```
Number of sequences: 1036
Average sequence length: 352.23
```
You can make a detailed check by running:
```
diff Data/myrep50.genomes-and-lengths.txt 1_Representative-Genomes/1.4_Building-Representative-Sets/Solutions/myrep50.genomes-and-lengths.solution.txt
```
Again, `diff` should return nothing if there is no difference between the two files.
