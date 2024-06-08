# Building "Representative Sets"

The reduction in cost of genetic-sequencing
has led to an explosion in the availability
of genetic data.
The public databases now contain roughly
1 Million genomes, and several Billion gene
and protein sequences. Given the huge abundance
of available sequence-data, there is thus a need
to extract subsets of tractable size that
summarize the diversity of available genomes,
gene, and protein sequences;
we call such a subset a "Representative Set".
We call a set of "Representative Genomes"
a "RepGen Set" for short, and a generic set
of representatives a "RepSet".

The goal of this first module is to construct
a set of representatives from a set of sequences
or genomes. However, we will first need to introduce
some basic concepts, such as working with the command-line
interface, the basic types of data used in bioinformatics,
and the use of a "Large Language Model" as a programming assistant.


## 1.1 Tab Separated Value Files

### TSV-1 - Learning to use Grimoire
Become familiar with ChatGPT and interacting
with the Grimoire AI model

### TSV-2 - Python Datatypes
Learn the different types of data that are supported by the Python
programming language

### TSV-3 - Working with TSVs
Create a program that reads a "tab-separated value" (TSV) data-file
and extracts selected data-fields from it

### TSV-4 - Modifying Code
Take your previous program and adjust it to support a new use case

### TSV-5 - Fixing Syntax Errors
Find out how to fix your code when it breaks

### TSV-6 - Documentation and its Uses
Write down how your program works and how to optimally interact with it 

## 1.2 FASTA Files

### FASTA-1 - What is FASTA Format?
Learn how bioinformatic sequence-data is stored 

### FASTA-2 - TSV to FASTA conversion
Mapping between the two data formats that we've used so far

### FASTA-3 - Reading FASTA Files
Create a program that reads and operates on FASTA-formatted sequence-data

### FASTA-4 - DNA to Protein Translation
Learn how DNA sequences are related to protein sequences,
and how to apply the genetic code that maps DNA into
proteins.

## 1.3 Kmers and Jaccard Similarities

### Kmers-1 What Are Kmers? -
Using DNA subsequences as proxies for an entire sequence

### Kmers-2 Jaccard Similarities -
Using Kmers to estimate how similar two sequences are

### Kmers-3 - Finding the Nearest Reference-Sequence
Given a set of "query" and "reference" sequences,
for each query find the closest reference sequence

### Kmers-4 - Protein vs DNA Kmers
How are protein and DNA similarities related?

## 1.4 - Building Representative Genome-Sets