# Kmer Exercise 1 - What are Kmers?

Objective: Learn about the concept of "Kmers" in sequences.

A DNA sequence is a sequence composed from the 4-character "alphabet" {a,c,g,t}. Each character is an abbreviation for one of the "nucleic acids" that are the subunits of a DNA molecule: Adenine, Cytosine, Guanine, and Thymine. These characters are also called "bases" by biologists (short for "nucleobases"). The four "bases" form pairs, Adenine with Thymine and Cytosine with Guanine, and each pair forms a "rung" on the twisting double-helix "ladder" of a DNA strand --
see [this picture](https://en.wikipedia.org/wiki/DNA#/media/File:DNA_chemical_structure.svg),
[this picture](https://en.wikipedia.org/wiki/DNA#/media/File:DNA_Structure+Key+Labelled.pn_NoBB.png)
and [this picture](https://en.wikipedia.org/wiki/DNA#/media/File:DNA_animation.gif).
Since in DNA, Adenine always pairs with Thymine and Cytosine always pairs with Guanine, it is sufficient to only list one of the two strands in the "double-helix", since given one strand the other (or "complementary") strand can be trivially determined via base-pairing.

Depending on the organism, the DNA sequences composing its genome can be quite long; a virus will typically have a genome of between a thousand and 100 thousand base-pairs, a bacterium or archaeum will have a genome of between several hundred thousand to several million base-pairs, and a large multicellular organism such as a mouse or a human will have a genome containing billions of base-pairs.

Proteins are also represented by sequences of characters, but they instead use a 20-character "alphabet", with each of the 20 characters standing for an "amino acid" instead of a "nucleic acid".

In order to make searching for patterns in DNA or proteins more manageable, bioinformaticians will often break up long sequences into a set of short sequences called "Kmers". A "Kmer" is a sequence or string of length "K"; for example, `atgcgt` would be an example of a DNA "6-mer". 

In this course, we will often be working with DNA 20-mers and protein 8-mers.
* Review-question: Remember when we did FASTA files? Can you remember which file-extensions were commonly used for DNA vs Proteins?

Kmers have many uses in bioinformatics; not only can they be used to simplify searching for patterns in DNA or protein sequences, but they can be used to quickly identify sequences and to define a measure of "similarity" between two sequences. In this lesson, we will explore how bioinformaticians construct a set of Kmers from a sequence, and examine some of their basic use-cases.


## Materials: 

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)
<!--
* FIG-Bioinformatics-Course/
    * Definitions.html
    * 1_Representative-Genomes/
        * 1.3_Kmers-and-Jaccard-Similarities/
            * Kmer-Exercise-1_What-are-Kmers.md
            * Solutions/
                * Kmer-Exercise-1_Solutions.md
-->

```
FIG-Bioinformatics-Course/
├── Definitions.html
└── 1_Representative-Genomes/
    └── 1.3_Kmers-and-Jaccard-Similarities/
        ├── Kmer-Exercise-1_What-are-Kmers.md  (You are here)
        ├── Code/     (where you will save your own code)
        ├── Data/
        │   ├── test_dna.fna
        │   └── test_protein.faa
        └── Solutions/
            ├── Kmer-Exercise-1_Solutions.md
            └── extract_kmers_from_fasta_solution.py
```

## Exercises:

Exercise Setup -- Please enter the following into Grimoire's "Message Grimoire" box:
  ```
    I am going to give you some definitions in an attached file;
    you don't need to respond to them, just learn them.
    I'm then going to ask you some questions.
  ```
Then click on the "Paperclip" icon and select the file `Definitions.html`, and hit the "Send Message" icon (the "Up-arrow" icon at the right of the "Message Grimoire" box). Grimoire should respond to your entering the above with something like "Understood" or "Got it". `Definitions.html` contains a list of definitions of terms and concepts specific to this course. We will be using the preceding method to "pre-prompt" Grimoire regarding the "context" for a set of exercises throughout the remainder of this course; you are welcome to look at the file `Definitions.html`, but some of the terms, definitions, and concepts may not make sense in advance of the exercise they are relevant to.

The following questions should be entered as individual prompts (or as Grimore puts them, "Messages"), but should be entered during the same session in which you uploaded the above definitions, since Grimoire will need these definitions to establish its "context", and it cannot remember the contents of one session in a different dession, but will instead fall back on its "default" definitions unless "re-prompted".

1. Ask Grimoire "What is a Kmer, and how are Kmers used in bioinformatics?"

2. "How can I know how many Kmers are in a given sequence?"
    * (Grimoire should answer with an explanation or formula that calculates the number of Kmers given the length of the sequence and `K`)

3. "Let 'cttgtatagtattgtcttgtgtatc' be a DNA sequence; please tell me the length of this sequence, and then list all of the 20-mers within this sequence."
    * NOTE: Do please check how many characters are in Grimoire's answers, because sometimes it will miscount the number of characters in a sequence -- e.g., it may get the length of the sequence wrong, or return a list of "19-mers" instead of "20-mers". Grimoire is not perfect! But it's pretty good at what it does do well, which is write code and explain concepts. If Grimoire miscounts the length of the entire sequence (which in this case should be '25') or the length of the Kmer, ask it "Are you sure about those lengths?" and it will usually correct its error.

4. "Is the Kmer 'atagctcga' within the sequence 'cttgtatagtattgtcttgtgtatc'?"

5. "Is the Kmer 'agtattgtc' within the sequence 'cttgtatagtattgtcttgtgtatc'?"

6. Ask Grimoire to write a program named `extract_kmers_from_fasta.py` that satisfies the following specifications:
  * Mandatory Kmer-length command-line argument short-form `-K`, long-form `--Kmer`, which must have an integer value
  * Optional command-line argument short-form `-t`, long-form `--type`, whose values may either be `dna` or `protein`
  * If invoked with an argument `-h` or `--help`, the program should print a "help" message to STDERR describing the program, listing the argument names, showing whether an argument is mandatory or optional, and indicating that the data-file will be read from STDIN, and should then exit
  * The program shall read a FASTA-formated file from STDIN; you are allowed to use BioPython to read, write, and operate on sequences.
  * if the `--type` argument is specified, then if `--type` is `dna` the sequence should be converted to lower-case, else if `--type` is protein the sequence should be converted to upper-case.
  * For each sequence, the program should extract all possible Kmers of length given by the Kmer-length argument, and write to STDOUT a two-column tab-separated file consisting of the Kmer and the sequence-ID that it came from .
  * If the `--type` argument is specified, then the Kmer should be checked before output for whether it contains characters that are not in the DNA or Protein alphabets, respectively, and invalid Kmers and the sequence-id they came from should be written to STDERR instead of STDOUT.
  * On exit the program should report to STDERR the number of sequences read, the number of Kmers written to STDOUT, and the number of invalid Kmers found.

Use VScode to save `extract_kmers_from_fasta.py` in the `Code/` subdirectory.

7. Run `extract_kmers_from_fasta.py` on the following files, which will be found in subdirectory `1.3_Kmers-and-Jaccard-Similarities/Data/`. These files are all small, so you should be able to figure out what the output should be in your head. Each of the files will contain one "bad" sequence, i.e., a sequence that contains invalid characters for that file-type.

* test_dna.fna
* test_protein.faa

Experiment with different values for `-K`. Then, see how the output changes  when you specify the `--type` argument.
* BONUS: What do you think will happen if you use the 'dna' type for the protein file, or the 'protein' type for the DNA file? Try it, and see if you have guessed correctly.

## Solution Check instructions:



## NOTES:
* Once we show Grimoire as fallible, should we start the pseudocode?
* We need a self check for terms like: kmer, dna, protein
* tool that reads a fasta and spits out all the kmers that appear in the file
* write a tool that reads a FASTA file and writes a two column file with the kmer and the sequence it came from
  - Spit out to STDERR any kmers that have a bad character (not inside agct)
