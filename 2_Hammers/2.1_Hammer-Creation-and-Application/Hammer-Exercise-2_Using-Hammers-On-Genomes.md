# Hammer Exercise 2 - Using Hammers to find the Nearest Representative

Objective: Given a set of representative genomes,  a set of associated "Hammers", and an unknown genome, determine which representative the unkonwn genome is most similar to.

In this exercise we will generate code that will search the "Mystery Genome" to see if it contains any "hammers". If the "Mystery Genome" contains many of the same "hammer" barcodes as one of the representative genomes, and very few "hammers" for any other representative genome, then this indicates that there is a high probability that the "Mystery Genome" and the representative genome that contributed the largest number of "hammer hits" are related.

## Materials

```
FIG-Bioinformatics-Course/
├── 2_Hammers 
│   └── 2.1_Hammer-Creation-and_Application/
│       └── Hammer-Exercise-2_Using-Hammers-on-Genomes.md (you are here)
└── Data/
    ├── MysteryGenome1.fna
    └── myrep10.PheS.hammers.tbl
```

## Exercises

1. We need to compare the DNA sequences in `MysteryGenome1.fna`
to the set of hammers that you generated in the previous exercise.
Ask Grimoire to write a program named `hammer_compare.py` that will:

    * Accept a mandatory TSV "hammers" filename argument,
    short name `-H`, long name `--hammers`.

    * Accept a mandatory TSV genome-names filename argument,
    short name `-G`, long name `--genome-names`.

    * Skip the header-line of the hammer-file and then read
    the first and second columns as a `hammer` and `feature_id`,
    respectively. A `feature_id` has the format 'fig|x.y.peg.z',
    where 'x', 'y', and 'z' are integers, and 'fig|' and '.peg.'
    are literal substrings. The portion 'x.y' is the `genome_id`
    for the 'feature_id'; extract the `genome_id` using a regex,
    and build a dictionary mapping each hammer to its genome_id.

    * Skip the header-line of the genome-names file and then read
    the first and second columns into a dictionary as a `genome_id`
    and `genome_name`, respectively.

    * Determine the Kmer-length `K` of the hammers from the hammers dictionary.
    (NOTE: all the hammers in the file will have the same length.)

    * Use BioPython to read the sequences of the genome from `STDIN`.

    * For each sequence, extract all possible Kmers, and if a Kmer is a hammer,
    increment the score for its associated `genome_id`; then repeat this operation on the reverse-complement of that sequence, since a gene can face in either direction.

    * Print to STDOUT a TSV file of the genome_ids found
    and their associated genome_name and score
    sorted by decreasing score.
    Please handle missing genome-names gracefully;
    if a genome_id does not have an associated genome-name,
    display the genome_name as 'Unknown sp.' in the output TSV file,
    and send a warning to STDERR that the name of genome_id was not in
    the genome-names file.
    
Finally, please ask Grimoire to translate the python code into pseudocode.

2. Once Grimoire is done, please paste its pseudocode and code into `Code/hammer_compare.py` and save it as usual.

3. Run your script on the `MysteryGenome1.fna`. The result should be a table of representative genome-IDs followed by how many hammers from that representative hit the "Mystery Genome". Take the genome-ID with the top number of hits as the most probable representative for the "Mystery Genome".


Congratulations! You have constructed a set of hammers for the PheS SOURs in the `myrep10` set of representative genomes, and then used them to identify the `myrep10` genome `511145.12` as the most likely representative for your "Mystery Genome"!

## Caveats

In this exercise we have presented a greatly simplified "toy version" of the the hammer-creation process to illustrate the basic concepts behind hammers. The "production" version of the code is significantly more complex; we use the sequences for the 20 SOURs that we have found to yield the most reliably informative hammers, and we also cross-check the entire genome in addition to the SOUR sequences to make certain that the hammer does not occur somewhere else within the genome, since by definition a hammer can only occur once within one representative genome and must be contained within a SOUR. We have avoided these additional complexities within this exercise in order to keep this exercise simple, but in a later exercise, you will refine your prompt to take these additional complexities into account.


## Self Check

Below are the list of genomes in the representative genome set. Use the counts for the mystery genome to check your work.

The "Mystery Genome" should have `511145.12` "Escherichia coli str. K-12 substr. MG1655" as its highest-scoring representative genome.

### NOTE: Table of counts is still missing !!!