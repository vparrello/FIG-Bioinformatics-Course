# Hammers Exercise 3 - Using Hammers on a Metagenomic Sample

Objective: Using Hammers to estimate which genomes are present in a metagenomic sample by finding which representatives contributed hammers that hit the sample.

In the last lesson, we use a set of hammers to determine which representative an unknown "Mystery Genome" was closest to.

In this exercise, we will give you a sample containing multiple genomes.
Such mixed samples are often called "metagenomes".
By tabulating which representatives contributed hammers that hit this mixture of genomes, you will be able to estimate which genomes are present in the sample, in a process analagous to taking a "poll" or "survey" of the members of a population.


## Materials

```
FIG-Bioinformatics-Course/
├── 2_Hammers
│   └── 2.1_Hammer-Creation-and-Application/
│       └── Hammer-Exercise-3_Using-Hammers-on-Samples.md (you are here)
└── Data/
    └── MysteryGenome2.fna
    └── rep10.list.tbl
    └── rep10.seqs.tbl
    └── myrep10.hammers.tbl - output from Hammer-Exercise-1
```

## Exercises

1. Recall your `Code/cmd_tsv_select_columns.py`, `Code/hammer_compare.py`, and `Code/hammer_creator.py` tools. Review them for their use case and inputs. 

2. Since we have already created the hammer set, we do not necessarily need to use it again. Using `Code/hammer_compare.py` exclusively allows you to save time in the process. Use this program to get an output for `MysteryGenome2.fasta` to help review your process from the previous lesson.

3. Next we will use this same process to find the hammer counts for a mystery sample from BV-BRC. Review how to use their commmand line tool to download the sample "SRR19064441". Be sure to move it into the `Data` folder so that your programs can more easily access it.

4. Run your `Code/hammer_compare.py` program against the contents of the sample `SRR19064441.fasta`. This will have multiple genomes inside of it so do not eliminate any genomes from consideration when looking at your output.

4. Currently we have a few counts of how many hammers we are expecting to be hit from a sample. But without doing this hundreds of times, it would be very hard to have an accurate picture of how many hammer hits you need to verify its presence. Ask Grimoire to explain to you how a genetic barcode might be inaccurate or inconclusive.

5. Tuning our barcodes so that they only hit our target is crucial to having a useful tool. Without this, we might be casting a wider area than necessary, hitting genomes that are not close, or classifying not enough genomes to be accurate. Ask Grimoire to describe some ways that we can tune our hammers to accurately depict our genome space.

In the next unit, we will explore tuning the hammers to be more precise and accurate by determining how worthy a hammer is.

### Self Check

TODO put answers here
2. Mystery Genome: 
3. Table of genomes in the sample
