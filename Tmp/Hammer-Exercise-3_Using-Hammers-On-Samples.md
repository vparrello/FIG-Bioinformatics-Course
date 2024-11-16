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
    ├── MysterySample1.fna
    ├── myrep10.hammers.tbl   (output from Hammer-Exercise-1)
    └── myrep10.genomes.tbl   (output from RepGen-Exercise-1)
```

## Exercises

1. If necessary, please review your `Code/hammer_creator.py` tool
for its inputs and usage.
If you'd like, you can run it on `Data/MysteryGenome2.fasta`
to help review the process covered in the the previous lesson;
you can find the answer below in the `Self-Check` section.

2. Next we will use this same process to find the hammer counts
for a "Mystery Sample" containing 3 genomes.
Run your `Code/hammer_compare.py` using the hammers from `myrep10`:

```
python Code/hammer_compare.py -H Data/myrep10.hammers.tbl -G Data/myrep10.genomes.tbl < Data/MysterySample1.fna > Data/MysterySample1.rep10.hammer-hits.tbl
```

To check your results, see the `Self-Check` for a table of the Top 5 hits.

Note that 3 representatives yielded hundreds of hits, as expected for a sample known to contain 3 genomes. However, there were also other representatives reported,
albeit only with a small number of hits. We briefly discussed why such weak hits might exist in the previous exercise, but their existance illustrates the need for a minimum threshold of hits before we can be confident that a relative of a representative is present in a sample.

3. Currently you have only seen a few instances of how many hammers hit a single genome or a sample containing multiple genomes. You have seen that some representives hit a given genome or sample with many hammers, while others hit only a few times, and we have argued that the representatives that yield the largest number of hits are most likely to be close relatives of a single genome or a sample-member. However, without the experience acquired by performing hammer-analyses on hundreds of genomes or samples, it may be hard to form an accurate picture of how many hammer-hits you need to relibly infer that a close relative of a representative is present in a sample. Please upload the `Definitions.html` file to Grimoire as usual, and then ask it to explain to you why a 20-character genetic-signature such as a hammer might yield inaccurate or inconclusive evidence for the presence of a close relative of a representatives within a sample, and how one might refine using signatures such as hammers to yield more reliable inferences.

4. Tuning our hammers so that they only hit their intended targets is crucial to having a useful tool. Without this, we might be casting a wider area than necessary, hitting genomes that are not close, or classifying not enough genomes to be accurate. Ask Grimoire to describe some ways that we can tune our hammers to accurately depict our genome space.

In the next unit, we will explore tuning the hammers to be more precise and accurate by determining how worthy a hammer is. We will also discuss how to obtain supporting evidence for the presence of a genome by requiring that the sample contain hits by hammers from more than one SOUR from the same representative.

### Self Check

TODO put answers here
2. Mystery Genome: 
3. Table of genomes in the sample
