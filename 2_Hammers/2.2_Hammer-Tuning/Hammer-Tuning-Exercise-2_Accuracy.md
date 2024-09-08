# Hammer Accuracy Exercise 2

Objective: Understand the concept of hammer accuracy by incorporating more SOURs in the creation process.

When looking at the accuracy of the hammers, we are trying to determine which genomes are inside of the sample with as little doubt as possible. We do this by looking at the universal k-mers between the representative genome and its neighbors. If we use more SOURs to create the hammers, we will have more universal k-mers to confirm the identity of the representative genome for that sequence. This will help us determine which genomes are inside of the sample with as little doubt as possible.

This exercise will help you to understand how we can guarantee the accuracy of the hammers. Nothing in biology is certain, but we can get as close as possible using an abundance of evidence in the form of hammers.

## Materials:
[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

FIG-Bioinformatics-Course/
├── 2_Hammers
│   └── 2.2_Hammer-Tuning/
│       └── Hammer-Tuning-Exercise-2_Accuracy.md (you are here)
└── Data/
    ├── rep10.hammers.tbl (from previous exercise)
    └── rep10.seqs.tbl


## Exercises:

1. Let's start by understanding the relationship between the number of SOURs used and hammer accuracy. Ask Grimoire to explain how increasing the number of SOURs in the hammer creation process can potentially improve the accuracy of genome identification. Also, ask about any potential drawbacks or limitations of using a large number of SOURs.

2. Not every SOUR is the same. Some are easier to find and are more common in the genome, while others are harder to find and are less common. Ask Grimoire to explain what the more useful SOURs are and why they are more useful.

3. TODO: Add more information here that talks about the different algorithms for choosing the SOURs without actually needing to do the process. Since I cannot find the data for all the sequences of the sours.

