# RepGen Exercise 2 - Finding the nearest RepGen

Objective: Fetch the PheS for a "Mystery Genome" from the [Bacterial and Viral Bioinformatic Resource-Center (BV-BRC)](https://www.bv-brc.org/), and find out which Representative Genome it is closest to.

## Materials

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

[The BV-BRC Toolkit](https://www.bv-brc.org/docs/cli_tutorial/index.html)

```
FIG-Bioinformatics-Course/
├── 1_Representative-Genomes/
│   └── 1.4_Building-Representative-Sets/
│       └── RepGen-Exercise-2_Finding-the-nearest-RepGen.md (you are here)
└── Data/
    └── rep10.seqs.seed_protein.faa
```

## Overview:

In research and medicine, we will often be confronted with a genetic sample that contains one or more unidentified genomes.
We can gain knowledge about which genomes are present in a sample by comparing the sample to a "Set of Representative Genomes" (RepGen Set).
The set `rep200` is our highest-resolution "production" RepGenSet; it allows identification of genomes down to the "species" level. Similarly, `rep100` allows identification down to the "Genus" level.
Our coarsest set of representatives is `rep10`, which allows identification to within broad families of genomes.

In this exercise, you will first install the BV-BRC Command-Line Toolkit, which will allow you to fetch the data and metadata for nearly a million sequenced genomes.
You will used the CLI to fetch the PHeS sequence for a "Mystery Genome" and then compare it to the `rep200` data using the tool `find_nearest_reference.py` that you constructed in Kmer-Ex-3. You will then fetch the actual identity of this genome from BV-BRC and compare it to the nearest representative that your tool selected.

## Exercises:

1. In order to access the data at BV-BRC, you will first need to register. Please go to
https://www.bv-brc.org/docs/quick_references/registration.html
and follow the instructions.

2. Once you have registered, you will need to download and install the Command-Line Interface toolkit as documented here:<br>
https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html#cli-installation

3. 