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
├── Code/
│   └── find_nearest_reference.py 
└── Data/
    └── rep10.seqs.seed_protein.faa
```

## Overview:

In research and medicine, we will often be confronted with a genetic sample that contains one or more unidentified genomes.
We can gain knowledge about which genomes are present in a sample by comparing the sample to a "Set of Representative Genomes" (RepGen Set).
The set `rep200` is our highest-resolution "production" RepGenSet; it allows identification of genomes down to the "species" level. Similarly, `rep100` allows identification down to the "Genus" level.
Our coarsest set of representatives is `rep10`, which allows identification to within broad families of genomes.

In this exercise, you will first install the BV-BRC Command-Line Toolkit (also referred to as the "P3 Commands" for historical reasons), which will allow you to fetch the data and metadata for nearly a million sequenced genomes.
You will used the CLI to fetch the PhlseS sequence for a "Mystery Genome" and then compare it to the `rep200` data using the tool `find_nearest_reference.py` that you constructed in Kmer-Ex-3. You will then fetch the actual identity of this genome from BV-BRC and compare it to the nearest representative that your tool selected.

## Exercises:

1. In order to access the data at BV-BRC, you will first need to register. Please go to
https://www.bv-brc.org/docs/quick_references/registration.html,
then click on the appropriate link within the left-sidebar or scroll down to the section for your operating-system (i.e., macOS, LINUX, or Windows) and follow the installation instructions within that section.

2. Once you have registered, you will need to download and install the Command-Line Interface toolkit as documented here:<br>
https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html#cli-installation

3. Once installed, under `macOS` and `Windows` launch the `BV-BRC.app` by double-clicking on the app icon. (Under LINUX, the app commands will be installed on your default `PATH`, so just open a shell window and start typing commands)

4. Login to BV-BRC using this command:
```
    p3-login yourBVBRCusername
```
which will prompt you for your BV-BRC password;
once logged in, you can issue any P3-command.

5. Get the PheS sequence for "Mystery Genome 1491.662" using the following incantation:
```
    p3-echo PheS | p3-find-features --attr patric_id,product --eq genome_id,1491.662 gene | p3-get-feature-sequence --col feature.patric_id > Data/mystery_PheS.faa
```
(Again, remember that the above command should be pasted in as a single line, even though it may appear to wrap around onto several lines on the screen.) We will break down this complex command later in this exercise.

6. Use your command `find_nearest_reference.py` to compare the sequence you just fetched to `rep10`:
```
    pyhthon3 Code/find_nearest_reference.py -K 8 -R Data/rep10.seqs.seed_protein.faa < Data/mystery_PheS.faa
```
Which representative genome was reported in the RepSet-description column?

7. Now let's find out the true identity of the "Mystery Genome":
```
    p3-echo 1491.662 | p3-get-genome-data --attr genome_name
```
Does the "Mystery Genome" have the same genus and species (first and second name) as the RepGen genome that `find_nearest_reference.py` found?<br>
Congratulations! you have correctly identified the genus and species of the "Mystery Genome"!

The above is a greatly simplified "cartoon version" of how one may go about identifying a completely new genome using a RepGen Set. The actual real-world procedure is more involved; one would first have to sequence and "assemble" the mystery-genome, annotate its genes, extract the sequence for the "Mystery PheS", and then perform the sequence comparison, but the basic concept is the same.

In a later exercise, you will learn about a type of "Genetic Barcode" called a "Hammer", which allows one to skip the intermediate steps to go directly from raw sequencer-data to a genome identification.

## Breaking down the steps in the "P3 incantations"
