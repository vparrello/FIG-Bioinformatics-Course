# Hammer Worthiness Exercise 3

Objective: Understand the concept of hammer worthiness by auditing hammers against representative and peer genome sequences.

When evaluating the worthiness of hammers, we're assessing how well they do the job of representing both the representative genome and its peer genomes inside its representative space. A hammer is considered 100% worthy when it appears in the representative genome and in 100% of the peer genomes it's meant to represent. This can be useful for creating antibiotics that wish to target a specific bacterium, for example. However, a 1% worthy hammer set can also be useful as this will create a more diverse set of genomes that can be represented by that representative genome. 

This exercise will help you understand how to assess hammer worthiness and how it can be used to create more accurate and diverse genomic representations.

## Materials:
[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

```
FIG-Bioinformatics-Course/
├── 2_Hammers
│   └── 2.2_Hammer-Tuning/
│       └── Hammer-Tuning-Exercise-3-Worthiness.md (you are here)
└── Data/
    ├── rep10.seqs.tbl (containing representative genome sequences)
    └── rep10.hammers.tbl (list of hammers for each representative genome)
```

## Exercises:

1. Start by asking Grimoire to explain the concept of hammer worthiness. What makes a hammer 100% worthy, and why is this important in genomic analysis?

2. Ask Grimoire to write you a script in Python that will help you calculate the worthiness of the hammers. It needs to incorporate the following requirements in the file `hammer_worthiness.py`:
    - Read the `rep10.hammers.tbl` file.
    - It should take as input the genome ID of the representative genome you are interested in
    - It should read a couple of peer genome files for that representative genome from the Data folder.
    - For that representative genome, create an output file with the name `{repgen_name}.hammerworthiness.tbl`
    - The output file should contain the hammers for that representative genome and the percentage of peer genomes that contain each hammer.

    **Note: The peer genome file that is reference in this prompt is yet to be created. You will have to create it using the `p3-genome-fasta` command. which we will get into shortly.**

3. The next thing we need to do is to find the peer genomes for the representative genome that you are interested in making the hammers worthy for. Optimally, we would do this exercise for all the representative genomes in the dataset. However, for this exercise, we will only focus on one representative genome in the dataset to learn the process. If you are not sure which genomes are peers of the representative genome, look back at the Representative-Genomes unit as a refresher.

4. Once you have a list of peer genomes, you can use the p3-genome-fasta command to download their contigs. The basic syntax is:
   ```
   p3-genome-fasta --contig [genome_id1,genome_id2,...] --output peer_contigs.fasta
   ```

5. Once you have the peer contigs, you can feed them into your hammer worthiness script to get the worthiness of the hammers for that representative genome. Keep in mind that the more peer genomes you have, the more accurate the worthiness assessment will be.

6. Once you have the worthiness of the hammers for that representative genome, you can follow the same steps for the other representative genomes in the dataset to get the worthiness of the hammers for all the representative genomes. This will require changing your `hammer_worthiness.py` script to loop over all the representative genomes in the dataset and download the peer genomes for each representative genome. Ask Grimoire to help you as this is a major change and can require multiple steps.

7. Keep in mind that the worthiness of the hammers for a representative genome does not always have to be 100%. Sometimes, you may want to create a more diverse set of genomes that can be represented by the representative genome. This can be useful for creating antibiotics that can target a specific bacterium, for example. Once you create the worthy hammers for 100%, test them on a sample that you have used before in previous exercises. See how different the results are from the previous exercise.

