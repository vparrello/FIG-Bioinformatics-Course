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

In this exercise, you will first install the `BV-BRC Command-Line Interface` Toolkit (`BV-BRC CLI`, also referred to as the "P3 Commands" for historical reasons), which will allow you to fetch the data and metadata for nearly a million sequenced genomes.
You will use the CLI to fetch the PheS sequence from a "Mystery Genome" and then compare it to the `rep200` data using the tool `find_nearest_reference.py` that you constructed in Kmer-Exercise-3. You will then fetch the actual biological name of the "Mystery Genome" from BV-BRC and compare it to the biological name of nearest representative genomne that your tool selected.

## Exercises:

1. In order to access the data at BV-BRC, you will first need to register. Please go to
https://www.bv-brc.org/docs/quick_references/registration.html,
and follow the istructions.

2. Once you have registered, you will need to download and install the Command-Line Interface toolkit as documented here:<br>
https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html#cli-installation<br>
BV-BRC supports three operating-systems; detailed installation-instructions for each OS may be found here:
    * [macOS](https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html#installation-on-macos)
    * [LINUX](https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html#installation-on-debian-ubuntu-mint-linux)
    * [Windows](https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html#installation-on-windows)

3. Once installed, launch the `BV-BRC.app` under `macOS` and `Windows` by double-clicking on the app's icon. It should look very similar to a terminal icon. (Under LINUX, the app commands will be installed on your default `PATH`, so just open a new terminal window and start typing commands.)

4. Login to BV-BRC using this command:
```
    p3-login your_BVBRC_username
```
which will prompt you for your BV-BRC password;
once logged in, you can issue any P3-command.

5. Get the PheS sequence for "Mystery Genome-ID 1491.662" using the following incantation:
```
    p3-echo PheS | p3-find-features --attr patric_id,product --eq genome_id,1491.662 gene | p3-get-feature-sequence --col feature.patric_id > Data/mystery_PheS.faa
```
(Again, remember that the above command should be pasted in as a single line, even though it may appear to wrap around onto several lines on the screen.) We will break down this complex command later in this exercise.

6. Use your command `find_nearest_reference.py` to compare the sequence you just fetched to `rep10`:
```
    python3 Code/find_nearest_reference.py -K 8 -R Data/rep10.seqs.seed_protein.faa < Data/mystery_PheS.faa
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

Earlier, we used the following "incantation":
```
p3-echo PheS | p3-find-features --attr patric_id,product --eq genome_id,1491.662 gene | p3-get-feature-sequence --col feature.patric_id > Data/mystery_PheS.faa
```
The above complex command is an example of a "command pipeline": A sequence of commands with modifying arguments separated by the symbol `|`, which is pronounced as "pipe" because it connects (or "pipes") the `STDOUT` from one command to the `STDIN` of the following command.
You can use a "pipeline" to construct an arbitrarily complex command or filter operation out of simple commands.

Let's break tis incantation down step-by-step to see what it does.

The first step is `p3-echo PheS`. If you type this command by itself, you will see the following output:
```
id
PheS
```
The above illustrates two points:
* P3-commands all generate "Tab-separated output with header-line", which is sent to STDOUT 

* the `p3-echo` command accepts one (or more) values as arguments,
and "echoes" them one after the other following the column-name.
The utility of `p3-echo` is to create TSV data that will be sent
to subsequent P3-commands

In the absence of an argument, `p3-echo` uses a "default" column-name of `id`. In cases where a different column-name is required, it can be specified using the `-t` argument (which you may think of as short for "column-title"), so for example:
```
p3-echo -t feature_id 'fig|1491.662.peg.1715'
```
will generate:
```
feature_id
fig|1491.662.peg.1715
```
It is possible (but cumbersome) to generate multiple-column TSV files using `p3-echo`, but this capability is rarely needed. Most often, you will use `p3-echo` to create single-column output as follows:
```
p3-echo -t columnName value1 value2 value3
```
which results in:
```
columnName
value1
value2
value3
```

The next step is to "pipe" the output of `p3-echo` to the input of `p3-find-features` as follows:
```
p3-echo PheS | p3-find-features --attr patric_id,product --eq genome_id,1491.662 gene
```
which produces:
```
id	feature.patric_id	feature.product
PheS	fig|1491.662.peg.1715	Phenylalanyl-tRNA synthetase alpha chain (EC 6.1.1.20)
```
Note that the output has acquire two new columns named `feature.patric_id`	and `feature.product`; the portion `feature` before the dot specifies the `type` of the data-column, while the portion after the dot specifies the "attributes" of the data-type that were selected in the command by the argument `--attr patric_id,product`.

This mechanism of appending columns to the end of each line allows us to incrementally build up arbitrarly complex data-tables via a sequence of commands. 

the final step is:
```
p3-echo PheS | p3-find-features --attr patric_id,product --eq genome_id,1491.662 gene | p3-get-feature-sequence --col feature.patric_id
```
which yields:
```
>fig|1491.662.peg.1715 Phenylalanyl-tRNA synthetase alpha chain (EC 6.1.1.20)
MRQKLEDIKNSAINELKTTLSKDQLEAIRVKYLGKKGELTQILRGMGALSQEERPIVGKVANEVRSYIEETIKEAFSDIKNKEKSIRLENETIDITMPGKKQAVGKRHPLDLTLESMKDIFISMGFTIEEGPEVELDKYNFEALNIPKNHPARGEQDTFYINDNLVLRTQTSPIQIRTMENQKPPIKMIAPGKVYRSDSVDATHSPIFYQMEGLVVDKGITFSDLKGTLELFAKRMFGDKVKTKFRPHHFPFTEPSAEMDATCFVCNGEGCKVCKGSGWIELLGCGMVHPQVLRNCNIDPEVYSGFAFGFGVDRMVMMKYGIDDIRLLYESDMRFLNQF
```
The above illustrates a final point, which is that in general,
most P3-commands will use one of the columns as the "key"
that it will use to look up data-values, and the "key" column-name
can be specied using the `--col` (short for "column") command-argument. In this case, we used `--col feature.patric_id` to specify that the
`feature.patric_id` column-value should be used as the data lookup-key
to be used by `p3-get-feature-sequence`.

In the absence of a `--col` argument, P3-commands default to using
the _last_ column of their input-table as their "key", which often simplifies incrementally building up an arbitrarily complex table by using the last attribute looked up as the "key" to look up the next data-item.

Finally, note that every P3-command accepts a `--help` argument
that will list and describe the arguments that the command accepts,
and that most commands also accept a `--fields` argument that lists
the data-attributes that they can access.

A detailed tutorial on how to use P3-commands can be found at:<br>
https://www.bv-brc.org/docs/cli_tutorial/cli_getting_started.html

A catalogue of the P3-commands may be found at:<br>
https://www.bv-brc.org/docs/cli_tutorial/command_list/index.html<br>
Clicking on a command-name will go to the `--help` output for that command.

