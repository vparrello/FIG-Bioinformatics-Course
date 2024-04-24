#### FASTA Exercise 1 - What is FASTA Format?

Objective: Become familiar with the FASTA file-format. 

FASTA is one of the most basic file-formats for storage and exchange of biological sequences. We use them constantly when looking at data in the bioinformatics space. In order to use them, we first need to understand how they came to be and what information is available in them. A sample FASTA file has been put into the Data folder for you to explore. 

#### Materials: 

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

* FIG-Bioinformatics-Course/
    * 1_Representative-Genomes/
        * Data/
            * Sample1.fasta
        * 1.2_FASTA-Files/
            * FASTA-Exercise-1_What-is-FASTA-format.md


#### Exercise: 

The FASTA file-format is defined as follows:

* A FASTA-file consists of one or more "records".

* A record has the following format:
    - Each record must begin with "record header" whose first character is a "greater-than" sign `>` and is immediately followed by a sequence-identifer. Sequence-identifiers must not contain spaces.

    - Any text on the record-header following the first space-character through the end of the line is considered to be a "sequence description" or "comment". "Sequence descriptions" are optional.
    "Sequence descriptions" are often used to describe the function of a sequence, or to indicate which genome the sequence came from.

    - The lines following the "record header" contain sequence-data.
    Sequence-data may be entered all on one line, or may be broken into multiple lines; however the sequence-data must not contain blank lines.

1. Ask Grimoire to explain the three important types of biological sequences: DNA, RNA, and protein.
(If Grimoire uses any words or terms that you are not familiar with, please ask it to explain these terms.)

2. Ask Grimoire to explain the "alphabets" that are used to represent DNA, RNA, and protein sequences.

4. Ask Grimoire to tell you the history of the FASTA format. 

5. Ask Grimoire to list the preferred file-extensions conventionally used for FASTA-format files, and explain what types of sequence-data each preferred file-extension is used to represent and store.

6(Suggestion). Pick one of the sequence-data types that appeals to you. Ask Grimoire to show you what data inside a file of that type would look like. Then pick another and do the same. Notice the differences between the data. 
    Ask Grimoire what situation or use case would need both types of data

6(Suggestion). Pick one of the sequence-data types that appeals to you.  Ask Grimoire what a scientist could learn from that type of data. 

6(Suggestion). Come up with some sample sequence of data that you think fits one of the extension types. Send Grimoire your own sample data and ask it to fix it according to each of the file extensions and how it would look different. 

## Solution Check instructions:


# NOTES:
    1st and 2nd suggestion to be tested in Grimoire
    Put bullets for all file structures
    Add file extension why the file extensions matter to the user and not code 
    DONE: How the two alphabets are different between protein and DNA
