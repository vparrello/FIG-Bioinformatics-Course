#### Tab-Separated Value File Exercise 1

 Objective: Become familiar with chatgpt and learn the context of Bioinformatics
 
ChatGPT is an AI model that can help you learn the basics of many topics and write computer-programs without needing to look for the information or learn a computer-language yourself. ChatGPT comes in many specialized versions; the version that we will be using in these exercises is called "Grimoire".
(However it also can make up information because it has been trained to be "helpful", so you should be cautious about always taking what it says as accurate.)

The following exercises are intended to get you used to interacting with "Grimoire", while also introducings some basic types and formats of bioinformatic data.

#### Materials: 

[chat.openai.com](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)
Data/bindict.tbl
Week 1/tsv_reader.py

#### Exercise: 

1. Ask Grimoire to explain what a tab-separated file with a header-line is.

2. Ask Grimoire to explain "STDIN" (standard in) and "STDOUT" (Standard Out) in the context of a command line. If you are unfamiliar with the concept of a command line, have Grimoire explain that as well.

3. Ask Grimoire to write a program that will list the names in a TSV-file's header-line columns, and then explain to you how the program works "line-by-line"

    **Note: Make sure that you explicitly use the term "line-by-line" as Grimoire will not explain everything you need without it.

4. Use Grimoire's "clipboard" icon in its code-window to copy the program to your clipboard.Under the Directory "1 - Representative Genomes/Week 1", you will see a file called "tsv_reader.py" which is empty if you open it.  Then paste the code you have on your clipboard into that file. 

6. Run "tsv_reader.py" on the file "1 - Representative Genomes/Data/bindict.tbl". 

## Solution Check instructions:
If you are successful, you will have the output that matches the 5 columns in the data file.
"genome_id	genome_name	RepGen.200	RepGen.100	RepGen.50"