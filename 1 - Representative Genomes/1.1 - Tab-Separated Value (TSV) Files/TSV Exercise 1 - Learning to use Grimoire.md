#### TSV Exercise 1 - Learning to use Grimoire

 Objective: Become familiar with using ChatGPT, and learn about basic bioinformatics, common data-formats, and use of command-line tools.
 
ChatGPT is an AI model that can instruct you on the basics of many topics and to write computer-programs, without needing to look up the information or to learn a computer-language yourself. ChatGPT comes in many specialized versions, but the version that we will be using in these exercises is a "Code-Wizard" called ["Grimoire"](https://chat.openai.com/g/g-n7Rs0IK86-grimoire/).

ChatGPT applications such as Grimoire are very powerful and convenient --- however, do please that sometimes it will "make things up" because it has been trained to be "helpful" even when it doesn't actually know the answer So you should be cautious about assuming that what ChatGPT applications tell you is "100%  accurate".
So, "Trust, but verify"! :-)

The following exercises are intended to get you used to interacting with "Grimoire", while also introducings some basic types and formats of bioinformatic data.

#### Materials: 

[chat.openai.com](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)
Data/bindict.tbl
Week 1/tsv_reader.py

#### Exercise: 

1. Ask Grimoire to explain what a "tab-separated value" (TSV) file with a header-line is.

2. Ask Grimoire to explain "STDIN" (Standard Input) and "STDOUT" (Standard Output) in the context of a command- line tool. If you are unfamiliar with the concept of a "command line", have Grimoire explain that as well.

3. Ask Grimoire to write a program that will list the names in a TSV-file's header-line columns, and then have it explain to you how the program works "line-by-line"

    **Note: Make sure that you explicitly use the term "line-by-line", as Grimoire may not give a detailed explaination of everything in the program without it.

4. Use Grimoire's "clipboard" icon at the upper-right of its code-window to copy the program to your clipboard. Launch VScode, click on "Open Folder" under the "File" menu, and navigate down to "FIG-Bioinformatics-Course", then "1 - Representative Genomes" within the course directory, then "1.1 - Tab-Separated Files", and finally a directory named "bin" (the traditional name for a directory that contains executable code  --- it originally stood for "binary"). Within this last directory, you will see a file named "tsv_reader.py", which will be empty when you open it. Paste the code within your clipboard into that file, then click on "Save" under the "File" menu to save your script to disk. 

5. Click on "New Terminal" under the VScode "Termimal" menu to open a terminal-window within VScode, and then run "tsv_reader.py" on the file "1 - Representative Genomes/Data/bindict.tbl". (Grimoire should have shown you how to run the program when it created it, so if you are uncertain how to run your program, refer back to your Grimoire-session, and ask questions if there is something you don't feel you understand yet.)

## Solution Check instructions:
If you are successful, you will have the output that matches the 5 columns in the data file.
"genome_id	genome_name	RepGen.200	RepGen.100	RepGen.50"