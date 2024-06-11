# TSV Exercise 1 - Learning to use Grimoire

Objective: Become familiar with using ChatGPT, and learn about basic bioinformatics, common data-formats, and use of command-line tools.
 
ChatGPT is an AI model that can instruct you on the basics of many topics without needing to look up the information online. ChatGPT can also help you to write computer-programs without first needing to learn a computer-language yourself, and it can often explain to you what went wrong with a program and how to fix it. ChatGPT comes in many specialized versions, but the version that we will be using in these exercises is a "Code-Wizard" called ["Grimoire"](https://chat.openai.com/g/g-n7Rs0IK86-grimoire/).


ChatGPT applications such as Grimoire are very powerful and convenient --- however, do please note that sometimes it will "make things up" because it has been trained to be "helpful" even when it doesn't actually know the answer! (Such "made-up answers are often called "hallucinations".) So you should be cautious about assuming that what ChatGPT variants tell you is always "100% accurate", and "fact-check" claims that it makes that seem implausible to you, particularly with regard to non-programming-related issues.

Thus, "Trust, but verify"! :-)
One of the skills that you will learn during this course will be how to test and debug the software that Grimoire generates for you. We will also provide solutions in each exercise module that you check your results against. And if you find yourself completely wedged, please feel free to contact the course-developers via the Discord Channel.

The following exercises are intended to get you used to interacting with "Grimoire", while also introducing some basic types and formats of bioinformatic data.

## Materials: 

"Grimoire", at <https://chat.openai.com/g/g-n7Rs0IK86-grimoire>

Below is a visual schematic representation of the "Directory tree" for this section of the course. (See the first exercise below if you are unfamiliar with the concept of a "Directory Tree".) Indentation is used to represent subdirectory-levels, and items that end in `/` mean "directory that has this name". In addition to the file you are currently reading, you will also need the files `bindict.tbl` and `tsv_headers.py`,
which are respectively in the `Data` and `Code` directories underneath the directory
`1_Representative-Genomes/`:

<!--
* FIG-Bioinformatics-Course/
    * 1_Representative-Genomes/
        * 1.1_Tab-Separated-Value_(TSV)_Files/
            * TSV-Exercise-1_Learning-to-use-Grimoire.md

        * Data/
            * bindict.tbl
        * Code/
            * tsv_headers.py
-->
```
FIG-Bioinformatics-Course/
└── 1_Representative-Genomes/
    ├── 1.1_Tab-Separated-Value_(TSV)_Files/
    │   └── TSV-Exercise-1_Learning-to-use-Grimoire.md
    ├── Data/
    │   └── bindict.tbl
    └── Code/
        └── tsv_headers.py
```

## Exercises: 

1. Ask Grimoire to explain the concepts of "Directories", "Directory Trees", and "File Paths" to you. In particular, ask it to explain the concept of the directories "dot" and "double-dot" if it did not already do so.

2. Ask Grimoire to explain what a "tab-separated value (TSV) file" with a header-line is.

3. Ask Grimoire to explain what `STDIN` (Standard Input), `STDOUT` (Standard Output), and `STDERR` (Standard Error) mean within the context of a command-line tool. If you are unfamiliar with the concept of a "command line", have Grimoire explain that as well.

4. Ask Grimoire to write a program that will list the names in a TSV-file's header-line columns, and then have it explain to you how the program works "line-by-line".
    * Note: Make sure that you explicitly use the term "line-by-line", as Grimoire may not give a detailed explanation of everything within the program without it.

    * Similarly, if Grimoire's explanation of a particular line of code
    still confuses you, ask it to break down that line "step-by-step".

The key-phrases "line-by-line" and "step-by-step" trigger a particular mode of reasoning that "Large Language Models" (LLMs) have been extensively trained for; in "line-by-line" or "step-by-step" mode, an LLM is likely to reason more clearly, and is less likely to jump to conclusions, make mistakes, or "hallucinate".

5. Use Grimoire's "clipboard" icon at the upper-right of its code-window to copy the program to your clipboard. To the left of VSCode, there are icons in a vertical stack. Click on the top icon that looks like a stack of two pieces of paper that is labelled "File Explorer". Navigate to `1_Representative-Genomes/` within the course directory, then to "1.1 - Tab-Separated Files", and finally to the directory named `Code/`. Within this `Code/` directory, you will see a file named `tsv_headers.py` that will be empty when you open it. Paste Grimoire's code that you saved to your clipboard into the file `tsv_headers.py`, then click on "Save" under the "File" menu to save your script to disk. 

6. Click on "New Terminal" under the VScode "Terminal" menu to open a terminal-window within VScode, and then run `tsv_headers.py` on the file `1_Representative-Genomes/Data/bindict.tbl`. (Grimoire should have shown you how to run the program when it created it, so if you are uncertain how to run your program, refer back to your Grimoire-session, and ask questions if there is something you don't feel you understand yet.) Under `macOS` or `LINUX`, the runtime syntax should be:

``` python3 tsv_headers.py ../Data/bindict.tbl ```

NOTE: The above should be entered as a single line, even if your browser might have wrapped it onto two lines.

Under `Windows`, you should instead replace the 'slashes' with 'backslashes', like this:

``` python3 Code\tsv_headers.py ..\Data\bindict.tbl ```

## Solution Check instructions:
If you are successful, the program should return output that matches the 5 columns in the data file.

```genome_id	genome_name	RepGen.200	RepGen.100	RepGen.50```
