# TSV Exercise 1 - Learning to use Grimoire

Objective: Become familiar with using ChatGPT, and learn about basic bioinformatics, common data-formats, and use of command-line tools.
 
ChatGPT is an AI model that can instruct you on the basics of many topics without needing to look up the information online. ChatGPT can also help you to write computer-programs without first needing to learn a computer-language yourself, and it can often explain to you what went wrong with a program and how to fix it. ChatGPT comes in many specialized versions, but the version that we will be using in these exercises is a "Code-Wizard" called ["Grimoire"](https://chat.openai.com/g/g-n7Rs0IK86-grimoire/).


ChatGPT applications such as Grimoire are very powerful and convenient --- however, do please note that sometimes it will "make things up" because it has been trained to be "helpful" even when it doesn't actually know the answer! (Such "made-up answers are often called "hallucinations".) So you should be cautious about assuming that what ChatGPT variants tell you is always "100% accurate", and "fact-check" claims that it makes that seem implausible to you, particularly with regard to non-programming-related issues.

Thus, "Trust, but verify"! :-)
One of the skills that you will learn during this course will be how to test and debug the software that Grimoire generates for you. We will also provide solutions in each exercise module that you check your results against. And if you find yourself completely wedged, please feel free to contact the course-developers via the Discord Channel.

The following exercises are intended to get you used to interacting with "Grimoire", while also introducing some basic types and formats of bioinformatic data.

## Materials: 

"Grimoire", at <https://chat.openai.com/g/g-n7Rs0IK86-grimoire>

Below is a visual schematic representation of the "Directory tree" for this section of the course. (See the first exercise below if you are unfamiliar with the concept of a "Directory Tree".) Indentation is used to represent subdirectory-levels, and items that end in `/` mean "directory that has this name"; otherwise, the item is a program or data file. In addition to the file you are currently reading, you will also need the files `bindict.tbl` and `tsv_headers.py`,
which are respectively in the `Data/` and `Code/` directories underneath the main directory `FIG-Bioinformatics-Course/`:

```
FIG-Bioinformatics-Course/
├── Code/
│   └── tsv_headers.py
├── Data/
│   └── bindict.tbl
└── 1_Representative-Genomes/
    └── 1.1_Tab-Separated-Value_(TSV)_Files/
        ├── TSV-Exercise-1_Learning-to-use-Grimoire.md  (you are here)
        └── Solutions/
            └── tsv_headers_solution.py
```

## Exercises: 

1. Ask Grimoire to explain the concepts of "Directories", "Directory Trees", and "File Paths" to you. In particular, ask it to explain the concept of the directories "dot" and "double-dot" if it did not already do so.

2. It is conventional to end a filename with a "extension" that indicates to the user what kind of file it is and what programs can open it. Ask Grimoire to explain the concept of a "file extension" to you.

NOTE: In this course, most files will have extensions of either `.py` indicating that they are Python-language programs, or `.tbl` indicating that they are a "table" of some type, although we will introduce some other special-purpose file-extensions later in the course.

3. Ask Grimoire to explain what a "tab-separated value (TSV) file" with a header-line is.

4. Ask Grimoire to explain what the terms `STDIN` (Standard Input), `STDOUT` (Standard Output), and `STDERR` (Standard Error) mean within the context of a command-line tool. If you are unfamiliar with the concept of a "command line", have Grimoire explain that as well.

5. Ask Grimoire to write a Python program that will read a tab-separated-value data-file with header-line from `STDIN`, print the field-names in that TSV-file's header-line columns to `STDOUT`, and then exit; any error-messages should be printed to `STDERR`. Then ask Grimoire to explain to you how this program works "line-by-line".
    * Note: Make sure that you explicitly use the term "line-by-line", as Grimoire may not give a detailed explanation of everything within the program without it.

    * Similarly, if Grimoire's explanation of a particular line of code
    still confuses you, ask it to break down that line "step-by-step".

The key-phrases "line-by-line" and "step-by-step" trigger a particular mode of reasoning that "Large Language Models" (LLMs) have been extensively trained for; in "line-by-line" or "step-by-step" mode, an LLM is likely to reason more clearly, and is less likely to jump to conclusions, make mistakes, or "hallucinate".

6. Use Grimoire's "clipboard" icon at the upper-right of its code-window to copy the program to your clipboard. Launch VScode, and click on "Open Folder" under the "File" menu, which opens the "File Explorer". Select the folder `FIG-Bioinformatics-Course/` then click `Open`. Within this folder near the bottom of the file-explorer window, you will find a directory named `Code/`. Click on `Code/` to expand this directory-listing. Within this `Code/` directory, you will see a number of files; scroll down to the file named `tsv_headers.py` and click on it, which will open that file in the file-editor. You will see that  `tsv_headers.py` is empty except for the following "comment":
```
# Paste Grimoire's code below this line,
# then select and click "Save" under the "File" menu
```
Paste Grimoire's code that you saved to your clipboard into the file `tsv_headers.py` below the comment, then click on "Save" under the "File" menu as instructed, which saves your program to disk. 

7. You are now in a position to run your program.
Grimoire has probably already shown you an example of how to run the program during its summary-discussion of the program, but since Grimoire is ignorant of the details of your operating-system, directory-structure, and the fact that you are using the VScode development-environment (unless it has explicitly been told these things), its example probably will not work verbatim. So, to run the program in your environment, please follow the steps below: 
* Go back to the "File-Explorer" window and select the top-level directory `FIG-Bioinformatics-Course/`.'

* Click on "New Terminal" under the VScode "Terminal" menu to open a terminal-window within VScode. Click on the "Terminal" window that VScode just opened, which will switch VScode's "focus" to that window. You can confirm that you are in the correct directory by entering `pwd` (short for "print working-directory") and hitting the `return` key; your computer should respond with a directory-path ending in `FIG-Bioinformatics-Course`.

* To run `tsv_headers.py` on the file `Data/bindict.tbl` under `macOS` or `LINUX`, the runtime syntax should be something like this:

```
    python3 Code/tsv_headers.py < Data/bindict.tbl
```

This command asks the `python3` program-interpreter to run the program `Code/tsv_headers.py`, reading its data from the file `Data/bindict.tbl`.

```
NOTE: The above command should be entered as a single line, even if your browser might have wrapped it onto two lines.
```

* Under `Windows`, you should instead replace the 'slashes' with 'backslashes', like this:

```
    python3 Code\tsv_headers.py < Data\bindict.tbl
```

## Solution Check instructions:
If you are successful, the program should return output that matches the 5 columns in the data file.

```genome_id	genome_name	RepGen.200	RepGen.100	RepGen.50```

In the `Solutions/` subdirectory of this module, we have included the program `tsv_headers_solution.py` that Grimoire wrote for us; this program may perhaps differ in precise detail from the program that Grimoire wrote for you, but it should be functionally equivalent.
