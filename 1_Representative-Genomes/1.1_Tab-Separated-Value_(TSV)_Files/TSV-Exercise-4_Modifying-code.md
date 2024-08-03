# TSV File Exercise 4 - Modifying Code

Objective: Revise an existing program to incorporate new use-cases and inputs

Editing your tools/programs to support various different use-cases makes your own code more versatile and effective. It also helps you to improve your own work and to apply what you have learned to your own work. This exercise focuses on revising your program from TSV Exercise 1.1 to incorporate some of the concepts you have learned regarding use of the command-line interface, as well as those for prompting Grimoire (ChatGPT). This exercise is more "hands off" in that you will mostly be applying a "Review" process to your own code.

## Materials: 
[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

```
FIG-Bioinformatics-Course/
├── 1_Representative-Genomes/
│   └── 1.1_Tab-Separated-Value_\(TSV\)_Files/
│       ├── TSV-Exercise-4_Modifying-code.md   (you are here)
│       └── Solutions/
│           └── tsv_headers_revised_solution.py
├── Code/
│   └── tsv_headers.py
└── Data/
    ├── bindict.tbl
    ├── data.tbl
    └── rep200.list.tbl
```

## Exercises:

1. Ask Grimoire to tell you how `tsv_headers.py` works by "attaching" the file to your question regarding how the code works before submitting it. (You "attach" a file to a prompt by saying _"I am going to attach a program that I would like to work on with you",_ then click on the "paperclip" icon at the left of the "Messages" box (AKA the "prompt"), select the program to upload, and finally click the "enter" icon at the right of the prompt.) If attaching the file does not work, try copying and pasting the contents of the file into the "Messages" box. Then, ask Grimoire to translate the attached code into "pseudocode", and explain the program to you "line-by-line". Insert the pseudocode and Grimoire's explanation into the beginning of `tsv_headers.py` as a "block-comment" for future reference.
* NOTE: "single-line comments" consist of anything that follows a `#` character through the end of the current line.
* Multiline "block-comments" are most easily constructed by placing three single-quotes or double-quotes by themselves on the lines preceding and following the text of the "block-comment", like this:
```
        """
        Comment-line 1
        Comment-line 2
        [...]
        Comment-line N
        """
```

2. In your previous exercises, your program used what are known as "positional arguments", which work like this:

    ``` python3 program_name arg1 arg2 arg3 ... ```

    Python programs also support "named arguments", which work like this:

    ``` python3 program_name -a argA -b argB -c argC ... ```

    Named arguments have the advantage over positional arguments that they can be either optional or mandatory, and unlike positional arguments, their order doesn't matter.

    Named arguments can have a "long form" as well as a short form; "long form" arguments look like this:
    
    ``` python3 program_name --nameA argA --nameB argB --nameC argC ... ```

    Ask Grimoire to tell you more about "short form" and "long form" named arguments, and ask it to give you some examples; then ask it any questions that you might have about named arguments.


3. In this exercise you are going to revise the `tsv_headers.py` program to add named arguments and a new use-case. Use what you learned about Command Line Arguments and Grimoire prompts to make the following improvements to your code. 
    * The program should accept its input datafile-name via the command-line argument `-i` (short for "Input"), e.g. `-i datafilename`.

    * The program should extract the header-names from the first line of the input datafile.

    * The program should accept an optional argument `-n` to specify the "total number of selected data-columns" that will be printed to STDOUT, e.g. `-n 4` means "Print a total of 4 selected columns". If this optional argument is not specified, then the program should print all of the selected columns in the input datafile.

    * The program can take an optional `-m` argument defining a "skip-factor" or "step-factor" by which it will skip over the unselected data-columns, e.g., `-m 10` means "print every 10th column", so that columns [1, 11, 21, 31...] will be selected, up to `-n` total columns.
    
4. Save the modified program as in previous exercises.

5. Once you have finished with these revisions, the program should be able to take the following prompt from the terminal:
    
    ``` python3 Code/tsv_headers.py -i Data/data.tbl -n 3 ```
    
    and it should return the output:

    ``` sample  1033731.3   1034345.3 ```

## NOTES:

1. The file `data.tbl` has over 2000 columns, and so can take a long time to load in applications. Bonus points if you can perform the entire program-revision without opening up the program file to edit it yourself.

2. Grimoire may choose to implement a "vertical" rather than a "horizontal" data-layout. For the purposes of this exercise, we are not worried about the layout.

3. Sometimes Grimoire will misinterpret the `-n` argument to mean the maximum column-index to be printed rather than the total number of columns


## Solution Check instructions:

If you are successful at revising your program, you should see the following output from each of the following commands. (NOTE: once again, commands should be entered as a single line, even if they appear to wrap over multiple lines of the screen):

* ``` python3 Code/tsv_headers.py -i Data/data.tbl -n 7 ```

    ``` sample	1033731.3	1034345.3	1042156.4	1105031.3	1118060.3	1121370.3 ```

* ``` python3 Code/tsv_headers.py -i Data/data.tbl -n 4 -m 20 ```

    ``` sample 1507.3 203120.7 40545.1270 ```

* ``` python3 Code/tsv_headers.py -i Data/rep200.list.tbl -n 5 -m 2 ```

    ``` genome_id   domain  species score ```

    NOTE: only 4 columns are printed, because the file `Data/rep200.list.tbl`
    only has 8 columns, so the script ran out of data before it had a chance to emit 5 columns

* ``` python3 Code/tsv_headers.py -i Data/bindict.tbl ```

    ``` genome_id	genome_name	RepGen.200	RepGen.100	RepGen.50 ```

We have included our own Grimoire-generated solution-code for these exercises as
`1.1_Tab-Separated-Value_\(TSV\)_Files/Solutions/tsv_headers_revised_solution.py`.
Our solution may differ in detail from the one that Grimoire found for you,
since Grimoire sometimes makes different choices regarding variable-names and programming-style, but your code and our code should be __functionally-equivalent__.
