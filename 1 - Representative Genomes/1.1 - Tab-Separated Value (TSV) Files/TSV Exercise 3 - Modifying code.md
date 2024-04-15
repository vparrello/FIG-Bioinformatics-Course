#### TSV File Exercise 3

 Objective: Adjust former programs for new use cases and inputs

 Editing your tools/programs into different use cases makes your own code more versitile and efficient. It also helps you to edit your own work and help you apply what you have learned to your own work. This exercise focuses on adjusting your first program to apply some of the concepts of the command line as well as those for prompting Grimoire (Chatgpt). This exercise is more hands off to help you apply the review process to your own code.

#### Materials: 
[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)
```
    FIG-Bioinformatics-Course/
        1 - Representative Genomes/
            Data/bindict.tbl
            Data/data.tbl
            Data/rep200.list.tbl
            bin/tsv_reader.py
```

#### Exercise: 

1. Ask Grimoire to tell you how `tsv_reader.py` works by attaching it as a file. (You "attach" a file to a prompt by saying "I am going to attach a program that I would like to work on with you", and then click on the "paperclip" icon at the left of the "Messages" box (AKA the "prompt"), select the program to upload, and then click the "enter" icon at the right of the prompt.) If that does not work, try copying and pasting the contents of the file into the "Messages" box. Then, ask Grimoire to explain the program to you "line-by-line". Insert the most important topics into the beginning of the file as a "comment" for future reference.\
NOTE: "single-line comments" consist of anything that follows a "#" character.\
Multiline "block-comment" are most easily constructed by placing three single-quotes or double-quotes by themselves on the lines preceding and following the text of the "block-comment", like this:
```
        """
        Comment-line 1
        Comment-line 2
        [...]
        Comment-line N
        """
```

2. In you previous exercises, your program used what are known as "positional arguments", which work like this:

    ``` python3 program_name arg1 arg2 arg3 ... ```

    Python programs also support "named arguments", which work like this:

    ``` python3 program_name -a argA -b argB -c argC ... ```

    Named arguments have the advantage over positional arguments that they can be either optional or mandatory, and unlike positional arguments, their order doesn't matter.

    Named arguments can have a "long form" as well as a short form; "long form" arguments look like this:
    
    ``` python3 program_name --nameA argA --nameB argB --nameC argC ... ```

    Ask Grimoire to tell you more about "short form" and "long form" named arguments, and ask it to give you some examples; then ask it any questions that you might have about named arguments.


3. In this exercise you are going to adjust the 'tsv_reader.py' program to add named arguments and a new use case. Use what you learned about Command Line Arguments and Grimoire prompts to make the following improvements on your code. 
    * The program should read the input-filename from command-line argument ```-i```
    * The program should extract the header-names from the input-file
    * The program can take an optional argument ```-n``` to specify how many of the first columns will be printed to standard output. If this argument is not specified, then the program prints all columns
    * The program can take an optional ```-m``` argument defining a "skip-factor" by which it will skip over the columns. (e.g., print every 2nd, 4th, or 8th column, for example).
    
3. Once you have finished, it should be able to take this prompt from the terminal:
    
    ``` python3 tsv_reader -i ../Data/data.tbl -n 8 -m 4 ```
    
    And come back with
    ``` 1042156.4 1121445.4 ```
    This data file has over 2000 columns and can take a long time to load in applications. Bonus points if you can do this entire program without opening the file yourself.

## Solution Check instructions:
If you are successful, you will have the following output for the related commands

``` python3 tsv_reader -i ../Data/data.tbl -n 7 ```

``` sample	1033731.3	1034345.3	1042156.4	1105031.3	1118060.3	1121370.3 ```

``` python3 tsv_reader -i ../Data/data.tbl -n 4 -m 20 ```

``` 1496.3893 203120.7 40545.1270 563192.3 ```

``` python3 tsv_reader -i ../Data/rep200.list.tbl -n 5 -m 2 ```

``` genome_name genus rep_id distance ```

``` python3 tsv_reader -i ../Data/bindict.tbl ```

``` genome_id	genome_name	RepGen.200	RepGen.100	RepGen.50 ```
