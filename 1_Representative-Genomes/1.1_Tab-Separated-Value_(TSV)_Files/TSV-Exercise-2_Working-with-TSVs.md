#### TSV Exercise 2 - Working with TSV Files

 Objective: Use Grimoire to create a program that reads Tab Separated Files from the Command Line.

 This exercise focuses on creating a program that reads in a list of column-names as command-line-arguments, extracts from an input TSV-file only those columns whose header-names match the names in the argument-list, and writes those columns to output. (For example, suppose that we have a TSV-file that has 15 columns, but we only care about the data contained in 3 of theose columns.)

 You will be going through the normal development process of a programmer, which is similar to the process of a writer, namely, Outline, Draft, Edit, Review, and then Publish. We will start with a description of the program, express that description as "Pseudo-code", Create the program, Debug for errors, Test output quality, and then Publish. See if you can tell which step you are in as you do the exercise.

#### Materials: 

"Grimoire" at <https://chat.openai.com/g/g-n7Rs0IK86-grimoire>

You will also need the files "rep200.list.tbl",
"command_line_kung.py" and "command_line_fu.py";
indentation is used to represent directory-levels: 

* FIG-Bioinformatics-Course/
    * 1_Representative-Genomes/
    * Data/
        * rep200.list.tbl
    * Code/
        * command_line_kung.py
        * Code/command_line_fu.py

#### Exercise: 

1. Ask Grimoire to explain what "arguments of a command-line program" means.

2. Ask Grimoire to write a program that will accept a list of keywords as its command-line arguments. Then ask Grimoire to explain the code to you "line-by-line" if it did not do so.

3. You will now prompt Grimoire to create a custom Python program to accomplish the following tasks. Create your own prompt to make Grimoire tell you this list of items specifically, then you can accept the code that it creates and paste it into the file "command_line_kung.py". This time also grab the pseudocode Grimoire gives you. (If Grimoire skips the step of giving you pseudocode, you can explicitly ask Grimoire to translate its code into pseudocode, which will be more readable and will help you to better understand what the real code is doing.)

    * The program should be written in python.
    * The program should read a tab-separated file from STDIN with a header-line as the first line.
    * The program should accept a list of keywords as command-line arguments.
    * The program should extract the columns whose headers match the keywords, and write those columns to the standard output in the same order that they were listed.
    * The program should warn the user (via STDERR so as not to pollute STDOUT) if any keyword does not match the header-line, and then exit.
    * The program should print something specific (again via STDERR) to let you know that it has completed it's task.

4. Take your prompt from above and reorder the steps you gave Grimoire. See how Grimoire changes the code and pseudo code. Paste both into the file "command_line_fu.py"

5. Run one of the programs and check to see no errors come up. If there is an error, type the error message into Grimoire and see what solution it can come up with to help you.
    * Note: Grimoire does not know your file system so it might not catch every detail. (For example, under Windows, the file-path separator is "\\" instead of "/".) If you come across a "File Not Found" error, that is because it cannot find the data file. You would need to use the "Absolute Path" of the file. Grimoire can teach you how to find the Absolute Path of the file, and instruct you on what any error-messages that you get mean, and suggest possible fixes. You got this!!

6. In order to run the program, go to the Terminal Menu at the top of VS code and click on New Terminal to open a new terminal underneath your code. Then you want to type in the following to run your command:

    ``` python3 Code/command_line_kung.py genome_id genome_name domain genus species rep_id score distance < Data/rep200.list.tbl ```

    Also try it with the second program but do not worry if this one does not succeed:
    
    ``` python3 Code/command_line_fu.py genome_id genome_name domain genus species rep_id score distance < Data/rep200.list.tbl ```

7. If your command_line_kung program did not succeed, go back to step 3. If it did, then try it again with this command. This one should throw a warning message and then stop the program completely
    
    ```python3 Code/command_line_kung.py genome genome_name representative_id score distance < Data/rep200.list.tbl ```
    
    * Bonus: Can you write your own similar test for `command_line_fu.py`? 

8. Once you have verified that the program succeeded in both 6 and 7, go ahead and save it as "cmd_tsv_select_columns.py" to remind you that this one uses the ComManD line. Congratulations on making your own program using Grimoire!

9. BONUS: Write down the adjustments that you made to your Grimoire prompt to make it more successful.You can use these tricks in future prompts to get to your solution more easily.
    * What problems should you ask Grimoire to take into account in the future? 
    * What words were confusing to Grimoire? 
    * What words made the pseudo code more understandable?

## Solution Check instructions:
If you are successful, you will have the output that matches the following:
* Step 5: ``` genome_id genome_name domain genus species rep_id score distance```
* Step 6: a warning message like: ```Warning: genome is not in the list of column names```

## Process - Step: 
Use this structure to help guide you through the programming process for future prompts
* Pseudo Code - Step 3
* Draft - Step 3 & 4
* Edits - Step 5
* Review - Step 6 & 7
* Publish - Step 8