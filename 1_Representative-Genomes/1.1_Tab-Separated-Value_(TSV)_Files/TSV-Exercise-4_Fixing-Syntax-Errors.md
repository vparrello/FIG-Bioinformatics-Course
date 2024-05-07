#### TSV File Exercise 4 - Fixing Syntax Errors

Objective: Find and eliminate syntax errors inside of your code

In every program that is written, there is always an error or two that comes up. VSCode (and other IDEs) often have signals to the user on what might be potential errors inside the code according to the rules of the programming language. It is important to not only recognize these signals, but to also understand the error messages that you might recieve while writing code. 

#### Materials: 
[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

* FIG-Bioinformatics-Course/
    * 1_Representative-Genomes/
        * Data/
            * bindict.tbl
            * data.tbl
            * rep200.list.tbl
        * Code/
            * error_message1.py
            * error_message2.py
            * error_message3.py


#### Exercise: 

1. Ask Grimoire to describe to you the following errors. Ensure that you ask follow up questions and examples for any type of error that you think would be hard to identify.
* Syntax errors
* Runtime errors
    * Name errors
    * Type errors
    * Index errors
    * Attribute errors
* Logical errors

2. Inside the Code directory, there are a few files that are set up to throw error messages. They are all labelled 'error_messageN.py' where N is a number. Open up the first file 'error_message1.py' to find the 3 Syntax Errors present. In order to get more information, run the program and see what error message comes out. 

3. Copy the error and paste it into Grimoire. Then ask Grimoire to explain to you the error and how to fix it. Apply this fix to your code and run it again. If you got another error message, follow the same flow by pasting it into Grimoire and asking it to help you fix it. Feel free to try and fix them on your own as an extra challenge.

4. Each file has been set up to throw a different error. Use the table below to Keep track of which errors appear and how to fix them. You can use this as reference in future programs.

| Error Type      | File Name | Possible Fix |
| :-------------- | --------- | ------------ |
| Syntax Error    |           |              |
| Name Error      |           |              |
| Type Error      |           |              |
| Index Error     |           |              |
| Attribute Error |           |              |
| Logical Error   |           |              |


## Solution Check instructions:
If you are successful, you should see the following output for each associated command:

error_message1.py: Add a colon to the end of lines 5, 12, and 16
error_message2.py: Add 'filename = sys.argv[1]' to line 10
error_message3.py: Remove line 11 "y='5'"
error_message4.py: Remove the counter and Next Up Message (Lines 10, 21, 22)
error_message5.py: Remove Line 20 'header.reverse()'
error_message_final.py: Put line 12 after line 14.
    remove "reverse()" from line 12
    Add a " to line 20 at the end before the paranthesis
    Add a : to the end of line 4