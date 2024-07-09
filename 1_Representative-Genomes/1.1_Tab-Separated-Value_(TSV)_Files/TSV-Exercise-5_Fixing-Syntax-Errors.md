# TSV File Exercise 5 - Fixing Syntax Errors

Objective: Find and eliminate syntax errors inside of your code

In almost every program that is written by a human, there will usually be an error or two that will be accidentally incorporated, such as spelling or minor syntax errors. VSCode (and other IDEs) may be able to detect and signal the program-developer that there are potential errors within their code based on the syntax rules of the programming language. It is important to not only recognize these signals, but to also understand the error messages that you might receive while writing or running code. 

## Materials: 
[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

```
FIG-Bioinformatics-Course/
├── 1_Representative-Genomes/
│   └── TSV-Exercise-5_Fixing-Syntax-Errors.md   (you are here)
├── Code/
│   ├── error_message1.py
│   ├── error_message2.py
│   ├── error_message3.py
│   ├── error_message4.py
│   ├── error_message5.py
│   └── error_message_final.py
└── Data/
    ├── bindict.tbl
    ├── data.tbl
    └── rep200.list.tbl   
```

## Exercise: 

1. Ask Grimoire to describe to you each of the following errors. Make sure that you ask follow up questions and request examples for each type of error that you think would be hard to identify.
* Syntax errors
* Runtime errors
    * Index errors
    * Name errors
    * Attribute errors
    * Type errors

2. Inside the Code directory, we have include a set of files that are designed to throw specific error-messages. Each filename has the form `error_messageN.py` where `N` is a number. Open up the first file `error_message1.py`, which contains 3 syntax errors; run this program to see which error messages are reported. 

3. Copy each error-message and paste it into Grimoire. Then ask Grimoire to explain each associated error to you and to recommend how you can fix it. Apply these fixes to your code and run it again. If you get another error message, follow the same workflow by pasting the message into Grimoire and asking it to help you fix it. Feel free to try and fix the errors on your own as an extra challenge.

4. Each file has been constructed to throw a different error. Most of them require a data file argument when being called from the command line. If you get a message that looks like:
```
Usage: python tsv_headers.py <filename>
```
then add `Data/data.tbl` to the end of your command to see if the program now works. Use the table below to keep track of which errors appeared and how to fix them. You can use this table as "quick reference" when debugging future programs.

** NOTE: If you are getting an error that says "[Errno 2] No such file or directory" for the error_message.py file you are working on. Try killing your terminal and opening a new one to ensure the error is coming from the program and not from your operating system.

| Error Type      | File Name | Possible Fix |
| :-------------- | --------- | ------------ |
| Index Error     |           |              |
| Name Error      |           |              |
| Attribute Error |           |              |
| Type Error      |           |              |
| Syntax Error    |           |              |


## Solution Check instructions:
If you are successful, you should see the following output for each associated command:

* error_message1.py: Add a colon to the end of lines 5, 12, and 16
* error_message2.py: Add `filename = sys.argv[1]` to line 10
* error_message3.py: Remove line 11 `y='5'`
* error_message4.py: Remove the `counter = 1` and `Next Up` message (Lines 10, 21, 22)
* error_message5.py: Remove Line 20 `header.reverse()`
* error_message_final.py: Put line 12 after line 14.
    * remove `reverse()` from line 12
    * Add a `"` to line 20 just before the closing paranthesis near the end
    * Add a `:` to the end of line 4