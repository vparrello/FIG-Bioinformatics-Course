#### TSV File Exercise 4 - Documentation and its uses

 Objective: Create documentation for your program for reference

Now that we have created our program, we will need to adjust it for outside use.  This means making sure that anyone who uses it can easily understand how it works as well as making sure it doesn't break in the process. We can accomplish this by catching any incorrect uses of the program and giving warning messages against it and creating documentation around how the program works. This exercise walks you through both of these processes while also adding functionality the user will be expecting.


#### Materials: 

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)
```
    FIG-Bioinformatics-Course/
        1 - Representative Genomes/
            Data/bindict.tbl
            Data/data.tbl
            Data/rep200.list.tbl
            Data/TSVReaderDocumentation.txt
            bin/tsv_headers.py
```

#### Exercise: 

1. Before we work on documentation for our code, we must figure out how to make sure the user does not break it.  This means solving for use cases that you might not have thought of. Give Grimoire a new copy of your program. Then tell it that you need to find 5 ways that a user might break the program (create an error message.)

2. Ask Grimoire to give recommendations on how it woud solve for each of these scenarios inside the program that it has. Be sure to ask it to explain how it solves the problem "line by line" to ensure you understand all of its solutions. 

3. Once your program has been adjusted this way, we now need to document what your program does and how people can use it. There are three places to add documentation.
    1. Inside the code. This is comments in strategic places to make sure people reading, editing, or working on the code can understand it easily.
    2. Inside the usage of the code. Think of this as the help menu in applications. We need to build this menu for the user to be able to access from inside the command line.
    3. In a completely separate document. This is the instruction manual on how the code works, what it is expected to do, and the possible use cases of the code.

4. Grimoire already has your code. Ask it to add comments in key places inside your code to explain what it does. You can also attempt to do this by hand to increase your understanding

5. Inside your terminal in VScode, try accessing the help menu by using the following command:

    ``` python3 bin/tsv_headers.py --help ```

    This should display a message of what every variable does and which ones are optional. Copy and paste the output into Grimoire and tell it that you need to adjust the aforementioned code to include help messages inside this output. Use that prompt to adjust your command line help menu.

6. Lastly we need to write our own documentation inside our own file. Look for a file called "TSVReaderDocumentation.txt" inside the Data directory. Using Grimoire, ask it to write documentation for the code that describes the following things.

    * Purpose of the program
    * Functions of the different variables
    * Troubleshooting recommendations
    * Author credits


## Solution Check instructions:
If you are successful, you will have documentation of your program in TSVReaderDocumentation.txt, the ```python3 tsv_header --help``` menu, and at least 10 comments inside of your code. 

