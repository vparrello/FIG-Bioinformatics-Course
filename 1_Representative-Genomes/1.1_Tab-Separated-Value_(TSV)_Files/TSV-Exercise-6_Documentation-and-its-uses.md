#### TSV File Exercise 6 - Documentation and its uses

 Objective: Create documentation for your program for reference

Now that we have created our program, we will need to adjust it for outside use.  This means making sure that anyone who uses the program can easily understand how to invoke it, as well as making sure that it doesn't break due to incorrect inputs. We can accomplish this by catching incorrect use of the program and giving the user warning messages, and creating documentation that shows the user how to invoke  the program. This exercise walks you through the error-checking and documentation processes, and also adds a capability that the user may be expecting.


#### Materials: 

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)


* FIG-Bioinformatics-Course/
    * 1_Representative-Genomes/
        * Data/
            * bindict.tbl
            * data.tbl
            * rep200.list.tbl
            * TSVReaderDocumentation.txt
        * Code/
            * tsv_headers.py


#### Exercise: 

1. Before we work on documentation for our code, we must figure out how to make sure the user does not break it.  This means giving careful consideration to possible incorrect invocations that you did not intend. Tell Grimoire that you are attaching a program (see "file attachment" in the previous exercise), and that you wish to protect against possible invocation errors such as missing mandatory arguments, nonexistent data-files, etc. Ask Grimoire to give recommendations on how to protect against each of the possible error-producing scenerios that it finds within the program. Be sure to ask it to explain how it has solved each problem "line by line" to ensure that you understand all of its suggested improvements. 

2. Once your program has been adjusted to protect against invocation errors, we now need to document what your program does and how people can use it. There are three places to add documentation:
    * Inside the code, by inserting comments at strategic points to make sure that people reading, editing, or working on the code can understand it easily.
    * Inside the usage or "help message" for the code. Think of this as being like the "help menu" in an application. We need to build a menu of command options that the user can request via the command line.
    * In a completely separate document. This is the instruction manual on how the code works, what it is expected to do, and the possible use cases of the code.

4. Grimoire already has your code. Ask it to insert comments at key points within your code to explain what each section does. You can also attempt to do this by hand to increase your understanding.

5. Ask Grimoire to insert a "help block" for your code that will be returned if the code is invoked with a `-h` or `--help` argument.

6. Lastly we need to write our own documentation inside our own file. Look for a file called "TSVReaderDocumentation.txt" inside the Data directory. Using Grimoire, ask it to write documentation for the code that describes the following things.

    * Purpose of the program
    * Functions of the different variables
    * Troubleshooting recommendations
    * Author credits


## Solution Check instructions:
If you are successful, you will have documentation of your program in TSVReaderDocumentation.txt, an added `--help` option that when invoked will explain to a user what the program does and how to use it, and a set of comments inserted into your code that will explain to another programmer what each section of the code does.
