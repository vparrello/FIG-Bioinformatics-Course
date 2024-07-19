# Tab Separated Value File Exercise 2 - Python Datatypes

Objective: Recognize Python Datatypes

The programs that you will be creating in this course will be written in the `Python` computer-language.
In this exercise, we will explore the different types of data that Python uses within its code. Most of these data types are common to all programming languages.

The program `python_datatypes_interactive_exercise.py` has been created to help you explore these data types while also practicing the use of the command line. Before you begin, ensure that you have a Terminal window open that uses the `bash` profile (macOS and LINUX) or `GitBash` (Windows) profile.
If you are unsure of how to find this setting, please refer to the "Download Instructions" for Git to make it the default setting for VSCode.

## Materials: 
[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)


```
FIG-Bioinformatics-Course/
├── Code/
│   └── python_datatypes_interactive_exercise.py
├── Data/
│   └── rep200.list.tbl
└── 1_Representative-Genomes/
    └── 1.1_Tab-Separated-Value_(TSV)_Files/
        ├── TSV-Exercise-Exercise-2_Python-Datatypes.md  (you are here)
        └── Solutions/
            └── tsv_headers_solution.py
```


## Exercises: 

1. The following is the list of data types that we will be discussing in this lesson. Ask Grimoire to explain each type and its uses to you. Make sure you specify that they are Python Datatypes.
    * Boolean
    * Integer
    * Float
    * String
    * List
    * Set
    * Tuple
    * Dictionary

2. A program called `python_datatypes_interactive_exercise.py` has been provided that will return the above list, and that for each datatype will provide a short definition for that type of data, followed by a code example that uses that type of data within a bioinformatics contex. You can run this program at any time during the course if you need to refresh you memory regarding a particular datatype and how it is used within the code. If the example provided within any section of this program is not clear to you, you can paste that section into Grimoire and ask it to explain the code to you "line-by-line". Likewise, if what any line does is not clear to you, ask Grimoire to explain that line to you "step-by-step".

3. Three of the most basic data-types are `boolean`, `integer`, and `string`. More complex data-types can be constructed from these basic data-types. Run the program `python_datatypes_interactive_exercise.py` for each of these three basic datatypes to see examples of the data that they contain:

* `python Code/python_datatypes_interactive_exercise.py boolean`

* `python Code/python_datatypes_interactive_exercise.py integer`

* `python Code/python_datatypes_interactive_exercise.py string`

NOTE: The words that follow a program's name are called `arguments`. The `arguments` of a program are passed to the program to control its behavior or to provide it with necessary information. If you are unfamiliar with the concept of a "program argument", please ask Grimoire to explain it to you.

4. Run the program `python_datatypes_interactive_exercise.py`, with an argument of `list`.
    * Hint: You can access and move back and forward through the previous commands you have entered onto the command line by using the up-arrow and down-arrow keys on your keyboard, and you can edit a previous command before resubmitting it. This trick will save you lots of typing!

Notice that one of the list examples has a duplicated data item,
and that both instances of the data-item are printed out. 
    

5. Next invoke the program with an argument of `set`.
In "Example 1", see if you can spot the difference between a `list` and a `set`. 
Notice that this time, no data is duplicated.
Notice also that the `list` is displayed with square brackets,
while the `set` is displayed with curly brackets.
"Example 2" shows that a list can be converted to a set using the `set()` function.

6. Ask Grimoire to explain to you the difference between a `list`, a `set`, and a `tuple`. The tuple's main feature is that it is immutable. Ask Grimoire to explain to you what "immutable data" is and why it is important to use within programming code.

7. Print the tuple by calling the `python_datatypes_interactive_exercise.py` program on the datatype `tuple`.

8. The last datatype is called a `dictionary`. This datatype acts like a traditional dictionary because you can "look up" data by using a "key" that will return an associated "value". For our purposes, the "key" in this sense is a string or an integer that provides a "name" or "identifier" that can be used to refer to its associated data-value. (Python does allow some data-types besides strings and integers to be used as "keys", but the rules governing keys are complicated, and will not be needed for this course.) The value associated with a key can be any python datatype, including sets, lists, or even another dictionary, which allows one to build up an arbitrariy complex data-structure. Ask Grimoire to explain to you what is meant by "key-value" pairs with a python dictionary. 

9. Print the dictionary by calling the `python_datatypes_interactive_exercise.py` program on the datatype `dictionary`. Notice how each "entry" (i.e., key-value pair) is unique.

IMPORTANT: Key-value pairs are entered into a `dictionary` by entering the "key" and the "value" separated by a colon. Each "key" can only occur in a dictionary once, but nothing restricts the associated values. For example, in a dictionary listing types of food, the keys "apple" and "pear" can both have the value of "fruit". Example: {"apple": "fruit", "pear": "fruit", "carrot": "vegetable"}.
Reentering an old key with a new value will update the key's value.

10. Because "dictionaries" are such an important tool for writing programs, we have included a separate interactive exercise on use of dictionaries. Please launch this exercise as follows:
```
    python3 Code/dictionary_interactive_exercise.py
```
The execise will welcome you, provide a brief description of what a "dictionary" is, and then prompt you regarding whether you wish to `add a key-value pair`, `delete a key`, or `quit`. After each operation that you enter, the program will display the new contents of the dictionary. Please experiment by adding an assortment of key-value pairs, deleting keys (and their associated values),
and re-adding old keys with new values, and observe how the contents of the dictionary changes after each operation. Once you think you have a good understanding of how the "add" and "delete" operations alter the contents of the dictionary, you can leave the exercise by typeing `q` or `quit`.

We encourage you to open `Code/dictionary_interactive_exercise.py` and have a look inside. if there is any part of the code that you feel you don't understand, you can paste it into Grimoire and ask Grimoire to explain the code to you "line-by-line". If there are parts that you feel you still don't understand, ask Grimoire to explain those lines "step-by-step".

## Solution Check instructions:
If you are successful, the program should return output that matches the following.

3. ```
    True
    41
    194439.7        Chlorobium tepidum TLS  Bacteria        256319  1097    194439.7        9999    0.00
    ```
4. ```
    ['194439.7', 'Chlorobium tepidum TLS', 'Bacteria', '256319', '1097', '194439.7', '9999', '0.00', 'Foo', 'Foo', 'Foo', '315', '1214']
    ```
5. ```
    {'Chlorobium tepidum TLS', '1097', '315', '256319', '0.00', 'Bacteria', '1214', '9999', '194439.7', 'Foo'}
    ```
7. ```
    ('194439.7', 'Chlorobium tepidum TLS')
    ```
9. ```
    {'511145.12': 'Escherichia coli str. K-12 substr. MG1655', '107806.10': 'Buchnera aphidicola str. APS (Acyrthosiphon pisum)', '685038.3': 'Escherichia coli O83:H1 str. NRG 857C', '1028307.3': 'Enterobacter aerogenes KCTC 2190', '585057.6': 'Escherichia coli IAI39', '99287.12': 'Salmonella enterica subsp. enterica serovar Typhimurium str. LT2', '568707.3': 'Bordetella bronchiseptica 253', '300267.13': 'Shigella dysenteriae Sd197', '160490.10': 'Streptococcus pyogenes M1 GAS', '585056.7': 'Escherichia coli UMN026', '386585.9': 'Escherichia coli O157:H7 str. Sakai', '1133852.3': 'Escherichia coli O104:H4 str. 2011C-3493', '71421.8': 'Haemophilus influenzae Rd KW20', '272947.5': 'Rickettsia prowazekii str. Madrid E', '393305.7': 'Yersinia enterocolitica subsp. enterocolitica 8081', '243274.5': 'Thermotoga maritima MSB8', '220341.7': 'Salmonella enterica subsp. enterica serovar Typhi str. CT18', '171101.6': 'Streptococcus pneumoniae R6', '210007.7': 'Streptococcus mutans UA159', '169963.11': 'Listeria monocytogenes EGD-e', '568814.3': 'Streptococcus suis BM407', '1208660.3': 'Bordetella parapertussis Bpp5', '208435.3': 'Streptococcus agalactiae 2603V/R', '312309.11': 'Vibrio fischeri ES114', '214092.21': 'Yersinia pestis CO92', '194439.7': 'Chlorobium tepidum TLS'}
    ```