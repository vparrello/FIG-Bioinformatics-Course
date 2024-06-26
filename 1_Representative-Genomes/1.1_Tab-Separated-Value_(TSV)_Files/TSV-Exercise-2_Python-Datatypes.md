# Tab Separated Value File Exercise 2 - Python Datatypes

Objective: Recognize Python Datatypes

The programs that you will be creating in this course will be written in the `Python` computer-language.
In this exercise, we will explore the different types of data Python uses within its code. Most of these data types are common to all programming languages, but a few are unique to Python itself. 

The program `data_types.py` has been created to help you explore these data types while also practicing the use of the command line. Before you begin, ensure that you have a Terminal window open that uses the `bash` profile (macOS and LINUX) or `GitBash` (Windows) profile.
If you are unsure of how to find this setting, please refer to the "Download Instructions" for Git to make it the default setting for VSCode.

## Materials: 
[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)


```
FIG-Bioinformatics-Course/
├── Code/
│   └── data_types.py
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
    * String
    * List
    * Set
    * Tuple
    * Dictionary

2. A program called `data_types.py` has been provided to help you see how each of these datatypes can be used with bioinformatics data. This program reads the file `rep200.list.tbl`, which contains "metadata" about more than 287,000 genomes, and extracts various elements from it into different python data-types. Open up the data-file to examine its format (WARNING: this file is big!), then look through the program to see what it does. If what any section of the program does is not clear to you, paste that section into Grimoire and ask it to explain the code to you "line-by-line". Likewise, if what any line does is not clear to you, ask Grimoire to explain that line to you "step-by-step".

3. Three of the most basic data-types are `boolean`, `integer`, and `string`. More complex data-types can be constructed from these basic data-types. Run the program `data_types.py` for each of these three basic datatypes to see examples of the data that they contain:

* `python Code/data_types.py boolean`

* `python Code/data_types.py integer`

* `python Code/data_types.py string`

NOTE: The words that follow a program's name are called `arguments`. The `arguments` of a program are passed to the program to control its behavior or to provide it with necessary information. If you are unfamiliar with the concept of a "program argument", please ask Grimoire to explain it to you.

4. When you run the program `data_types.py`, it reads the file that it is given and turns that file into a list. That list is then looped through one by one and separated into its own individual columns. We can do this because the file is in the Tab Separated Value (TSV) format. We also added some extra data at the end to help show the difference between lists and sets. Print that list by calling the `data_types.py` program with an argument of `list`.
Notice that the list has a duplicate data point. 
    * Hint: You can access previous commands in the command line by using the arrow keys on your keyboard.

5. Next the program takes that list and inserts it into a `set`. See if you can spot the difference between a `list` and a `set`. Print the set by calling the `data_types.py` program using the datatype-argument `set`, just as you did for the previous three datatypes. Notice that at this point, no data is duplicated.

6. Ask Grimoire to explain to you the difference between a `list`, a `set`, and a `tuple`. The tuple's main feature is that it is immutable. Ask Grimoire to explain to you what "immutable data" is and why it is important to use within programming code.

7. Print the tuple by calling the `data_types.py` program on the datatype `tuple`.

8. The last datatype is called a `dictionary`. This datatype acts like a traditional dictionary because you can "look up" data by using a "key" that will return an associated "value". The "key" in this sense is a string or integer that provides a "name" or "identifier" referring to its associated data. The value can be any other datatype, including another dictionary. Ask Grimoire to explain to you what is meant by "key-value" pairs with a python dictionary. 

9. Print the dictionary by calling the `data_types.py` program on the datatype `dictionary`. Notice how each "entry" (i.e., key-value pair) is unique.

IMPORTANT: Each "key" can only occur in a dictionary once, but nothing restricts the associated values. For example, in a dictionary listing types of food, the keys "apple" and "pear" can both have the value of "fruit". Example: {"apple": "fruit", "pear": "fruit}

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