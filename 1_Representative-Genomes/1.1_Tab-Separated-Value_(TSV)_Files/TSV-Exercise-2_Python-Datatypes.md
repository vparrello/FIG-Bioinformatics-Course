#### Tab Separated Value Files Exercise 2 - Python Datatypes

Objective: Recognize Python Datatypes

Python is a programming language that lends itself to the large consumption of data. In this exercise, we will explore the different types of data python uses within its code. Most of these data types are common among all programming languages but a few are unique to Python itself. 

The program "data_types.py" has been created to help you explore these data types while also practicing the use of the command line. Before you begin, ensure that you have a Terminal window open of the GitBash profile. If you are unsure of how to find this setting, please refer to the "Download Instructions" around Git to make it a default setting for VSCode.

#### Materials: 
[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

* FIG-Bioinformatics-Course/
    * 1_Representative-Genomes/
        * 1.1_Tab-Separated-Value_(TSV)_Files/
            * Code/
                * data_types.py
        * Data/
            * rep200.list.tbl

#### Exercise: 

1. The following is the list of data types that we will be discussing in this lesson. Ask Grimoire to explain each type and its uses to you. Make sure you specify that they are Python Datatypes.
    * Boolean
    * Integer
    * String
    * List
    * Set
    * Tuple
    * Dictionary

2. A program called "data_types.py" has been provided to help you see how each of these datatypes can be used with bioinformatics data. This program takes the file "rep200.list.tbl", which contains information about representative genomes, and converts pieces of it into different data types. Open up the file and look through the code to see what it does.

3. Call the program on each of the following datatypes to see examples of the data populate. These three are the most basic types and all subsequent data types will use a combination of these to compose the rest of the datatypes.
    * `python 1_Representative-Genomes/1.1_Tab-Separated-Value_\(TSV\)_Files/Code/data_types.py boolean`
    * `python 1_Representative-Genomes/1.1_Tab-Separated-Value_\(TSV\)_Files/Code/data_types.py integer`
    * `python 1_Representative-Genomes/1.1_Tab-Separated-Value_\(TSV\)_Files/Code/data_types.py string`

4. When you call the program, it reads the file that it is given and turns that file into a list. That list is then looped through one by one and separated into its own individual columns. We can do this because the file is in the Tab Separated Value (TSV) format. We also added some extra data at the end to help show the difference between lists and sets. Print that list by calling the "data_types.py" program on the datatype "list". 
    *** Hint: You can access previous commands in the command line by using the arrow keys on your keyboard.

5. Next the program takes that list and inserted it into a set. See if you can spot the difference between a list and a set. Print that set by calling the "data_types.py" program on the datatype "set".

6. Ask Grimoire to explain to you the difference between a list, a set, and a tuple. The tuple's main feature is that it is immutable. Ask Grimoire to explain to you what "immutable data" is and why it is important to use while programming code.

7. Print the tuple by calling the "data_types.py" program on the datatype "tuple".

8. The last datatype is called a dictionary. This datatype acts like a traditional dictionary because you can "look up" data by using a key that will always return a specific value. The "key" in this sense is a string or integer that is unique to that data. The value can be any other datatype including a dictionary. Ask Grimoire to explain to you what the key value pairs are inside of a python dictionary. 

9. Print the dictionary by calling the "data_types.py" program on the datatype "dictionary". Notice how each key and value are unique.

## Solution Check instructions:
If you are successful, the program should return output that matches the following.

3. ```
    True
    15
    910454.3        Prochloron didemni P2-Fiji      Bacteria        1214    1216    1094894.3       315     0.09

    ```
4. ```
    ['910454.3', 'Prochloron didemni P2-Fiji', 'Bacteria', '1214', '1216', '1094894.3', '315', '0.09', 'Foo', '315', '1214']
    ```
5. ```
    {'Prochloron didemni P2-Fiji', '315', '1094894.3', '1214', 'Foo', '910454.3', '1216', 'Bacteria', '0.09'}
    ```
7. ```
    ('910454.3', 'Prochloron didemni P2-Fiji')
    ```
9. ```
    {'511145.12': 'Escherichia coli str. K-12 substr. MG1655', '107806.10': 'Buchnera aphidicola str. APS (Acyrthosiphon pisum)', '685038.3': 'Escherichia coli O83:H1 str. NRG 857C', '1028307.3': 'Enterobacter aerogenes KCTC 2190', '585057.6': 'Escherichia coli IAI39', '99287.12': 'Salmonella enterica subsp. enterica serovar Typhimurium str. LT2', '568707.3': 'Bordetella bronchiseptica 253', '300267.13': 'Shigella dysenteriae Sd197', '160490.10': 'Streptococcus pyogenes M1 GAS', '585056.7': 'Escherichia coli UMN026', '386585.9': 'Escherichia coli O157:H7 str. Sakai', '1133852.3': 'Escherichia coli O104:H4 str. 2011C-3493', '71421.8': 'Haemophilus influenzae Rd KW20', '272947.5': 'Rickettsia prowazekii str. Madrid E', '393305.7': 'Yersinia enterocolitica subsp. enterocolitica 8081', '243274.5': 'Thermotoga maritima MSB8', '220341.7': 'Salmonella enterica subsp. enterica serovar Typhi str. CT18', '171101.6': 'Streptococcus pneumoniae R6', '210007.7': 'Streptococcus mutans UA159', '169963.11': 'Listeria monocytogenes EGD-e', '568814.3': 'Streptococcus suis BM407', '1208660.3': 'Bordetella parapertussis Bpp5', '208435.3': 'Streptococcus agalactiae 2603V/R', '312309.11': 'Vibrio fischeri ES114', '214092.21': 'Yersinia pestis CO92', '194439.7': 'Chlorobium tepidum TLS'}
    ```
