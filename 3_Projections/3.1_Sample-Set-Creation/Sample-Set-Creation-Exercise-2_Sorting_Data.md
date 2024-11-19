# Sample Set Creation Sorting Data 

Objective: Sort your gathered data into diseased and control groups.

Now that we have gathered all the data we need, it's time to make the data useable for our machine learning model. In order to create a classifier, 
the machine learning model will need to know which data are from diseaed patients
and which data are from the controls.

## Materials

* [Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)

* [NCBI](https://www.ncbi.nlm.nih.gov/)

* [SRAToolkit Download FASTA Documentation](https://www.ncbi.nlm.nih.gov/books/NBK242621/)

```
FIG-Bioinformatics-Course/
    ├── 3_Projections
    │   └── 3.1_Sample-Set-Creation/
    │       └── Sample-Set-Creation Exercise-2_Sorting_Data.md (you are here)
    └── Data/
```
    


## Exercise

1. Despite the similarities in their names, the `FASTQ` and `FASTA` file-formats
are not related, and they serve different purposes.
`FASTA` files only store sequence-data, while `FASTQ` files store 
both sequence data and information about the "quality" (degree of confidence)
of the sequence data, and hence the "Q" in its name.
Ask Grimoire to explain the difference between the FASTQ and FASTA file-formats,
and the respective use-cases for each format.

2. In this lesson, you will be sorting your data into two categories, and will need places to put all that data you will sort. You will be putting your data into two directories within the `Data` folder for this course: one folder for 'Diseased' and one folder for 'Controls'. 
To create these directories using VScode, click on the `Data` directory in the File Explorer,
which will open it out to show the files and subdirectories within it (the data directory should also be highlighted). At the top of the File Explorer display, there are two buttons for adding something to the selected directory: new file and new folder, symbolized by "Document plus" and "Folder plus" icons.
You will click on New Folder and name it 'Diseased', then create another folder and name it 'Control'. Make sure these two folders are under the Data directory! Whichever directory is highlighted is where the new folder will go!
(If you accidentally create a file or folder in the wrong directory, you can drag it to the correct directory.)

* Alternatively, you can also create folders using the command-line. Open your terminal within VSCode. Start by asking Grimoire how to create a directory named 'Diseased' within the directory named 'Data' using the VScode command line,
and enter its recommendation into the terminal window.
Create the 'Controls' directory using an analogous command.
To check that these directories were created, use the command `ls Data`, and your new folders should now appear!

3. Now that you have downloaded your data, and created the folders they will go in, you need to start sorting them!
By now, all your samples should be stored in the Data folder of this course. Open your samples inside VSCode. If you don't already have the extension Rainbow CSV, go to your Extensions tab and install it; this will make combing throught the file much eaiser!

# TODO
Create two directories in the data folder of the course
    use command line (or GUI) (Q: what is the command line prompt to create folders? Ask Grim???)
Separate into control and diseased groups
    Tag samples from the metadata
    Also could be in the name when doing the search
    Look inside PRJDB9292(BioProject) • DRP008133(SRA Study) • All experiments • All runs( This is what we want)
        MetaData button in here - Consistent data for most projects
        Side Note: They all get the same name when doanloaded. Rename them ASAP so that you don't lose your data from ambiguous names. MENTION IT AGAIN
        Mention word wrap in the view menu and rainbow CSV extension 
Roughly 1 GB per 10 samples (Have students check their space on their computer)
Need minimum of 100 samples, 50 from each group. (Can we do this from one paper?)

Use FASTQ (explain the difference between the two files)