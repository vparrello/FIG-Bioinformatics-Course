# Sample Set Creation Sorting Data 

Objective: Sort your gathered data into diseased and control groups.

Now that we have gathered all the data we need, it's time to make the data useable for our machine learning model. In order to create an accurate (can't remember the word), the machine learning model will need to know what data is from diseaed patients and what data is a control.

## Materials
[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)
[NCBI](https://www.ncbi.nlm.nih.gov/)
[SRAToolkit Download FASTA Documentation](https://www.ncbi.nlm.nih.gov/books/NBK242621/)

```
FIG-Bioinformatics-Course/
    ├── 3_Projections
        └── 3.1_Sample-Set-Creation/
            └── Sample-Set-Creation Exercise-2_Sorting_Data.md (you are here)
    Data/
```
    


## Exercise

1. Explain the difference between fastq and fasta

2. In this lesson, you will be sorting your data into two categories, and will need places to put all that data you will sort. You will be putting this data into two directories in the data folder of your course: one folder for 'Diseased' and one folder for 'Controls'. 
To create these directories, click on the data directory so you can see all the files and subdirectpries within (the data directory should also be highlighted). At the top of (what is this called, the directory thingy?), there are two buttons for adding something to your directory: new file and new folder. 
You will click on New Folder and name it 'Diseased', then create another folder and name it 'Control'. Make sure these two folders are under the Data directory! Wherever your highlight is, that is where the new folder will go under!

2. (Command line way) Open your terminal within VSCode. Start by asking Grimoire how to create folders using the command line. Use the prompts (he?it?) gives you to create a 'Diseased' and a 'Controls' folder. To check that they were created, use the command ls and your new folders should appear!
(possible Grimoire prompt: in visual studio code, how can you create folders using the command line?)

3. Now that you have downloaded your data, and created the folders it will go in, you need to start sorting it! By now, all your samples should be stored in the Data folder of this course. Open your samples inside VSCode. If you don't already have the extension Rainbow CSV, go to your Extensions tab and install it; this will make combing throught the file much eaiser!

# TODO
Create two directories in the data folder of the course
    use command line (or goey) (Q: what is the command line prompt to create folders? Ask Grim???)
Separate into control and diseased groups
    Tag samples from the metadata
    Also could be in the name when doing the search
    Look inside PRJDB9292(BioProject) • DRP008133(SRA Study) • All experiments • All runs( This is what we want)
        MetaData button in here - Consistent data for most projects
        Side Note: They all get the same name when doanloaded. Rename them ASAP so that you don't lose your data from ambiguous names. MENTION IT AGAIN
        Mention word wrap in the view menu and rainbow CSV extension 
Roughly 1 GB per 10 samples (Have students check their space on their computer)
Need minimum of 100 samples, 50 from each group. (Can we do this from one paper?)

Use fastq (explain the difference between the two files)