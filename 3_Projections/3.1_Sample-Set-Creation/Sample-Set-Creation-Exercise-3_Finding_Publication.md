# Sample Set Creation Exercise 3: Finding Publications

Objective: Find a publication that is relevant to your topic.

Sometimes finding data is the hardest part of getting a projection started. Especially since machine learning requires a lot of data to be effective. Therefore, we will be trying to find a publication that is relevant to your topic and then use the data from that publication to create a model.  

Many times a run of experiments are created to answer a theoretical question, or even fill out a gap in the current knowledge. This means that the data does not come from a publication, has less review, and does not necessarily need to have accurate or complete metadata before it is cataloged in the ncbi database. However publications are trying to defend their hypothesis, and therefore are more likely to have accurate and complete metadata. This is why we will be looking for an academic publication that is relevant to your topic because the data is more likely to be defined as diseased or control.

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
    └── Controls/
    └── Diseased/
```

This lesson will accomplish the same goal as the previous lesson, but it will be done in a different way. Therefore if your data is fully sorted and you have a sufficient amount of data, then feel free to move on to the next lesson.


## Exercise

*ALWAYS RESET YOUR PATH* 

Type `source ~/.bashrc` into your command line to reset your path to the Course directory before starting each exercise.

1. Travel to the http://www.ncbi.nlm.nih.gov/ and search for your topic. This time, click on the PubMed link.

2. A list of papers will now show up on your screen. There are a few things to look out for and keep in mind when looking at these papers. 
    1. Many authors will write multiple papers on the same topic. When they do this, they have a tendency to reuse old data. Be sure to check against your own findings to make sure that any data you find is new. Also avoid looking into papers that share the same authors. 
    2. Not every paper shared in this search is free to the public. It is entirely possible that you will click on a link and find that you need to pay for the paper. If this is the case, move on to the next paper. 
    3. Not every paper will have a Data Availability statement. If you find a paper that does not have a Data Availability statement, you'll have to read more deeply into the paper to see if they mention a BioProject number, SRA Experiment ID, or other information mentioned from our SRA search.

3. As an example, let's look at a paper with a Data Availability statement so we know how to find it. Go to this paper. https://pmc.ncbi.nlm.nih.gov/articles/PMC9663292/

4. This paper written by Zachary Wallen among others mentions that the gut microbiome is a biomarker for Parkinson's Disease. If you scroll down through the paper, you will notice that it has a Data Availability statement. This means that the authors have made data, but it does not guarantee that the data is in the NCBI database. 

5. You'll see that the Data Availability statement mentions a BioProject ID explicitly. It also mentions the Sequence Read Archive (SRA) and gives the project number with a link attached, This is the best case scenario. It gives us multiple sources of information we can use to track down the data and its metadata. Click on the PRJNA number.

6. Scroll down to find the "Project Data" section. This gives you the number of samples in the project. Click on the number in the SRA experiments row to see the SRAs. 

7. This will give you a list of all the samples from that project in a search. Notice that this author gave us the courtesy of labelling each sample with the disease state in its name. Not every author will do this, so I am going to continue to find the metadata. Plus it will make our lives easier when we are sorting the data later. 

8. Click on a single sample and see the Study, Sample, Library, and Runs sections. Click on the "all runs" link in the Study section to grab all the samples from that study.

9. This page has everything we want. The metadata for the process, information on collection, and even all the metadata in a list. Go ahead and find the Select box and download the metadata using the "Metadata" button. 

Notice that the table you just downloaded is also pictured on this website. You can see it has columns like "Age_at_collection", "Bases", "Bytes", and "Case_status". This Case_status column is what we want. It tells us the disease state of the sample. So even if the author had not named each sample with its disease status, we can use this column to help us sort the data.

10. Go ahead and open your new SraRunTable.csv file. Notice that you now have 700+ samples that you could use to run your model. This is the power of finding a publication that is relevant to your topic. They will do a lot of the work of sorting your data for you. However, keep in mind that very few publications will have a Data Availability statement or data that is available in the NCBI database. Some authors don't even have their publications out there for free. So this method of finding data is sometimes frustrating and can take a lot of time to get the results you want. 
