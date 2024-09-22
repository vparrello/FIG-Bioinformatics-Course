# Sample Set Creation by Gathering Data

Objective: Learn how to gather data from the NCBI database to create a sample set of genomes for a specific research topic.

Bioinformatics research is about establishing patterns in genetic material that can then be used to make predictions about the behavior of other genomes. Traditionally, this is done by creating a model that uses features of the genomes to make a prediction. The features are often derived from the genetic material itself, such as k-mers, or from other sources, such as gene expression data. Since these patterns can be so small, it is often difficult to find them using human analysis, so we use bioinformatics tools to help us find them.

One of these tools is the representative genome set that you created in the first unit. A second one can be the hammer set that you created in the second unit. Another can be the bv-brc.org and the multitude of tools available in the sratoolkit you have also installed throughout this course.We are now going to use these tools to create an even more powerful tool: a machine learning model that we can use to either diagnose a disease or create a treatment to target the disease. To do this, we need to gather data from the NCBI database that is relevant to the research topic we are interested in.

## Materials
[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)
[NCBI](https://www.ncbi.nlm.nih.gov/)
[SRAToolkit Download FASTA Documentation](https://www.ncbi.nlm.nih.gov/books/NBK242621/)

FIG-Bioinformatics-Course/
├── 3_Projections
    └── 3.1_Sample-Set-Creation/
        └── Sample-Set-Creation Exercise-1_Gather_Data.md (you are here)

## Exercise

1. Ask Grimoire to describe the different types of machine learning models that can be used to predict the presence of a disease. Specifically mention classification, random forest, and sample set creation.

2. The machine learning model that we are going to create is a classification model. This model basically takes a set of samples that have been categorized into diseased or healthy and uses that as a key to determine what the pattern for each group is. Then is uses the genomes found within a mystery sample to make a prediction about the health status of the sample. Kind of like we were doing in a smaller scale in the previous units, we will be doing it on a larger scale and with a lot more data. This increases the accuracy of the model and allows us to make more accurate predictions. Below are a list of research topics that can be applicable to our objective of classification. Choose one of the topics below. If you would like to choose a different topic, do a little research first to ensure that the gut microbiome has been studied for that specific topic. Otherwise it is possible that your machine learning model will not be able to distinguish between the two groups and will not be able to make accurate predictions.

   1. **Obesity**
   2. **Type 2 Diabetes**
   3. **Inflammatory Bowel Disease (IBD)**
   4. **Mental Health**
   5. **Dietary Allergies**
   6. **Autoimmune Diseases**
   7. **Alzheimer's Disease**
   8. **Parkinson's Disease**

3. Now that we have chosen a research topic, we need to gather data from the NCBI database. Please go to https://www.ncbi.nlm.nih.gov/ and search for the topic you chose. This will give you a few different things that are relevant to our objective. 

![alt text](NCBISearch.png)

Firstly, we have the Literature section. This section will give you research papers and articles that are relevant to your topic. It can give you a good idea of what is currently known about the topic and can help you understand the current state of research. It can also lead us to where data is available within the Data sections of each academic paper published on the topic. If you want to find data this way, click on "PubMed" and a search with all of the papers relevant to your topic will be published. Make sure you pay attention to the data availability section of each paper as well as the authors. Often times an author will use the same data set for multiple publications and if they do you might end up with duplicate data if you use both of their papers.

Secondly we have the Genomes section. This section will give you any and all genomes that have been tagged with your specific topic. Depending on your need, you can choose to look at the BioProject which has all the samples from a multiple research project, or a BioSample which contains individual samples from one experiment, or the SRA which contains all of the raw data from a single sample. Ask Grimoire to describe any other sections of this page that you find interesting.

Click on the SRA heading to see all the data we are looking for.
![alt text](SRA.png)

4. Inside of this search, you will see a list of different SRA samples that are relevant to your topic. Each one has a link, a description of metadata, and an Accession number underneath it. Ask Grimoire to describe what "spots" and "bases" are referring to in terms of a metagenome sample. 

5. Good Data is important. It keeps our classifier accurate and our predictions useful. However, having an overwhelming abundance of data is the only way to ensure that our model can identify enough of a pattern to make an accurate prediction. In this way, we need to find a balance between the number of samples that we use and the quality of the data that we gather. If you click on any of the links inside of your search, you will find a page with all the metadata associated with that sample. 

![alt text](SRAMetaData.png)

