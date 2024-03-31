<!--
Created by: Victoria Parrello
Last Updated: 2/8/2024
-->

# Introduction

This course is meant to give an introduction to AI-assisted code-development via the exploration of three major topics in applied bioinformatics. Each topic has been organized into a module consisting of a number of related subtopics and associated exercises. You will find the exercises inside each topic's folder, along with any data that you might need inside that topic's Data-folder, and related code inside that topics "bin" folder. If you need any help or are not sure how to continue, you are encouraged to use AI tools to help you throughout your learning journey.

## Content

### 1. Representative Genomes  
1. Construct and learn to use command-line tools to manipulate tab-separated files.
2. Construct and use command-line tools to read, process, and write DNA sequence-data in FASTA format
3. K-mers and Jaccard Similarities - Use your knowledge of genome space to construct K-mers that can compare two different sequences
3. Singly-Occuring Universal Roles (SOURs) - Understand that a set of SOURs can act as unique genome signatures that can be used to pick out the Representative Genomes
4. Create a Representative Genome Set - Compare the genomes in your universe to find a subset of "Representative Genomes" that "covers" the full set.
### 2. Hammers  
1.  Strike the Role! - Choose a single repgen genome and discover how to find out what the names of all its neighboring peer genomes are. How can you prove they are all of the same represented set?
2.  Create your own hammer set - Use the role you found to create a "hammer set" of signatures used to find the nearest repgen genome. Experiment with using many roles to create generalized "hammer sets".
3. Use the hammers - Identify the genomes contained within a sample using the hammers from your repgen set. Which repgen set made the hammers most effective?
### 3. Projections  
1.  Sample Set Creation - Create the sample report needed to feed the hammer output into an "X-matrix"
2.  X-Matrix Battleship - Create the X-matrix necessary to feed into a classifier
3.  Creating a Classifier - Create the ancilary data to make a classifier and explore the inner workings of what happens at each stage.
4.  Tuning the Classifier - Validation data ensures that the classifier should work well on real-world data. Find out what quality testing you can do to make the system more accurate. 

## Tools to Use

Feel free to use any of the following tools to aid you on your journey.   
[VSCode Download](https://code.visualstudio.com/download)
[Python Installation instructions](https://github.com/PackeTsar/Install-Python)  
[Bv-brc.org](https://www.bv-brc.org/)  
[Chatgpt](https://chat.openai.com/)  
[Git Bash Installation](https://git-scm.com/downloads)
**Note: Git is only required if you are planning on using the versioning functionality of the repository

## Getting Started

Follow these directions to download and use this course.

1. Login to github using the "sign in" button in the top right hand corner of the page. 

2. Then go to the git repository "https://github.com/vparrello/FIG-Bioinformatics-Course/tree/master" (which is this one) and look for the green "<> Code" button. It can be found on the righthand side of the page between the "Add File" button and the "About" section. Once you click on the button, you will find a "Download zip" button that contains all the files in the repository. Click on "Download Zip" to download the repository.

3. Find the file in your downloads folder and upzip it into a Project folder. It will create a folder in that space called "FIG-Bioinformatics-Course-master". Take note of where this folder is as you will need to navigate there in vscode.

4. If you have already downloaded VSCode, then skip to step 7. Otherwise, download vscode from "https://code.visualstudio.com/download" according to your operating system.

5. Find the file you just downloaded in your downloads folder and execute the file to start the installation.

6. Follow the instructions for the installation wizard until it goes to finish. All the default values used are acceptable for our purposes. If it does not automatically start VSCode when it finishes, search for it in your programs and apps to launch it.

7. When vscode launches, find "File" in the far left hand corner of the application. I'm talking all the way in that corner. Right next to the logo in windows or "Code" on Mac machines.

8. Choose Open Folder. Then find the "FIG-Bioinformatics-Course-master" folder you just put onto your computer and select it. This is the file from Step 3. It will prompt you to ask if you trust the authors. I would recommend "yes" because I trust myself, but that's up to you. ;-)

9. Happy Dance! You are done with your setup! Now go to your first exercise under 1 - Representative Genomes to get started.         
