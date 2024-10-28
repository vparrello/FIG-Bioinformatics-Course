These are instructions on how to download the course through github. If you have not done so already, follow the instructions video for your computer's operating system and follow along. 

1. Download the required programs for this course. There are 4 in total. Follow the links below to get your computer's operating system's version of the necessary programs.

    * [VSCode Download](https://code.visualstudio.com/download) - [Setup Video](https://code.visualstudio.com/docs/setup/setup-overview)
    * [Python Installation instructions](https://github.com/PackeTsar/Install-Python) 
    * [Git Bash Installation](https://git-scm.com/downloads) 

2. Once you have downloaded these programs onto your computer, unzip one of the following videos to guide you through the installation of each program. Default settings are recommended for all installations.

    * Windows11-initial-setup-v0_1.zip
    * MacOS-initial-setup-v0_1.zip (Comming Soon)
    * Linux-initial-setup-v0_1.zip (Coming Soon)

3. Launch VSCode. If you have not already chosen a folder, look in the "File" menu at the top left of VSCode and click on "Open Folder". Then use the prompts to choose an empty folder.

4. Open a Terminal inside VSCode by going to the menu in the top left corner and following File, Edit, Selection, and so on until you find Terminal. Click on it and select New Terminal.

5. Inside the terminal in the upper right hand corner, there is an icon of a + and a dropdown menu. Click inside the dropdown menu and select "Select Default Profile". A prompt will show up at the top of your screen asking which profile you wish to choose. Click on GitBash.
    **Note: If GitBash is not an option at this point, please review the Git Bash Installation and set up video instructions. Your Git Bash installation may be misconfigured

6. Click in the new terminal and type `git --version` to verify your git bash installation.

7. Next to download the course material, type `git clone https://github.com/vparrello/FIG-Bioinformatics-Course.git` inside your terminal.

8. Once that finishes, type `git branch ` and then your name. This will be the name of your branch.  

9. Next type `git checkout ` and then the name of your branch. This creates a personal copy of the material for you.

git branch --set-upstream-to=origin/master <branch>

9. Happy Dance! You are done with your setup! The course syllabus and how to update the course material are located in the folder `0_Getting-Started` along with this document. Go to exercise `TSV-Exercise-1_Learning-to-use-Grimoire` under `1_Representative-Genomes/1_Tab-Separated-Value_(TSV)_Files` to get started on the course.