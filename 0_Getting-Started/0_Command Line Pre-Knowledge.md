# Command Line Pre-Knowledge

This will be your hub for basics of the command line that you will need to know as you traverse this course. Use this as your "cheat sheet" when working with the command line. If this is your first interaction with such a thing, make sure you keep this handy!

## Basics to know when working in a command line
1. There is no spell check in the command line; this means if there is anything misspelled, the terminal will not catch it and your command will likely thrown an error!
2. Commands are case-sensitive; your commands will throw an error if the wrong case is used.
3. Much like spelling and capitals, spaces are also sensitive when using the command line. A missed space or an added space will likely mess up your command.
4. Within the bash shell, you only use the forward slash (/) and never the back slash.
5. Viewing files in VSCode in the editing mode can be hard to read. To change the viewing of the file, use the shortcut 'shift + command + v' on Mac and 'shift + control + v' on Windows.
6. If there is a command you don't understand, just ask Grimoire! Many times Grimoire will give you a very long explanation of a command, but the first couple paragraphs will be the most helpful.
7. Any command that does not have a file path (we will explain what a file path is later in this document) will assume that you are already in the directory you need to be in. If the file or files you need for your command are in a different directory, your command will need a file path to get there.

# Basic commands to use
While practice is important, some examples of these commands contain real files and paths. If you want to practice commands, create your own testing files; this is to keep you from removing or losing important files. 

## File and Directory Commands
1. List files and directories in the current directory
```
ls
```
List all files in a specific directory (we will explain what we mean by 'path' later in this document)
```
ls path/to/directory
```
*Example use of this command:*
```
ls 1_Representative-Genomes/1.1_Tab-Separated-Value_TSV_Files/Solutions
```
List files and directories in the current directory, including hidden files
```
ls -la
```
2. Navigate to the named folder
```
cd path/to/directory
```
*Example:*
```
cd FIG-Bioinformatics-Course/0_Getting-Started
```

3. Show full path to the current (working) directory
```
pwd
```

4. Move files to a different directory
```
mv file1 file2 file4 directory_destination
```
*Example:*
```
mv 0_Course-Syllabus.md 0_Download-Instructions.md FIG-Bioinformatics-Course
```
5. Rename a file or directory
```
mv oldname newname
```
*Example:*
```
mv data.tbl dataNew.tbl
``` 
6. Clear all text in the terminal
```
clear
```
7. Copy a file
```
cp filename duplicatefilename
```
*Example:*
```
cp data.tbl datacopy.tbl
```
To copy a file to another directory
```
cp path/to/directory
```
*Example:*
```
cp FIG-Bioinformatics-Course/Data/data.tbl
```

8. Delete (remove) a file
```
rm filename
```
*Example:*
```
rm data.tbl
```

9. List all commands entered during a session
```
history
```

## File Viewing and Manipulation
1. Display the contents of a file
```
cat filename
```
*Example:*
```
cat data.tbl
```
2. View file content one screen at a time
```
less filename
```
*Example:*
```
less data.tbl
```
3. View the first 10 lines of a file
```
head -n 10 filename
```
*Example:*
```
head -n 10 data.tbl
```

4. View the last 10 line of a file
```
tail -n 10 filename
```
*Example:*
```
tail -n 10 data.tbl
```

## Disk and System Monitoring
1. Display how much diskspace is in use
```
df -h
```
2. Display the size of a file or directory
```
df -sh path/or/directory
```
*Example:*
```
df -sh FIG-Bioinformatics-Course/0_Getting-Started
```


## Archiving and Compression
1. Compress a file and convert to a '.gz' file
```
gzip filename
```
*Example:*
```
gzip data.tbl
```

2. Decompress a file and convert to a different kind of file
```
gunzip filename.gz
```
*Example:*
```
gunzip data.gz
```

3. Decompress a STDOUT
```
gunzip -c filename.gz
```
*Example:*
```
gunzip -c data.gz
```

## Advanced Commands
1. Make a new directory
```
mkdir directory_name
```
*Example:*
```
mkdir FIG-Bioinformatics-Course
```

2. Remove an empty directory (directories containing files cannoy be removed until emptied)
```
rmdir directory_name
```
*Example:*
```
rmdir FIG-Bioinformatics-Course
```

# Creating File Paths
Within this course, you will eventually be asked to create your own file paths, and this could cause some confusion. In this section, we will go over what file paths are and how they are created.

## What are file paths
File paths function as a road map to a specific folder or file. In the same way directions on a map might show you landmarks on the way to your destination, a file path will lead you through each landmark in your file system.

File paths are structured like this:
```
home/user/Documents/file
```
Each forward slash acts as a seperator between landmarks. These are the stops on the way to your destination. In this template, 'home' would be the first jump, 'user' the next jump, 'Documents' the third jump, and 'file' being the destination. Let's see what a real file path looks like using the directory tree you are working with.

## How to use file paths
In order to get to the course syllabus, you have to drop down the menu of the course, then the 0_Getting-Started directory, then the course syllabus is right in there! This is how we reach it in the Explorer section of VSCode. Reaching it through the terminal is very similar.

*In this next section, follow along with the commands in your own terminal!*

To get somewhere, you have to know where you are starting. To know exactly where you are, you can use this command:
```
pwd
```
This will tell you where you are starting. Since this is where you are at, there is no need to include it in your file path; you will start with the next landmark. In the case of the syllabus, if we're starting in FIG-Bioinformatics-Course, we need to go down into 0_Getting-Started. The map part of your file path would look like this:
```
0_Getting-Started
```
That's it! Depending on your command, this is all you need! If you want to move down into the 0_Getting-Started directory, we use this command:
```
cd 0_Getting-Started
```
This will change your starting place to the directory of 0_Getting-Started. You'll notice this reflected at the beginning of your command line. Now you can freely access the course syllabus (or any other file in this directory) by using commands like 'cat' which shows the contents of a file.

Now, using this method of step by step is fool-proof, but tedious. The way to include all these landmarks in the same command is to use forward slashes between jumps. In the step by step, our starting place was FIG-Bioinformatics-Course, then we went to 0_Getting-Started, and finally landed on our destination of 0_Course-Syllabus.md. This is all the information we need to complete a file path which looks like this:
```
0_Getting-Started/0_Course-Syllabus.md
```
Now you can stick this file path onto a command like this, which will show the contents of the syllabus in your terminal:
```
cat 0_Getting-Started/0_Course-Syllabus.md
```
Remember that you can move into other directories, but you cannot move into files. You can view files, create files, remove files, and many other commands!

If you've run into any error while following along in your terminal, retrace your steps; if you still don't know what went wrong, ask Grimoire!
