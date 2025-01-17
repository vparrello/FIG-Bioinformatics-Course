# Command Line Foo

This will be your hub for basics of the command line that you will need to know as you traverse this course. Use this as your "cheat sheet" when working with the command line. If this is your first interaction with such a thing, make sure you keep this handy!

## Basics to know when working in a command line
1. There is no spell check in the command line; this means if there is anything misspelled, the terminal will not catch it and your command will likely error!
2. Commands are case-sensitive; your commands will error if the wrong case is used.
3. Much like spelling and capitals, spaces are also sensitive when using the command line. A missed space or an added space will likely mess up your command.
4. Within the bash shell, you only use the forward slash (/) and never the back slash.
5. Viewing files in VSCode in the editing mode can be hard to read. To change the viewing of the file, use the shortcut 'shift - command - v' on Macointosh and 'shift - control - v' on Windows.
6. If there is a command you don't understand, just ask Grimoire! Many times Grimoire will give you a very long explanation of a command, but the first couple paragraphs will be the most helpful.

## Basic commands to use
1. To list all files and directories
```
ls
```
2. To navigate between folders
```
cd [directory name]
```
3. To create a new directory
```
mkdir [directory name] (useful command) (move to)
rmdir [directory name]
```
4. To create a new file
```
touch [file name]
```
5. Move or rename a file or directory
```
mv [oldname.txt] [newname.txt]
```
6. Clear all past commands in terminal
```
clear
```
7. Send argument to standard output and generate collum header
```
p3 echo
```
8. 
```
echo [in appendix]
```
9. Print the working directory (directory you are already in)
```
pwd
```
10. Copy files or directories
```
cp file.txt /destination/
example
file to file and file to file in another directoryt
```
11. Remove a file
```
rm [file name]
```
12. Display file content
```
cat [filename.txt]
```
13. fhiewohqfg
```
head
```
tail
less
gzip compresses a file and gunzip uncompresses a file
gunzip/c filename.txt
df -h
wc - word count and line count, hammer exercise 4 quality control
diff



# Creating File Paths
Within this course, you will eventually be asked to create your own file paths, and this could cause some confusion. In this section, we will go over what file paths are and how they are created.

## What are file paths?
File paths function as a road map to a specific folder or file. In the same way directions on a map might show you landmarks on the way to your destination, a file path will lead you through each landmark in your file system.

File paths are structured like this:
```
home/user/Documents/file.txt
```
Each forward slash acts as a seperator between landmarks. These are the stops on the way to your destination. In this template, 'home' would be the first jump, 'user' the next jump, 'Documents' the third jump, and 'file.txt' being the destination. Let's see what a real file path looks like using the directory tree you are working with.

In order to get to the course syllabus, you have to drop down the menu of the course, then the 0_Getting-Started directory, then the course syllabus is right in there! This is how we reach it in the Explorer section of VSCode. Reaching it through the terminal is very similar.

To get somewhere, you have to know where you are starting. To know exactly where you are, you can use this command:
```
pwd
```
This will tell you where you are starting. Since this is where you are at, there is no need to include it in your file path; you will started with the next landmark. In the case of the syllabus, if we're starting in FIG-Bioinformatics-Course, we need to go down into 0_Getting-Started. The file path would look like this:
```
0_Getting-Started
```
That's it! Depending on your command, this is all you need! If you want to move down into the 0_Getting-Started directory, we use this command:
```
cd 0_Getting-Started
```
This will change your starting place to the directory of 0_Getting-Started. You'll notice this reflected at the beginning of your command line. Now you can freely access the course syllabus by using commands like 'cat' which shows the contents of a file.
Now, using this method of step by step is fool-proof, but tedious. The way to include all these landmarks in the same command is to use forward slashes between steps. In the step by step, our starting place was FIG-Bioinformatics-Course, then we went to 0_Getting-Started, and finally landed on our destination of 0_Course-Syllabus.md. This is all the information we need to complete a file path which looks like this:
```
0_Getting-Started/0_Course-Syllabus.md
```
Now you can stick this file path onto a command like this, which will show the contents of the syllabus in your terminal:
```
cat 0_Getting-Started/0_Course-Syllabus.md
```
Remember that you can move into other directories, but you cannot move into files. You can view files, create files, remove files, and many other commands!












(explain what a file path is? Difference between absolute and relative) is there any mention in the rest of the course about absolute and relative file paths? Do I need to go into this more? Create exercise for file paths; scavenger hunt??

command line to find file path
using cd and ls
using search feature for the gooey
(Files, folders, and paths) Name for the
Copy template into a file
cp command
  
  cd command
  mkdir command
  touch (file) command 
  look up touch command and file creation
  Mv command (google what command comes up first)
Mac: shift - command - v
Windows: shift - control - v

explain gitbash command ___________
