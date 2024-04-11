#### Kmer Exercise 1 - What are Kmers?

Objective: Learn about the concept of "Kmers" in sequences.

#### Materials: 

[Grimoire](https://chat.openai.com/g/g-n7Rs0IK86-grimoire)
``` 
FIG-Bioinformatics-Course/
    Basics.html
    1 - Representative Genomes/
        1.3 - Kmers and Jaccard Similarities/
            Kmer Exercise 1 - What are Kmers.md
            Kmer Exercise 1 - Solutions.md
```

#### Exercise:

The term `Kmers` is short for "Subsequence of length K", where `K` stands for some integer. We will be interested in Kmers from both DNA and protein sequences.
Typically, we will be working with DNA 20-mers and protein 8-mers.

1. Enter the following into Grimoire's "Message Grimoire" box:
```
  I am going to give you some definitions in an attached file;
  you don't need to respond to them, just learn them.
  I'm then going to ask you some questions.
```
Then click on the "Paperclip" icon, select the file `Basics.html`, and hit the "Send Message" icon (the "Up-arrow" icon at the right of the "Message Grimoire" box).
Grimoire should respond to entering the above with something like "Understood" or "Got it". `Basics.html` contains a list of definitions of terms and concepts specific to this course. We will be using the preceding method to "pre-prompt" Grimoire regarding the "context" for a set of exercises throughout the remainder of this course; you are welcome to look at the file `Basics.html`, but some of the terms, definitions, and concepts may not make sense in advance of the exercise they are relevant to.

The following questions should be entered as individual prompts (or as Grimore puts them, "Messages"), but should be entered during the same session in which you uploaded the above definitions, since Grimoire will need these definitions to establish its "context", and it cannot remember the contents of one session in a different dession, but will instead fall back on its "default" definitions unless "re-prompted".

2. Ask Grimoire "What is a Kmer, and how are Kmers used in bioinformatics?"

3. "How can I know how many Kmers are in a given sequence?"
(Grimoire should answer with an explanation or formula that calculates the number of Kmers given the length of the sequence and `K`)

4. "Let 'cttgtatagtattgtcttgtgtatc' be a DNA sequence; please tell me the length of this sequence, and then list all of the 20-mers within this sequence."
(NOTE: Do please check how many characters are in Grimoire's answers, because sometimes it will miscount the number of characters in a sequence -- e.g., it may get the length of the sequence wrong, or return a list of "19-mers" instead of "20-mers". Grimoire is not perfect! But it's pretty good at what it does do well, which is write code and explain concepts.)

5. "Is the Kmer 'atagctcga' within the sequence 'cttgtatagtattgtcttgtgtatc'?"

6. "Is the Kmer 'agtattgtc' within the sequence 'cttgtatagtattgtcttgtgtatc'?"
