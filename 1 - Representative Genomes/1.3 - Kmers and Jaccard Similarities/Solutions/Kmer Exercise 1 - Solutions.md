1. and 2. should hopefully be self-explanatory.
3. The number of Kmers in a sequence of length `n` is `(n−k+1)`.

4. The sequence 'cttgtatagtattgtcttgtgtatc' is 25 characters long; therefore, the number of 20-mers in this sequence is (25 - 20 + 1) = 6.

Notice the pattern in the set of 20-mers:
```
cttgtatagtattgtcttgt
ttgtatagtattgtcttgtg
tgtatagtattgtcttgtgt
gtatagtattgtcttgtgta
tatagtattgtcttgtgtat
atagtattgtcttgtgtatc
```
in each instance, the first character at the left of the preceding instance is removed, the remainder of the subsequence shifts over by one, and the next character from the original sequence is appended to the right end.

5. No, the Kmer 'atagctcga' does not occur within the sequence 'cttgtatagtattgtcttgtgtatc'.

6. Yes, the Kmer 'agtattgtc' does occur within the sequence 'cttgtatagtattgtcttgtgtatc', beginning at the 7th character.