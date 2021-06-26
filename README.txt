Approach-

For this task I used 20 random protein sequences from UniProts random protein file( <50%).

My program starts with a k, finds all the kmers per sequence, does a pairwise matching based on BLOSUM62. I maintain a count for every score using a dictionary. This dictionary is printed. I use this dictionary to draw a Histogram. I also draw Boxplots, QQplots, and conduct the Kolmogorov-Smirnov Test for normality.

I repeat this for all k values.

Next, I try shuffling the sequences using UShuffle. and conduct the same procedure once again.



I back translate these proteins to DNA sequences and repeat the same analysis. Match/Mismatch scheme is 1/0.

Algorithmic Analysis-

Say x is the average length per sequence in the fast file. n is the no. of sequences and k is the kmer number.

No. of kmers for 1 sequnece=x-k+1

No. of kmers for n sequences=n(x-k+1)

pairwise comparison of n(x-k+1) kmers = n(x-k+1)C2

Also, we will be iterating through each of these pairs.



Time Complexity =O (kn(x-k+1)C2)= 

This is cubic in k, meaning as K increases, the time complexity increases at a cubic rate.

Protein:

How to run my program?

python3 randomization.py randomp.fasta blosum62.txt 1

randomization.py is the python file, randomp.fasta contains the protein sequences, and 1 is the flag for proteins.


DNA:

How to run my program?

python3 randomization.py randomd.fasta blosum62.txt 2

randomization.py is the python file, randomd.fasta contains the nucleotide sequences, and 2 is the flag for nucleotide.

NOTE: the program takes in blosum.txt as a parameter but does NOT use it.
