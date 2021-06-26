import numpy as np
import matplotlib.pyplot as plt
import sys
from Bio import SeqIO
from ushuffle import shuffle, Shuffler
from matplotlib import pyplot as plt
from statistics import variance
import math
from scipy import stats
import pylab




"""
 Is a simple i.i.d. randomization a good model for proteins? Choose a random set of 50 protein sequences (from ncbi or uniprot).
 This problem can be combined with project 3, sequence randomization, for a total of 30-80 pts.
 For the test data set, calculate the distribution of word scores using words of various lengths (minimum: 8, 12, 16) using a Blosum scoring table (or any table read from an external file).
 Display as a frequency table (score vs count).
o Bonus: 10 pts. Display as a histogram plot
 Compare the score distributions generated for random sequences with different k-mer word
lengths (n=1..5) preserved. Use the code in project 3, above, or the program Ushuffle.
 Mean and standard deviation are sufficient for the comparison, although you may include more
sophisticated tests.
 Bonus: 20 pts. repeat for DNA using word lengths 4, 6, and 8.
 What k-mer model best fits the unrelated protein data?

"""

def count_kmers(read, k):

    # Start with an empty dictionary
    counts = []
    # Calculate how many kmers of length k there are
    num_kmers = len(read) - k + 1
    # Loop over the kmer start positions
    for i in range(num_kmers):
        # Slice the string to get the kmer
        kmer = read[i:i+k]
        # Add the kmer to the dictionary if it's not there
        counts.append(kmer)
        # Increment the count for this kmer

    # Return the final counts
    return counts


input=sys.argv[1]
scoringmatrix = sys.argv[2]
flag=sys.argv[3]
seenpairs=[]

class Matrix:
  def __init__(self, matrix_filename):
    self._load_matrix(matrix_filename)

  def _load_matrix(self, matrix_filename):
    with open(matrix_filename) as matrix_file:
      matrix = matrix_file.read()
    lines = matrix.strip().split('\n')

    header = lines.pop(0)
    columns = header.split()
    matrix = {}

    for row in lines:
      entries = row.split()
      row_name = entries.pop(0)
      matrix[row_name] = {}

      if len(entries) != len(columns):
        raise Exception('Improper entry number in row')
      for column_name in columns:
        matrix[row_name][column_name] = entries.pop(0)

    self._matrix = matrix

  def lookup_score(self, a, b):
        a = a.upper()
        b = b.upper()

        if a not in self._matrix or b not in self._matrix[a]:
            print(a)
            print(b)
            raise InvalidPairException('[%s, %s]' % (a, b))
        return self._matrix[a][b]



if int(flag)==1:
    BlosumMatrix = Matrix(scoringmatrix)
    BlosumMatrix._load_matrix(scoringmatrix)

    kmersizelist=[8,12,16]
    finaldict={}
    for z in kmersizelist:
        allseq=[]
        for seq_record in SeqIO.parse(input, "fasta"):
          seq=seq_record.seq
          kmerdict=count_kmers(str(seq),z)
          allseq.append(kmerdict)


        flat_list = [item for sublist in allseq for item in sublist]
        #print(flat_list)
        #print(len(flat_list))

        for i in range(len(flat_list)):
            for j in range(len(flat_list)):
                if(i>=j):
                    continue
                else:
                    score = 0
                    for a in range(len(flat_list[j])):

                        A=flat_list[i][a]

                        B=flat_list[j][a]

                        value = float(BlosumMatrix.lookup_score(A, B))
                        score += value

                    if score not in finaldict:
                        finaldict[score] = 0
                        # Increment the count for this kmer
                    finaldict[score] += 1

        print("K="+str(z))
        print(finaldict)

        sumtot = sum(k * v for k, v in finaldict.items())
        totalcounts = sum(finaldict.values())
        av=(sumtot / totalcounts)
        print("Average Score:"+str((av)))

        plt.bar(finaldict.keys(), finaldict.values(), color='pink')
        plt.axvline(sumtot / totalcounts, color='teal', linestyle='dashed', linewidth=1)
        plt.xlabel('Score')
        plt.title("Histogram for K="+str(z))
        plt.ylabel('Count')
        plt.show()

        data = [i for item, count in zip(finaldict.keys(), finaldict.values()) for i in [item] * count]
        var=variance(data)
        print("Variance in Score:" + str(var))
        print("Standard Deviation in Score:" + str(var**0.5))

        plt.boxplot(data)
        plt.title("Boxplot for K=" + str(z))
        plt.show()

        stats.probplot(data, dist="norm", plot=pylab)
        pylab.title("QQPlot for K=" + str(z))
        pylab.show()

        kstest_test = stats.kstest(data,'norm')
        print(kstest_test)




    kmersizelist = [1,2,3,4,5]
    finaldict = {}
    for z in kmersizelist:
        allseq = []
        for seq_record in SeqIO.parse(input, "fasta"):
            seq = seq_record.seq
            seq=seq.encode("utf-8")
            seq=shuffle(seq, 2).decode("utf-8")
            kmerdict = count_kmers(str(seq), z)
                # print(kmerdict)
            allseq.append(kmerdict)



        flat_list = [item for sublist in allseq for item in sublist]
        #print(flat_list)
        #print(len(flat_list))

        for i in range(len(flat_list)):
            for j in range(len(flat_list)):
                if(i>=j):
                    continue
                else:
                    score = 0
                    for a in range(len(flat_list[j])):

                        value = float(BlosumMatrix.lookup_score(flat_list[i][a], flat_list[j][a]))
                        score += value

                    if score not in finaldict:
                        finaldict[score] = 0
                        # Increment the count for this kmer
                    finaldict[score] += 1

        print("K="+str(z))
        print(finaldict)

        sumtot = sum(k * v for k, v in finaldict.items())
        totalcounts = sum(finaldict.values())
        av = (sumtot / totalcounts)
        print("Average Score:" + str((av)))

        plt.bar(finaldict.keys(), finaldict.values(), color='blue')
        plt.axvline(sumtot / totalcounts, color='red', linestyle='dashed', linewidth=1)
        plt.xlabel('Score')
        plt.title("Histogram for K="+str(z))
        plt.ylabel('Count')
        plt.show()

        data = [i for item, count in zip(finaldict.keys(), finaldict.values()) for i in [item] * count]
        var=variance(data)
        print("Variance in Score:" + str(var))
        print("Standard Deviation in Score:" + str(var**0.5))

        plt.boxplot(data)
        plt.title("Boxplot for K=" + str(z))
        plt.show()

        stats.probplot(data, dist="norm", plot=pylab)
        pylab.title("QQPlot for K=" + str(z))
        pylab.show()

        kstest_test = stats.kstest(data,'norm')
        print(kstest_test)



if int(flag) == 2:
    kmersizelist = [4,6,8]
    finaldict = {}
    for z in kmersizelist:
        allseq = []
        for seq_record in SeqIO.parse(input, "fasta"):
            seq = seq_record.seq
            kmerdict = count_kmers(str(seq), z)
            allseq.append(kmerdict)

        flat_list = [item for sublist in allseq for item in sublist]
        # print(flat_list)
        # print(len(flat_list))

        for i in range(len(flat_list)):
            for j in range(len(flat_list)):
                if (i >= j):
                    continue
                else:
                    score = 0
                    for a in range(len(flat_list[j])):
                        A = flat_list[i][a]

                        B = flat_list[j][a]

                        if( A==B):
                            score=score+1

                    if score not in finaldict:
                        finaldict[score] = 0
                        # Increment the count for this kmer
                    finaldict[score] += 1

        print("K=" + str(z))
        print(finaldict)

        sumtot = sum(k * v for k, v in finaldict.items())
        totalcounts = sum(finaldict.values())
        av = (sumtot / totalcounts)
        print("Average Score:" + str((av)))
        plt.bar(finaldict.keys(), finaldict.values(), color='pink')
        plt.axvline(sumtot / totalcounts, color='teal', linestyle='dashed', linewidth=1)
        plt.xlabel('Score')
        plt.title("Histogram for K="+str(z))
        plt.ylabel('Count')
        plt.show()

        data = [i for item, count in zip(finaldict.keys(), finaldict.values()) for i in [item] * count]
        var=variance(data)
        print("Variance in Score:" + str(var))
        print("Standard Deviation in Score:" + str(var**0.5))

        plt.boxplot(data)
        plt.title("Boxplot for K=" + str(z))
        plt.show()

        stats.probplot(data, dist="norm", plot=pylab)
        pylab.title("QQPlot for K=" + str(z))
        pylab.show()

        kstest_test = stats.kstest(data,'norm')
        print(kstest_test)

    kmersizelist = [1,2,3]
    finaldict = {}
    for z in kmersizelist:
        allseq = []
        for seq_record in SeqIO.parse(input, "fasta"):
            seq = seq_record.seq
            seq = seq.encode("utf-8")
            seq = shuffle(seq, 2).decode("utf-8")
            kmerdict = count_kmers(str(seq), z)
            # print(kmerdict)
            allseq.append(kmerdict)

        flat_list = [item for sublist in allseq for item in sublist]
        # print(flat_list)
        # print(len(flat_list))

        for i in range(len(flat_list)):
            for j in range(len(flat_list)):
                if (i >= j):
                    continue
                else:
                    score = 0
                    for a in range(len(flat_list[j])):
                        A = flat_list[i][a]

                        B = flat_list[j][a]

                        if( A==B):
                            score=score+1
                            
                    if score not in finaldict:
                        finaldict[score] = 0
                        # Increment the count for this kmer
                    finaldict[score] += 1

        print("K=" + str(z))
        print(finaldict)

        sumtot = sum(k * v for k, v in finaldict.items())
        totalcounts = sum(finaldict.values())
        av = (sumtot / totalcounts)
        print("Average Score:" + str((av)))

        plt.bar(finaldict.keys(), finaldict.values(), color='blue')
        plt.axvline(sumtot / totalcounts, color='red', linestyle='dashed', linewidth=1)
        plt.xlabel('Score')
        plt.title("Histogram for K="+str(z))
        plt.ylabel('Count')
        plt.show()

        data = [i for item, count in zip(finaldict.keys(), finaldict.values()) for i in [item] * count]
        var=variance(data)
        print("Variance in Score:" + str(var))
        print("Standard Deviation in Score:" + str(var**0.5))

        plt.boxplot(data)
        plt.title("Boxplot for K=" + str(z))
        plt.show()

        stats.probplot(data, dist="norm", plot=pylab)
        pylab.title("QQPlot for K=" + str(z))
        pylab.show()

        kstest_test = stats.kstest(data,'norm')
        print(kstest_test)

        # print(data)


