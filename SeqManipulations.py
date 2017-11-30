#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      I am


#
# Created:     02/11/2017
# Copyright:   (c) I am 2017
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import time
import numpy
import scipy
import sklearn
IUPAC_DICT={'W':['A','T'], 'S':['C','G'], 'M':['A','C'],\
 'K':['G','T'], 'R':['A','G'], 'Y':['C','T'],\
 'B':['C','G','T'], 'D':['A','G','T'], 'H':['A','C', 'T'],\
 'V':['A','C','G'], 'N':['A','G','C','T'] }

def ReplaceAmbiguousNucleotides(sequence):
    """Replaces ambiguous nucleotides with one of the nt that could be represented."""
    seq_array=numpy.fromstring(sequence.upper(), '|S1')
    for ambiguous_nt in IUPAC_DICT.keys():
        amb_ind=numpy.where(seq_array==ambiguous_nt)[0]
        replacement_nt=numpy.random.choice(IUPAC_DICT[ambiguous_nt], size=len(amb_ind), replace=True)
        seq_array[amb_ind]=replacement_nt
    return ''.join(seq_array)

def ParseSeqName(seq_name):
    level_1=seq_name.split("_")
    info_dict={}
    for field in level_1:
        info=field.split('=')
        if len(info)==2:
            info_dict[info[0]]=info[1]
    return info_dict

def CountKMERS(sequence, k=10 ):
    """Decomposes sequence into kmers."""
    kmerSet={}
    for x in range(0,len(sequence)-k+1):
        kmer=str(sequence[x:x+k])
        if kmerSet.has_key(kmer)==False:
            kmerSet[kmer]=0.
        kmerSet[kmer]+=1.
    return kmerSet
##        kmerSet.add(str(complement[x:x+k]))
    return list(kmerSet)

def SequenceToInt(seq):
    seq_array=numpy.fromstring(seq, '|S1')
    int_array=numpy.ndarray((len(seq_array),))
    int_array.fill(0)
    nt_array=['T', 'C', 'G']
    for i, nt in enumerate( nt_array):
        int_array[seq_array==nt]=i+1
    return int_array

def ComputeNeighborDistance(sequence,k=20, mismatches=1, max_sample=200):
    kmer_dict=CountKMERS(sequence,k)
    kmer_list=kmer_dict.keys()
    radius=float(mismatches)/len(kmer_list[0])
    start=time.clock()
    encoded_kmers=numpy.array( [SequenceToInt(kmer) for kmer in kmer_list ])
    ball_tree=sklearn.neighbors.BallTree(encoded_kmers, metric='hamming')
##    print "Build time = {0}".format(time.clock()-start)
    count_list=[]
    start=time.clock()
    indices=numpy.arange(len(encoded_kmers))
    if len(encoded_kmers)>max_sample:
        indices=[True]*max_sample+[False]*(len(encoded_kmers)-max_sample)
        numpy.random.shuffle(indices)
    for kmer in encoded_kmers[indices]:
        neighbors=ball_tree.query_radius(kmer[None,:], radius)[0]

        count_list.append(sum( [kmer_dict[kmer_list[n] ] for n in neighbors]))
##    print"Check time = {0}".format(time.clock()-start)
    S= numpy.mean( count_list)

    return S

def CorrectedNeighborDistance(sequence, k=20, mismatches=3, step=200):
    s1=ComputeNeighborDistance(sequence, k, mismatches)
    s2=ComputeNeighborDistance(sequence, k+step, mismatches)
    linreg=scipy.stats.linregress([k, k+step], [s1, s2])
    slope=linreg.slope
    intercept=linreg.intercept
    est_rpt_len=(-1* slope)**-1
    score=intercept
    return score, est_rpt_len



def main():
    pass

if __name__ == '__main__':
    main()
