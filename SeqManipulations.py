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


import numpy
import scipy

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
def main():
    pass

if __name__ == '__main__':
    main()
