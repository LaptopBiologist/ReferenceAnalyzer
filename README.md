# ReferenceAnalyzer

!!!! STILL IN DEVELOPMENT !!!!

This is a collection of Python functions for building and curating a set of repeat consensus sequences. 

ReferenceRepeatAnalyzer.py uses time-domain analysis to identify tandem repeats in a reference genome and generate consensus sequences. It includes function for computing autocorrelation and cross-correlation functions on DNA sequences, which are equivalent to ungapped sequence alignments.

FindHomology.py is primarily for identifying redundancies in a set of consensus sequences. It uses BLASTn to perform all pairwise alignment of the consensus sequences and then constructs a graph from the alignments. Each sequence is represented as a vertice in the graph and sequences are connect with edges if they share an alignment with percent identity greater than a specified cutoff. The graph is then partitioned using the Louvain community finding algorithm. Each community reflects a set of sequences that share homology and should be manually examined. A text file (communities.txt) is output which describes the identified communities and the sequences in each community are output to separate \*.fa files. Those sequences that do not share homology with other sequences are output to singletons.fa.

Usage:


python FindHomology.py -i <path to input fasta> -o <output directory to be created> -B <path to the BLAST directory> -C <percent identity cutoff>
  
For example:

python FindHomology.py -i c:/.../Repbase_19.06_DM_7-8-16.fa -o C:/.../test_communities -B "c:\NCBI\blast-2.5.0+" -C 80

Note: It is important that the directory containing the input file does not contain a subdirectory named "/temp", as the function creates and subsequently deletes a temp folder.
 

Dependencies

Statistics and scientific computing: Numpy, Scipy, Statsmodels, Scikit-Learn

Graph manipulations: Networkx, Python-Louvain

Bioinformatics: Biopython,  HTSeq 0.6.1

Plotting: Matplotlib, Seaborn

Most of these dependencies are included with the common Python bundles (eg Enthought Canopy). If you don't already have them installed, most are easy to install with pip.  The exceptions are HTSeq and Python-Louvain. On Windows the current version of HTSeq (0.9.1) is difficult to install, but version 0.6.1 can be installed with an executable found at https://pypi.python.org/pypi/HTSeq/0.6.1. Python-Louvain cannot be installed with pip. For simple installation instructions, see: https://pypi.python.org/pypi/python-louvain/0.3
