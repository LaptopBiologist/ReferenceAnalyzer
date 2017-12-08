# ReferenceAnalyzer

This is a collection of Python functions for building and curating a set of repeat consensus sequences. 

ReferenceRepeatAnalyzer.py uses time-domain analysis to identify tandem repeats in a reference genome and generate consensus sequences. It includes function for computing autocorrelation and cross-correlation functions on DNA sequences, which are equivalent to ungapped sequence alignments.

FindHomology.py is primarily for identifying redundancies in a set of consensus sequences. It use BLASTn to perform all pairwise alignment of the consensus sequences and then constructs a graph from the alignments. Each sequence is represented as a vertice in the graph and sequences are connect with edges if they share an alignment with percent identity greater than a specified cutoff. The graph is then partitioned using the Louvain community finding algorithm. Each community reflects a set of sequences that share homology and should be manually examined. A text file (communities.txt) is output which describes the identified communities and the sequences in each community are output to separate \*.fa files. Those sequences that do not share homology with other sequences are output to singletons.fa.

Dependencies
Python-Louvain cannot be installed with pip. See for instructions: https://pypi.python.org/pypi/python-louvain/0.3
