#-------------------------------------------------------------------------------
# Name:        Reference Repeat Analyzer
# Purpose:
#/
# Author:      I am
#
# Created:     07/09/2017
# Copyright:   (c) I am 2017
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy
import scipy
from Bio import SeqIO
from Bio import Seq
from Bio.Blast import NCBIXML

import matplotlib
##matplotlib.use('Agg')
from matplotlib import pyplot
from statsmodels.tsa import stattools
from statsmodels.tsa.stattools import acf
import sklearn
import sklearn.gaussian_process
import sklearn.kernel_ridge
import sklearn.metrics
import SeqManipulations
import sklearn.gaussian_process.kernels as kernels
from sklearn.base import BaseEstimator
from Bio import Seq
from matplotlib import pyplot
import seaborn
import itertools
import subprocess
import csv
import os
import collections
##import ReferenceAnalyzer.FindHomology
from FindHomology import ComputeKmerCompositionEntropy

import gzip

import sys
import shutil

import copy

pyplot.interactive(False)
#Global
MUSCLE_PATH=None
BLAST_PATH='c:/ncbi/blast-2.5.0+/bin/blastn.exe'
import matplotlib
seaborn.set_style('white')

class ReferenceRepeat():
    def __init__(self, sequence, left, right, tandem_flag):
        #The repeat's sequence as a string
        self.sequence=sequence

        #The repeat's left and right edges
        self.left=left
        self.right=right
        self.chrom=''
        self.length=right-left

        #Whether the repeat was found adjacent to another copy
        self.tandem_flag=tandem_flag

class RollingKernelAverage(BaseEstimator):
    def __init__(self,bandwidth=.5, weights=[]):
        self.bandwidth=bandwidth
        self.weights=weights
    def fit(self, X,Y):
        self.X=X
        self.Y=Y
        #build ball_tree
        self.KDTree=sklearn.neighbors.KDTree(self.X)

##        for x in X:
##            self.weights=Epanichenikov ( numpy.subtract.outer( self.X,self.X)/self.bandwidth)
##            weights/=weights.sum(0)


    def predict(self, X):
        mean_list=[]
        #find points within radius
        indices, distances=self.KDTree.query_radius(X, self.bandwidth, return_distance=True)
        for i in range(len( indices)):
##            print len(distances[i])
            weights=Epanichenikov (distances[i]  /self.bandwidth)
            weights/=weights.sum()

##        print weights.shape
            mean_list.append( (self.Y[indices[i]]*weights).sum(0))

        return numpy.array(mean_list)

    def score(self, X, Y):
        mean=self.predict(X)
        return -1*((mean-Y)**2).sum()




def MakeDir(newdir):
    if os.path.exists(newdir)==False:
        os.mkdir(newdir)

def splitFile(infile, tempDir, limit=500):
    MakeDir(tempDir)
    refSeq=GetSeq(infile)
    inDir='/'.join(infile.split('/')[:-1])
    root=infile.split('/')[-1]
    rootname='.'.join(root.split('.')[:-1])
    ext=root.split('.')[-1]
    fileCount=1
    outfile="{0}/{1}_temp_{2}.{3}".format(tempDir, rootname, fileCount, ext)
    outHandle=open(outfile, 'w')
    count=0
    print len(refSeq)
    for key in refSeq:
        if count>=limit:
            print count
            outHandle.close()
            fileCount+=1
            count=0
            outfile="{0}/{1}_temp_{2}.{3}".format(tempDir, rootname, fileCount, ext)
            print outfile
            outHandle=open(outfile, 'w')
        seq=str( refSeq[key])
        name= CleanName(key)

        outHandle.write('>'+name+'\n')
        outHandle.write(seq+'\n')
        count+=1

    outHandle.close()

def CountSequences(k=6, alphabet='ATGC'):
    """Enumerate all possible k-length products of the alphabet and then
    removes redundancies. Returns the cardinality."""
    sequences=itertools.product(alphabet, repeat=k)
    return list(sequences)


def OccupancyProbability(n, i, k):
    n=float(n)
    i=float(i)
    k=float(k)
    j=n-i
    n_j=n-j
    m=numpy.arange(n_j)
    coef_1=scipy.special.binom(n, j)
    coef_2=(-1)**(m)* scipy.special.binom(n-j, m)*(1.-(j+m)/n)**k
    return coef_1*coef_2.sum()


def Epanichenikov(u):
    k=.75*(1-u**2)*(abs(u)<=1)
    return k

def GetKMERS(sequence, k=10 ):
##    complement=str(Seq.Seq( sequence).reverse_complement())
    kmerSet=set()
    for x in range(len(sequence)-k+1):
        kmerSet.add(str(sequence[x:x+k]))
##        kmerSet.add(str(complement[x:x+k]))
    return list(kmerSet)




def GetSeq(ref, upper=False, rename=False, clean=True):

    """Reads a fasta, returns of a dictionary of strings keyed by entry name."""
    if ref.split('.')[-1]=='gz':
        handle=gzip.open(ref )
    else:
        handle=open(ref, 'r')
    lib=SeqIO.parse(handle, 'fasta')
    SeqLen={}
    for rec in lib:
##        if refName.count(rec.name)==0: continue
        if rename==True:
            name=rec.description.split(' ')[1].split('=')[1]
            print name
            SeqLen[CleanName(name)]=str( rec.seq)
        else:
            if clean==True:
                SeqLen[CleanName(rec.name)]=str( rec.seq)
            else:
                SeqLen[rec.name]=str( rec.seq)
        if upper==True: SeqLen[CleanName(rec.name)]=SeqLen[CleanName(rec.name)].upper()
    handle.close()
    return SeqLen

def CleanName(name):
    illegal=['|', '!','>', '<', '?','/','*', ':']
    for i in illegal:
        name='_'.join(name.split(i))
    return(name)

def TandemFinder(infile, outdir,muscle_path, threshold):
    SetMUSCLEPath(muscle_path)
    MakeDir(outdir)
    out_fasta='{0}/consensus_repeats.fa'.format(outdir)
    fasta_handle=open(out_fasta, 'w')
    out_annotation='{0}/genome_annotation.tsv'.format(outdir)
    annotation_handle=open(out_annotation, 'w')
    annotation_table=csv.writer(annotation_handle, delimiter='\t')
    header=['Chrom', 'Repeat Length', 'Copy number', 'Shannon information', 'Diversity', 'Complexity', 'Start', 'End','Per nucleotide information', 'Sequence']
    annotation_table.writerow(header)
    log_file='{0}/errlog.txt'.format(outdir)

    repeat_dir='{0}/repeats'.format(outdir)
    image_dir='{0}/images'.format(outdir)
    MakeDir(repeat_dir)
    MakeDir(image_dir)
    sequences=GetSeq(infile, rename=True)
    for key in sorted( sequences.keys(), reverse=True):
##        if key!='3R': continue
##        out_dir='/'.join(outfile.split('/')[:-1])
##        out_root='.'.join( outfile.split('/')[-1].split('.')[:-1])
##        outimage='{0}/{1}_{2}.png'.format(out_dir, out_root, CleanName( key))
        masked_seq=sequences[key]
        interval_dictionary= FindPeriodicity(sequences[key], '', key)
        #Sort intervals by length:
        interval_list=[]
        true_intervals={}
        for period in sorted( interval_dictionary.keys()):

            for interval in interval_dictionary[period]:
                l,r=interval
                interval_list.append([period, interval,abs( r-l)])
        sorted_intervals=sorted(interval_list, key=lambda x:x[2], reverse=True)

        for period, interval, interval_length in sorted_intervals:
            left, right=interval
            #Modified to return a dictionary
            if period<30:
                masked=MaskRepeatsOfSize(masked_seq[left:right], period, threshold)
                masked_seq=masked_seq[:left]+masked+masked_seq[right:]
                continue
            repeat_dict, masked=GetRepeatsInInterval(sequences[key][left:right], period, threshold, masked_seq[left:right])
            masked_seq=masked_seq[:left]+masked+masked_seq[right:]
            #Update repeat positions:
            repeat_count=0
            left_boundary, right_boundary=numpy.inf, 0
            for major_period in repeat_dict.keys():
                if true_intervals.has_key(major_period)==False:
                    true_intervals[major_period]=[]
                repeat_list=repeat_dict[major_period]
                left_boundaries, right_boundaries=[],[]
##                length_dict=
                for i in range(len(repeat_list)):
                    repeat_list[i].left+=left
                    repeat_list[i].right+=left
                    repeat_list[i].chrom=key

                    #Determine the correct boundaries of the array

                    if repeat_list[i].tandem_flag==True:
                        repeat_count+=1
                        left_boundaries.append( repeat_list[i].left)
                        right_boundaries.append( repeat_list[i].right)
                        if true_intervals.has_key(len(repeat_list[i].sequence))==False:
                            true_intervals[ len(repeat_list[i].sequence)]=[]
                        true_intervals[len(repeat_list[i].sequence)].append((repeat_list[i].left, repeat_list[i].right))
##                left_boundaries=numpy.array(sorted( list(set( left_boundaries))))
####                right_boundaries=numpy.array(right_boundaries)[sort_ind]
##                boundary_dist=numpy.diff(left_boundaries)
##                array_breaks=[0]+ list(numpy.where(boundary_dist>major_period*5)[0]+1)+[len(left_boundaries)-1]
##                if major_period>5000: print left_boundaries
####                print array_breaks
##                array_list=[]
##                for i in range(len( array_breaks)-1):
##                    left_edge=left_boundaries[ array_breaks[i]]
##                    right_edge=left_boundaries[ array_breaks[i+1]]+major_period
##
##                    array_list.append((left_edge, right_edge))
##                    true_intervals[major_period].append((left_edge, right_edge))
##                print array_list
                left_boundary=min(left_boundaries )
                right_boundary=max(right_boundaries)
                print left_boundary, right_boundary
##                name_root="Chr={0}_period={1}".format(key, period)
                if repeat_count==0: continue

                #Build consensus
                fasta_name='{0}/{1}_{2}_{3}_{4}'.format(repeat_dir, key.split('_')[0], major_period, left_boundary, right_boundary)
                cons_dict=SummarizeRepeats(repeat_list, fasta_name, log_file)

            #Write consensus
                for cons_name in cons_dict.keys():
                    chrom, period, left_boundary, right_boundary= fasta_name.split('/')[-1].split('_')
                    cons_num, cons_len, cons_count, left_boundary, right_boundary, mean_info, seq_div, complexity=cons_name.split('_')
                    if int(cons_count)==1: continue
                    seq, info= cons_dict[cons_name]
##                    true_intervals[major_period].append((left_boundary, right_boundary))
                    seq_name='{0}_{1}_{2}_{3}_{4}'.format(chrom,cons_len, cons_count,  left_boundary, right_boundary )
                    fasta_handle.write('>{0}\n'.format(seq_name))
                    fasta_handle.write('{0}\n'.format(seq))

                    row=[key, cons_len, cons_count, mean_info, seq_div, complexity,  left_boundary, right_boundary, info, seq]
                    annotation_table.writerow(row)
        PlotTandems(true_intervals)
        image_file='{0}/{1}.png'.format(image_dir, key)
        pyplot.savefig(image_file)
##        pyplot.yscale('symlog')
##        image_file='{0}/{1}_log.png'.format(image_dir, key)
##        pyplot.savefig(image_file)
        pyplot.close()
    annotation_handle.close()
    fasta_handle.close()
##        print key
def ExtractTandems(infile, outfile):
    sequences=GetSeq(infile, upper=True)#, rename=True)
    for key in sorted( sequences.keys()):
        out_dir='/'.join(outfile.split('/')[:-1])
        out_root='.'.join( outfile.split('/')[-1].split('.')[:-1])
        outimage='{0}/{1}_{2}.png'.format(out_dir, out_root, CleanName( key))

        p,w= FindPeriodicity(sequences[key], outfile, key)
        print key
##        edge=int( key.split('_')[-1].split('-')[0])
        edge=0
        if len(p)>0 or (numpy.array(p)>0).sum>0:
            pyplot.scatter(numpy.array(w)+edge,numpy.log10( p), c='r', alpha=.5)
            pyplot.ylabel('Repeat Length (Log10)')

##            pyplot.yscale('log')
            pyplot.xlabel('Chromosome Position')
            try:
                pyplot.savefig(outimage)
            except:
                a=1
            pyplot.close()

def AverageCorrAtPeriod(autocorr):
    mean_corr=[]
    for i in range(len(autocorr)/2):
        mean_corr.append(numpy.mean(autocorr[numpy.arange(len(autocorr))%i==0][1:]))
    return mean_corr


def FindBestPhase(seq, repeat_size, identity_threshold):
    char_array=numpy.fromstring(seq.upper(), '|S1')
    char_array=char_array[char_array!='N']
    char_mat=ArrayToMatrix(char_array, repeat_size)
    #Convert this into complex numbers for ease of calculation
    seq_array=numpy.array([0.]*len(char_array), 'complex')
    #A/T is the real plane
    seq_array[char_array=='A']+=1
    seq_array[char_array=='T']-=1

    #G/C is the imaginary plane
    seq_array[char_array=='G']+=1j
    seq_array[char_array=='C']-=1j

    #divide the sequence into non-overlapping kmers of the repeat size
    #storing these as an ndarray array where axis-0 represent kmers
    #and axis-1 represents
    shifted_mat=ArrayToMatrix(seq_array, repeat_size)

    #Compute percent identity between neighboring k-mers
    kmers,length=shifted_mat.shape
##    print char_mat.shape
##    print shifted_mat.shape
    stored_kmers=[]
    mode_list=[]
    nonrpt=[]
    len_list=[]

    for k in range( kmers-2):
        perc_ident=numpy.mean( shifted_mat[k,:]==shifted_mat[k+1,:])
        if perc_ident>identity_threshold:
            stored_kmers.append(char_mat[k,:])
        else:
            if len(stored_kmers)>0:
                stored_kmers.append(char_mat[k,:])
                mode=''.join(scipy.stats.mode(stored_kmers).mode[0])
                len_list.append(len(stored_kmers))
                mode_list.append(mode)
                stored_kmers=[]
            else:
                nonrpt.append(char_mat[k,:])
##    perc_ident=numpy.mean( shifted_mat[-2,:]==shifted_mat[-1,:])
##    if perc_ident>identity_threshold:
##        stored_kmers.append(char_mat[k,:])
    if len(stored_kmers)>0:
        stored_kmers.append(char_mat[-1,:])
        mode=''.join(scipy.stats.mode(stored_kmers).mode[0])
        len_list.append(len(stored_kmers))
        mode_list.append(mode)
        stored_kmers=[]
    else:
        try:
            nonrpt.append(char_mat[-1,:])
        except:
            a=1

    return mode_list, ''.join(numpy.array(nonrpt).flatten()), len_list



def FindPeriodicity(seq,outfile,seq_key='', window_size=50000, step_size=50000):

##    print len(seq)
    period_list=[]
    window_list=[]
    corr_list=[]
    mode_list={}
    len_dict={}
    windows_containing_period={}
    #Test over range of window sizes: 1000, 10000, 100000, 1000000
    for window_size in [10000,50000,100000,500000]:
        windows_containing_period[window_size]={}
        step_size=window_size/2
        max_repeat=window_size/2
        window_edges=numpy.arange(0,len(seq)-window_size+step_size, step_size)

        for i in range(len(window_edges)-1) :
            seq_slice=seq[window_edges[i]:window_edges[i]+window_size]

            for j in range(5):
                if len(seq_slice)<100: break
                try:
                    autocorr=Autocorrel(seq_slice,max_repeat)

                except:
                    break

                true_period=numpy.argmax(autocorr[1:])+1

                #Only interested in base periodicities>3
                if true_period>=3:
                    #Only going to look for periodicities>20
                    period=numpy.argmax(autocorr[3:])+3
                    #Note that this periodicity was identified




                    #Remove repeats of this length from the sequence
                    len_seq=len(seq_slice)
                    seq_slice=RemoveRepeatsOfSize(seq_slice, period, .8)
                    if len_seq==len(seq_slice): break
                    if windows_containing_period[window_size].has_key(period)==False:
                        windows_containing_period[window_size] [period]=[]
                    windows_containing_period[window_size] [period].append(i)
                    period_list.append(period)

                            #Remember how many of instances of this repeat were identified
                else: break
        #Consolidate intervals within window length
##
        for period in windows_containing_period[window_size]:
            windows_containing_period[window_size] [period]=[(window_edges[l],window_edges[r]+window_size) \
            for l,r in  LookForSplits(sorted( list( set( windows_containing_period[window_size] [period]))))]

    #Group intervals across windows by period:
    period_dict={}
    for window_size in windows_containing_period.keys():
        for period in windows_containing_period[window_size].keys():
            if period_dict.has_key(period)==False:
                period_dict[period]=[]
            for interval in windows_containing_period[window_size][period]:
                period_dict[period].append(interval)
##    return period_dict
    #Group overlapping windows:
    for period in period_dict.keys():
        period_dict[period]=GroupOverlappingIntervals(period_dict[period])

    return period_dict


def GroupOverlappingIntervals(intervals):
    sorted_intervals=sorted(intervals, key=lambda x:x[1])
    left_list, right_list=zip(* sorted_intervals)

    left_list=numpy.array(left_list)
    right_list=numpy.array(right_list)

    cluster_list=[]
    count=0
    while len(left_list)>1:
        count+=1
        ind=( left_list<=right_list[0])

        if sum(ind)>1:
            new_interval=(min(left_list[ind]), max(right_list[ind]))
            left_list=numpy.hstack((new_interval[0],  left_list[~ind]))
            right_list=numpy.hstack((new_interval[1], right_list[~ind]))
            sort_ind=numpy.argsort(right_list)
            left_list=left_list[sort_ind]
            right_list=right_list[sort_ind]
##            print new_interval

        else:
            cluster_list.append((left_list[0], right_list[0]))
            left_list=left_list[1:]
            right_list= right_list[1:]


    cluster_list.append((left_list[0], right_list[0]))
    return cluster_list


def PlotTandems(per_dict):
##    for scale in per_dict.keys():
    for rpt_len in per_dict.keys():
        if rpt_len<=0: continue
        for l,r in per_dict[rpt_len]:
            pyplot.plot((l,r), (rpt_len, rpt_len), c='black', alpha=1)
    pyplot.ylabel('Repeat Unit Length')
    pyplot.xlabel('Position')
##    pyplot.show()
##def PlotTandems(per_dict):
##    for scale in per_dict.keys():
##        for rpt_len in per_dict[scale].keys():
##            for l in per_dict[scale][rpt_len]:
##                pyplot.scatter(l, rpt_len, c='red', alpha=.8)
##    pyplot.show()
def RemoveRepeatsOfSize(seq, rpt_len, threshold):
    seq_array=numpy.fromstring(seq, '|S1')
    identity_signal=ShiftedIdentity(seq, rpt_len)



    high_identity_intervals=IdentifyHighIdentityRegions(identity_signal,threshold)

    for interval in high_identity_intervals:
        l,r=interval
        r+=rpt_len
        seq_array[l:r]='N'
    unmasked_ind=seq_array!='N'

    return ''.join(seq_array[unmasked_ind])


def MaskRepeatsOfSize(seq, rpt_len, threshold):
    seq_array=numpy.fromstring(seq, '|S1')
    identity_signal=ShiftedIdentity(seq, rpt_len)

    high_identity_intervals=IdentifyHighIdentityRegions(identity_signal,threshold)

    for interval in high_identity_intervals:
        l,r=interval
        r+=rpt_len
        seq_array[l:r]='N'


    return ''.join(seq_array)

    consensus_seq=[]
    info_list=[]
    out_dict={}
    for period in sorted( mode_list.keys()):
        out_dict[period]={}
        for rpt_length in mode_list[period]:
            out_dict[period][rpt_length]=[[],[]]
            if len(mode_list[period][rpt_length])>1:
                len_list=[max( len_dict[k].values()) for k in list(mode_list[period][rpt_length]) ]
##                print period, rpt_length
##                print len_list
                consensus, counts=ClusterSet(list(mode_list[period][rpt_length]),.8, len_list)
##                print counts
            else:
                size=sum([  max( len_dict[k].values()) for k in list(mode_list[period][rpt_length])])
                consensus, counts=list(mode_list[period][rpt_length]), [size]

##            consensus_seq+=(consensus)
            for j in range(len(counts)):
                out_dict[period][rpt_length][0].append(consensus[j])
                out_dict[period][rpt_length][1].append((period, rpt_length, counts[j]))
##                info_list+=[(period, rpt_length, counts[j])]
    outhandle=open(outfile, 'a')
    count=0
##    for i in range(len(consensus_seq)):
    for p in sorted( out_dict.keys()):
        for r in sorted( out_dict[p].keys()):
            info_list=out_dict[p][r][1]
            consensus_seq=out_dict[p][r][0]
            for i in range(len(out_dict[p][r][1])):
                count+=1
##                print info_list
                outhandle.write('>{0}_{1}_monomer={2}_period={3}_count={4}\n'.format(seq_key, count, info_list[i][0], info_list[i][1], info_list[i][2]))
                outhandle.write('{0}\n'.format(consensus_seq[i]))
    outhandle.close()
    return period_list, window_list

def ClusterFASTA(infile, outfile):
    sequences=GetSeq(infile, clean=False)
    length_dict={}
    count_dict={}
    chr_dict={}
    for key in sequences:
        length=len(sequences[key])
        if length_dict.has_key(length)==False:
            length_dict[length]=[]
            count_dict[length]=[]
            chr_dict[length]=[]
        length_dict[length].append(sequences[key])

        counts=float(SeqManipulations.ParseSeqName(key)['count'])
        chrom=key.split('_')[0]

        count_dict[length].append(counts)
        chr_dict[length].append(chrom)

    consensuses=[]
    counts=[]
    chr_list=[]
    for key in length_dict:

        consensus, count, chrom=ClusterSet(length_dict[key],.8, count_dict[key], chr_dict[key])
        consensuses+=consensus
        counts+=count
        chr_list+=chrom

    outhandle=open(outfile, 'w')
    count=0
    for j in range(len( consensuses)):
        c=consensuses[j]
        instances=counts[j]
        if instances*len(c)<300: continue
        count+=1
        name='>{0}_length={1}_counts={2}_chrom={3}\n'.format(count, len(c), instances, ';'.join(chr_list[j]))
        outhandle.write(name)
        outhandle.write('{0}\n'.format(c))

    outhandle.close()
##
##def ClusterSet(sequences, threshold=.8, len_list=[],chr_names=[]):
##    clusters=[]
##    return_chr=chr_names!=[]
##    terminate=False
##    count=0
##    count_list=[]
##    chr_list=[]
##    while len(sequences)>0 or count<100 :
##        count+=1
##        try:
##            clusters.append([numpy.fromstring(sequences[0], '|S1')])
##            count_list.append([len_list[0]])
##            chr_list.append(set([chr_names[0]]))
##        except:
##            break
##        cluster_ind=[0]
##        for s in range(1,len(sequences)):
##            s1,s2,perc_id=QuickMatch( sequences[s],sequences[0])
##            sc1,sc2,perc_id_comp=QuickMatch(str(Seq.Seq( sequences[s]).reverse_complement()),sequences[0])
##            if perc_id>=threshold or perc_id_comp>threshold:
##                if perc_id>perc_id_comp:
##                    clusters[-1].append(s1)
##                    cluster_ind.append(s)
##                    if chr_names!=[]:chr_list[-1].add(chr_names[s])
##                    if len_list!=[]:
##                        count_list[-1].append(len_list[s])
##                    else:
##                        count_list[-1].append(1)
##                else:
##                    clusters[-1].append(sc1)
##                    cluster_ind.append(s)
##                    if chr_names!=[]:chr_list[-1].add(chr_names[s])
##                    if len_list!=[]:
##                        count_list[-1].append(len_list[s])
##                    else:
##                        count_list[-1].append(1)
##        #Remove clusters from sequence
##
##        for s in sorted( cluster_ind, reverse=True):
##            del sequences[s]
##            del len_list[s]
##            del chr_names[s]
####        print len(sequences)
####        if count>3: return clusters
##    consensus_list=[]
##    counter=[]
##    for c in range(len( clusters)):
##        cluster_array=numpy.array(clusters[ c])
####        print cluster_array.shape
####        print scipy.stats.mode( numpy.array(c),0).mode
##        consensus_list.append(''.join (scipy.stats.mode( numpy.array(clusters[ c]),0).mode[0]))
##        counter.append(sum(count_list[c]))
##
##    #consensus_list =  list of consensus sequences
##    #counter = list indicating how many mononmers went into each consensus
##    if return_chr==False:
##        return consensus_list, counter
##    else:
##        return consensus_list, counter, chr_list



def ClusterSet(sequences, threshold=.8, len_list=[],chr_names=[]):
    clusters=[]
##    print len(sequences)
    return_chr=chr_names!=[]
    terminate=False
    count=0
    count_list=[]
    chr_list=[]
    while len(sequences)>0  :
        count+=1
##        try:
        clusters.append([numpy.fromstring(sequences[0], '|S1')])
        if len_list!=[]:
            count_list.append(len_list[0])
        else: count_list.append(1)
        if chr_names!=[]:
            chr_list.append(set([chr_names[0]]))
##        except:
##            break
        cluster_ind=[0]
        for s in range(1,len(sequences)):
            #Find the best alignment betwen the sequecnes using a heuristic
            s1,s2,perc_id=QuickMatch( sequences[s],sequences[0])
            sc1,sc2,perc_id_comp=QuickMatch(str(Seq.Seq( sequences[s]).reverse_complement()),sequences[0])

            #Check whether the percent identity of the match exceeds a threshold
            if perc_id>=threshold or perc_id_comp>threshold:
                if perc_id>perc_id_comp:
                    clusters[-1].append(s1)
                    cluster_ind.append(s)
                    if chr_names!=[]:chr_list[-1].add(chr_names[s])
                    if len_list!=[]:
                        count_list[-1]+=len_list[s]
                    else:
                        count_list[-1]+=1
                else:
                    clusters[-1].append(sc1)
                    cluster_ind.append(s)
                    if chr_names!=[]:chr_list[-1].add(chr_names[s])
                    if len_list!=[]:
                        count_list[-1]+= len_list[s]
                    else:
                        count_list[-1]+=1
        #Remove clusters from sequence
##        print count_list[-1],
        for s in sorted( cluster_ind, reverse=True):
            del sequences[s]
            if len_list!=[]:
                del len_list[s]
            if chr_names!=[]:
                del chr_names[s]
##        print len(sequences)
##        if count>3: return clusters
    consensus_list=[]
    counter=[]
    for c in range(len( clusters)):
        cluster_array=numpy.array(clusters[ c])
##        print cluster_array.shape
##        print scipy.stats.mode( numpy.array(c),0).mode
        consensus_list.append(''.join (scipy.stats.mode( numpy.array(clusters[ c]),0).mode[0]))
        counter.append(count_list[c])

    #consensus_list =  list of consensus sequences
    #counter = list indicating how many mononmers went into each consensus
    if return_chr==False:
        return consensus_list, counter
    else:
        return consensus_list, counter, chr_list

def QuickMatch(seq1, seq2, seeds=100):
    assert len(seq1)==len(seq2)
    if len(seq1)>20:
        kmer_set_1=set(GetKMERS(seq1))
        kmer_set_2=set(GetKMERS(seq2))
        seeds=numpy.random.choice(list(kmer_set_1), size=seeds, replace=True)
        seed_intersection=list( set(seeds)&kmer_set_2)
    ##    print len(set(seeds))
    ##    print len(seed_intersection)
        disp=[]
        if len(seed_intersection)==0:
            return [], [], 0
        for s in seed_intersection:
            disp.append(seq1.find(s)-seq2.find(s))
        mode=scipy.stats.mode(disp).mode[0]
    ##    print mode

    else:#If it is short, use brute force: check every rotation
        mode_list=[]
        perc_id=[]
        for i in range(len(seq1)):
            rotated_seq=seq1[i:]+seq1[:i]
            seq1_array, seq2_array=numpy.fromstring(rotated_seq, '|S1'),numpy.fromstring( seq2, '|S1')
            mode_list.append(i)
            perc_id.append( numpy.mean( seq1_array== seq2_array))
        mode=mode_list[ numpy.argmax(perc_id)]

    seq1_array, seq2_array=numpy.fromstring( seq1[mode:]+seq1[:mode], '|S1'),numpy.fromstring( seq2, '|S1')
    return seq1_array, seq2_array , numpy.mean( seq1_array== seq2_array)




def BoothsAlgorithm(s):
    #Booth's algoritm for finding the lexicographical minimum string rotation
    #From Wiki
    S = s*2      # Concatenate string to it self to avoid modular arithmetic
    f = [-1] * len(S)     # Failure function
    k = 0       # Least rotation of string found so far
    for j in xrange(1,len(S)):
        sj = S[j]
        i = f[j-k-1]
        while i != -1 and sj != S[k+i+1]:
            if sj < S[k+i+1]:
                k = j-i-1
            i = f[i]
        if sj != S[k+i+1]: # if sj != S[k+i+1], then i == -1
            if sj < S[k]: # k+i+1 = k
                k = j
            f[j-k] = -1
        else:
            f[j-k] = i+1
    length=len(s)

    return s[k-length:]+s[:k]



    return period_list, window_list, corr_list, mode_list

def PowerLaw(x,a,k,b):
    return a* numpy.exp(-k*x)+b

def SubstractSignal(signal, periodicity):
    ind=periodicity*numpy.arange(0, len(signal)/periodicity)[1:]
    signal_at_periods=signal[ind]
    pyplot.scatter(( ind), ( signal_at_periods))
    pyplot.show()
    opt=scipy.optimize.curve_fit(PowerLaw, ind, signal_at_periods,p0=(1,0,.5))
    print opt
    curve=PowerLaw(numpy.arange(len(signal)),*opt[0])
    pyplot.plot(signal, c='b')
    pyplot.plot(numpy.arange(len(signal)), curve)
    pyplot.show()

##def Autocorrel(seq, max_size, p=False):
####    seq_array=numpy.fromstring(seq.upper(), '|S1')
##
##    char_array=numpy.fromstring(seq.upper(), '|S1')
##    char_array[char_array!='N']
####    return LaggedIdentity(seq_array, max_size)
####    AT=numpy.array([0]*len(seq_array))
####    GC=numpy.array([0]*len(seq_array))
####    AT[(seq_array=='A')]=1
####    AT[(seq_array=='T')]=-1
####    GC[(seq_array=='G')]=1
####    GC[(seq_array=='C')]=-1
##
##    seq_array=numpy.array([0.]*len(char_array), 'complex')
##    #A/T is the real plane
##    seq_array[char_array=='A']+=1
##    seq_array[char_array=='T']-=1
##
##    #G/C is the imaginary plane
##    seq_array[char_array=='G']+=1j
##    seq_array[char_array=='T']-=1j
####    return seq_array
####    acf_AT, q_AT, p_AT=acf(AT,nlags=max_size,qstat=True, fft=True)
####    acf_GC, q_GC, p_GC=acf(AT,nlags=max_size,qstat=True, fft=True)
####    pyplot.plot(q_AT)
####    pyplot.show()
####    autocorr=numpy.sign(acf_AT* acf_GC)* ( acf_AT* acf_GC)**.5
####    autocorr=(acf_AT+acf_GC)/2.
##    autocorr=acf(seq_array,unbiased=False, nlags=max_size, fft=True)
##
####    autocorr_2=acf(autocorr,nlags=max_size, fft=True)
##
####    autocorr[(p_AT>.01)*(p_GC>.01)]=0.
##    return autocorr

def SeqToComplex(seq):
    char_array=numpy.fromstring(seq.upper(), '|S1')
    char_array[char_array!='N']
##    return LaggedIdentity(seq_array, max_size)
##    AT=numpy.array([0]*len(seq_array))
##    GC=numpy.array([0]*len(seq_array))
##    AT[(seq_array=='A')]=1
##    AT[(seq_array=='T')]=-1
##    GC[(seq_array=='G')]=1
##    GC[(seq_array=='C')]=-1

    seq_array=numpy.array([0.]*len(char_array), 'complex')
    #A/T is the real plane
    seq_array[char_array=='A']+=1
    seq_array[char_array=='T']-=1

    #G/C is the imaginary plane
    seq_array[char_array=='G']+=1j
    seq_array[char_array=='C']-=-1j
    return seq_array

def CCF(sig1, sig2):
    fft_1=scipy.fft(sig1)
    fft_2=scipy.fft(sig2)
    CCF=fft_1*numpy.conj(fft_2)
    return CCF


def ShiftedIdentity(target, repeat_len):
    """Compares a sequence with a shifted version of itself. Good for visualizing"""

    seq_array=numpy.fromstring(target, '|S1')

    matches=(seq_array[:-repeat_len] == seq_array[repeat_len:])
    #Ensure Ns are always treated as mismatches
    N_pos=( (seq_array[:-repeat_len]=='N')+ (seq_array[repeat_len:]=='N') )>0
    matches[N_pos]=0

    smoothing_kernel=[1./repeat_len]*repeat_len
    perc_id=scipy.signal.fftconvolve(smoothing_kernel, matches)
    return perc_id

def LookForSplits(indices):
    indices=numpy.array(list(indices))

    diff_ind=numpy.diff(indices)

    splits=[0]+ list(numpy.where(diff_ind>1)[0])+[len(indices)-1]
    intervals=[(indices[splits[0]], indices[splits[1]]) ]
##    print indices


    for i in range(1,len(splits)-1):
        intervals.append((indices[splits[i]+1], indices[splits[i+1]]))

    return intervals
def IdentifyHighIdentityRegions(identity, threshold):
    indices=numpy.where(identity>=threshold)[0]
    indices=numpy.array([-2]+list(indices)+[len(identity)+1])

    diff_ind=numpy.diff(indices)

    splits=numpy.where(diff_ind>1)[0]
    intervals=[]
    for i in range(len(splits)-1):
        intervals.append(( indices[splits[i]+1], indices[ splits[i+1]]))
    return intervals

def ExtractRepeatsOfSize(seq, rpt_len, threshold=.8):

    identity_signal=ShiftedIdentity(seq, rpt_len)
##    return identity_signal
    high_identity_intervals=IdentifyHighIdentityRegions(identity_signal, threshold)
    #The intervals identified go from the end of the first repeat in the array
    #to the beginning of the second repeat.
##    return high_identity_intervals
    boundaries=[]
    repeats=[]

    for interval in high_identity_intervals:

        l,r=interval
        if r-l<1: continue
        #The boundaries of the array need to be inferred, because regions matching
        #the repeat will be surrounded by slopes of decaying identity:
        #    _______
        #   /       \
        #  /         \
        #
        #The intervals identified contain some of the sloping regions. To remove
        #at least some of these regions, we first summarize the distribution
        #of percent identity within the interval.
        mean_PI=numpy.mean(identity_signal[l:r])
        std_PI=numpy.std(identity_signal[l:r])

        #We define a cutoff as identity more than two standards below the mean
        local_cutoff=mean_PI-2*std_PI

        #We refine the interval to the first and last positions where identity
        #is below the cutoff
        begin=l+ numpy.where(identity_signal[l:r]>=local_cutoff)[0][0]
        end=l+ numpy.where(identity_signal[l:r]>=local_cutoff)[0][-1]
        boundaries.append((begin, end))
        num_rpts=int( (end-begin)/rpt_len)
        #The we determine the largest tandem of complete repeat units that could
        #fit in this interval
        expected_array_size=rpt_len*num_rpts

        #The set of positions within this array could start is:
        #(begin, end - expected_array_size)
        #We assume each is equally likely and take the expected value as
        #the array's start.
        uncertainty=int(((end-begin)-expected_array_size)/2)
        repeat_begin= begin+ uncertainty

        #We extract the second repeat in the array. The above determination of
        #the starting point is likely to make some error, but we expect the true
        #juncion to be near the edge of the repeat.
        repeats.append(( seq[repeat_begin:repeat_begin+rpt_len],(repeat_begin,repeat_begin+rpt_len)  ))


    return repeats

def GetRepeatsInInterval(seq, period, threshold, masked=None):
    if masked!=None:
        candidate_rpts=ExtractRepeatsOfSize(masked, period, threshold)
    else:
        candidate_rpts=ExtractRepeatsOfSize(seq, period, threshold)
    candidate_rpts=list(set(candidate_rpts))
    rpt_dict={}
    masked_seq=copy.copy( seq)
    checked_rpts=set()
    identified_rpts=[0]*3
    last_len=len(candidate_rpts)+1
    while len( candidate_rpts)>0:
##    for rpt_info in candidate_rpts:
##    while len()

        rpt,interval=candidate_rpts[0]
        l,r=interval
        print "N prop=", float( masked_seq[l:r].count('N'))/(r-l)
##        if float( masked_seq[l:r].count('N'))/(r-l)>.1: continue

        if len( checked_rpts&set(rpt))>0:
            del candidate_rpts[0]
            continue
        checked_rpts.add(rpt)

        identified_rpts, masked, rpt_len=ExtractRepeatFromArray(seq, masked_seq, rpt, threshold=threshold)
    ##            del masked_seq
        masked_seq=''.join( masked)
        for rpt_len in identified_rpts.keys():
            if rpt_dict.has_key(rpt_len)==False:
                rpt_dict[rpt_len]=[]
            print rpt_len, len(identified_rpts[rpt_len])


            rpt_dict[rpt_len]+=identified_rpts[rpt_len]

        print len(rpt), "Masked sequence:", (numpy.fromstring(masked_seq, '|S1')=='N').sum()
        identified_rpts=[1]*3
        candidate_rpts=ExtractRepeatsOfSize(masked_seq, period, threshold)
        if len(candidate_rpts)>=last_len:
            break
        last_len=len(candidate_rpts)
##        candidate_rpts=list(set(candidate_rpts)-checked_rpts)
##        repeat_list+=identified_rpts
##        if len(identified_rpts)==0: break
    return rpt_dict, masked_seq


#m=(y'-y)/(x'-x)

def PickPeriods(autocorr):
    best_period=numpy.argmax(autocorr[1:])+1
    length_auto=len(autocorr)
    slope=(1*autocorr[ best_period])/(length_auto-best_period)
    print slope
    beta=slope*length_auto
    return test-( numpy.arange(length_auto)*slope+beta)

def FindAdditionalPeriods(identity_signal ,distance_list,exclude, threshold, cutoff=3):
    counter=collections.Counter(distance_list)

    period_list=[]
    count_list=[]
    period_set=set()
    for key in counter.keys():
        if list(period_set).count(key)>0: continue
##        if counter[key]<cutoff: continue
        lengths, identities=AnalyzeHitsAtPeriod(identity_signal, key, threshold=threshold)
        if len(lengths)<cutoff: continue
        period_set|=set( lengths)
        average_ident=numpy.median(identities)
        print "\tPeriod: ", key, " Stats:", numpy.sum(lengths), average_ident, numpy.sum(lengths) * average_ident

        count_list.append(numpy.sum(lengths) * average_ident)
        period_list.append(key)
    sort_ind=numpy.argsort(count_list)[::-1]

    pyplot.plot( identity_signal)
    pyplot.show()
    return list (numpy.array(period_list)[sort_ind])


def AnalyzeHitsAtPeriod(identity_signal, expected_periodicity, threshold ):

    hit_length=[]
    hit_identity=[]
    #This counts the number of consecutive repeats
    consecutive_repeats=0
    hits=numpy.where(identity_signal>=threshold)[0]
    if len(hits)==0:
        return
    distance_between_hits=numpy.diff(hits)
##        distance_between_hits=numpy.hstack((distances))
    dist_tracker=0
    for index, distance in enumerate( distance_between_hits):
        #If the distance to the next repeat is within 105% of the expected repeat length
        #this is part of a tandem array. Extract the repeat, mask it, and increment
        #the consecutive repeat counter
        tandem=False

        if distance<=expected_periodicity*1.15 and distance>expected_periodicity*.85:
            hit_length.append(distance)
            hit_identity.append(identity_signal[hits[index]] )
            consecutive_repeats+=1

        elif consecutive_repeats>0:
            hit_length.append(distance)
            hit_identity.append(identity_signal[hits[index]] )
            #stuff
            consecutive_repeats=0


    if consecutive_repeats>0:
        hit_length.append(expected_periodicity)
        hit_identity.append(identity_signal[hits[-1]] )
    return hit_length, hit_identity


def ExtractRepeatFromArray(sequence, masked_sequence, query, threshold=.8, min_length=30):

    #Note to self: Want to retain information about their locations

    #Only search for the query repeat
    identity_signal=LaggedUngappedAlignment(query, masked_sequence)
##    pyplot.plot(identity_signal)
##    pyplot.show()
    seq_array=numpy.fromstring(masked_sequence, '|S1')

    #It obvious from examination of identity_signals that heuristics could be
    #devised to identify indels. Such a heuristic should be incorporated here.

    hits=numpy.where(identity_signal>=threshold)[0]

    #If no hits, terminate
    if len(hits)==0:
        return [], masked_sequence,0

    #Compute the distance between adjacent hits--the actual periodicity of the repeat
    distances=numpy.diff(hits)

    #Determine the most frequent periodicity

    #If the most frequent periodicity is larger than the length of the query repeat
    #the query might be part of a higher-order repeat. However, if there aren't a lot
    #of hits, this might just reflect dispersed repeats or large insertions, and
    #we don't want to pick those up. So, impose a count threshold.
    #So, if the most frequent periodicity is longer than expected and occurs at least
    #3 times, update the expected periodicty to reflect that, otherwise determine
    #expected periodicity from the query repeat
##    periods_list=[len(query)]
##    periods_list=FindAdditionalPeriods(identity_signal, distances, len(query), threshold)
    repeat_list={}
##    if periods_list==[]: periods_list=[len(query)]
    possible_periods=ConstructModeTree(identity_signal,.8)
    if len( possible_periods)>0 and len(query)<300:
        period_list=OrganizePeriods(possible_periods)
        period_list.append((len(query), threshold,3))
    else: period_list=[(len(query), threshold,3)]
    for expected_periodicity,threshold,rpt_counts in period_list:
        if expected_periodicity>50*len(query): continue
        if rpt_counts<3: continue
        expected_periodicity=int(expected_periodicity)

        #This counts the number of consecutive repeats
        consecutive_repeats=0
        hits=numpy.where(identity_signal>=threshold)[0]
        if len(hits)<=1:
            continue
        distance_between_hits=numpy.diff(hits)
##        distance_between_hits=numpy.hstack((distances))
        for index, distance in enumerate( distance_between_hits):
            #If the distance to the next repeat is within 105% of the expected repeat length
            #this is part of a tandem array. Extract the repeat, mask it, and increment
            #the consecutive repeat counter
            tandem=False
            if distance<=min_length:
                #Really want to ignore very simple repeats. Skip AND mask.
                left, right=hits[index],hits[index]+distance
                seq_array[left:right]='N'
##                identity_signal[hits[index]]=0
                identity_signal[left:right]=0
                continue
            if distance<=expected_periodicity*1.15 and distance>expected_periodicity*.85:
                left, right=hits[index],hits[index]+distance
                tandem=True
                consecutive_repeats+=1

            elif consecutive_repeats>0:
                left, right= hits[index],hits[index]+expected_periodicity
                tandem=True
                consecutive_repeats=0
##
##            else:
##                left, right= hits[index],hits[index]+expected_periodicity
##                tandem=False
            if tandem==True:

                repeat_seq=sequence[left: right]
                if repeat_list.has_key(expected_periodicity)==False:
                    repeat_list[expected_periodicity]=[]

                repeat_list[expected_periodicity].append(ReferenceRepeat(repeat_seq, left, right,tandem))
                seq_array[left:right]='N'
                identity_signal[left:right]=0




##                    identity_signal[hits[index]]=0
        if consecutive_repeats>0:
                left, right= hits[-1],hits[-1]+expected_periodicity
                tandem = True
                consecutive_repeats = 0
                repeat_seq = sequence[left: right]
                if repeat_list.has_key(expected_periodicity)==False:
                    repeat_list[expected_periodicity]=[]

                repeat_list[expected_periodicity].append(ReferenceRepeat(repeat_seq, left, right,tandem))
                seq_array[left:right]='N'
                identity_signal[left:right]=0
##                identity_signal[hits[index]]=0
##    pyplot.plot(identity_signal)
##    pyplot.show()
    if len(repeat_list.keys())>0:
        period_keys=HierarchicalCluster(repeat_list.keys(),.1)

    else: return {}, seq_array, expected_periodicity


    repeat_dict={}
    for key_list in period_keys:
        mean_key=int(numpy.mean(key_list))
        repeat_dict[mean_key]=[]
        for key in key_list:
            repeat_dict[mean_key]+=repeat_list[key]
    masked_sequence=''.join(seq_array)
    return repeat_dict, seq_array, expected_periodicity

def LaggedUngappedAlignment(query, target):

    """"Computes the percent identity for every possible ungapped alignment
    between the query and target. It does this by encoding each nucleotide with
    a sequence of four numbers:

        A: 1000
        T: 0100
        G: 0010
        C: 0001

    And then computing the cross-correlation function of the two signals.
    Every fourth position in the CCF indicates the number of matching nucleotides
    in a window equal to the length of the query, because the cross-correlation
    of F and G is:

        Conv( F, G_) (t)=sum{for all T} F(T)*G_(T-t)
        where G_ is the time-reversed version of G.

    The product F(T)*G_(T-t) will only equal one if the nucleotides match.

    To get percent identity, the number of matches is divided by query length.
        """

    #To avoid memory errors if the target is longer than 100000, we split it and
    #stitch the results together
    if len(target)<200000:

        query_length=len(query)
        sig_1=SeqToExanded(query)
        sig_2=SeqToExanded(target)
        matches=scipy.signal.fftconvolve(sig_1, sig_2[::-1])[3::4]
        perc_id=matches/query_length
        return perc_id[query_length-1:-query_length+1][::-1]

    else:
        split_indices=numpy.arange(0, len(target), 100000)

        if max(split_indices)!=len(target):
            split_indices=numpy.hstack((split_indices, len(target) ))
##        print split_indices

        for i,v in enumerate(split_indices[:-1]):
##            print i, v

            if i==0:
                sliced_seq=target[0:split_indices[i+1]]
                percent_id=LaggedUngappedAlignment(query, sliced_seq)

            else:
##                print v-len(query),split_indices[i+1]
                sliced_seq=target[v-len(query):split_indices[i+1]]
                percent_id=numpy.hstack((percent_id,LaggedUngappedAlignment(query, sliced_seq)[1:] ))
        return percent_id

def SetMUSCLEPath(path):
    global MUSCLE_PATH
    MUSCLE_PATH= path

def RunMuscle(infile, outfile, maxiters, logfile):
    command=[MUSCLE_PATH, '-in', str( infile), '-out',str( outfile), '-maxiters', str(maxiters), '-diags']
    errhandle=open(logfile, 'a')
    process=subprocess.Popen(command, stderr=errhandle)
    process.communicate()
    errhandle.close()

def RunProfileProfileAlignment(infile1, infile2, outfile, logfile):
    command=[MUSCLE_PATH, '-profile', '-in1', str( infile1), '-in2', str(infile2), '-out',str( outfile)]
    errhandle=open(logfile, 'a')
    process=subprocess.Popen(command, stderr=errhandle)
    process.communicate()
    errhandle.close()

def MultipleSequenceAlignment(infile, outfile, maxiters, logfile):

    seq=GetSeq(infile)
    if len (seq.keys())<500:
        RunMuscle(infile, outfile, maxiters, logfile)

    else:
        print "\tSplitting *.fasta..."
        #Split the output into 1000 sequence files
        temp_dir='{0}_temp'.format( '.'.join( infile.split('.')[:-1]))
        splitFile(infile, temp_dir)

        #Get MSA for each file
        file_list=os.listdir(temp_dir)
        out_file_list=[]
        for f in file_list:
            print "\tRunning MSA for {0}...".format(f)
            outname='{0}/{1}_out.fa'.format(temp_dir, '.'.join( f.split('.')[:-1]))
            RunMuscle('{0}/{1}'.format(temp_dir, f), outname, maxiters, logfile )
            out_file_list.append(outname)

        #Combine the MSAs with profile-profile alignments
        print "\tCombining MSAs..."
        RunProfileProfileAlignment(out_file_list[0],out_file_list[1], outfile, logfile)
        if len(out_file_list)>2:
            for i,v in enumerate(out_file_list[2:]):
                RunProfileProfileAlignment(v,outfile, outfile, logfile)
        print "\tRemoving temporary files..."
        shutil.rmtree(temp_dir)
        print '\tDone.'

def WriteFasta(sequences, outfile):
    outhandle=open(outfile, 'w')
    for i,s in  enumerate( sequences ):
        outhandle.write('>{0}_{1}_{2}\n'.format(s.chrom, s.left, s.right))
        outhandle.write('{0}\n'.format(s.sequence))
    outhandle.close()


def SummarizeRepeats(sequences, outfile, log_handle):
    """Takes a list of sequences and generates a consensus
    """

    fasta_output='{0}_repeats.fa'.format(outfile)

    #Limit this to tandem repeats
    tandem_repeats=[s  for s in sequences if s.tandem_flag]
    min_len=min( [abs(s.length) for s in tandem_repeats ])
    if len(tandem_repeats)==0: return {}

    #Output sequences to *fasta
    WriteFasta(sequences, fasta_output)

    #If there are more than two sequences run an MSA
    #Run the multiple alignment
    msa_output='{0}_msa.fa'.format(outfile)
    if len(tandem_repeats)>2:
        MultipleSequenceAlignment(fasta_output, msa_output, 1, log_handle)
    else:
        #If the sequences are short, use a small word size
        if min_len<100: word_size=7
        else: word_size=11
        AlignWithBLAST(fasta_output, msa_output,  log_handle,word_size)

    #Otherwise, run a blast alignment:

    #To do: Account for the possibility that not all repeats are well described
    #by one consensus
    #   Use the multiple to cluster sequences
    #       Not implemented
    #   Run multiple alignments on each cluster
    #       Not implemented

    #Build consensus
    consensus_dict=GetConsensusFromFasta(msa_output)

    return consensus_dict

def AlignWithBLAST(infile, outfile,log_file, word_size=11, blastdir=BLAST_PATH):
    sequences=GetSeq(infile)
    seq_keys=sequences.keys()
    temp_1='{0}_1.fa'.format('.'.join( infile.split('.')[:-1]))
    temp_handle=open(temp_1, 'w')
    temp_handle.write('>{0}\n'.format(seq_keys[0]))
    temp_handle.write('{0}\n'.format( sequences[seq_keys[0]]))
    temp_handle.close()

    temp_2='{0}_2.fa'.format('.'.join( infile.split('.')[:-1]))
    temp_handle=open(temp_2, 'w')
    temp_handle.write('>{0}\n'.format(seq_keys[1]))
    temp_handle.write('{0}\n'.format( sequences[seq_keys[1]]))
    temp_handle.close()

    #Align the temporary files
    temp_out='{0}_temp.xml'.format('.' .join(outfile.split('.')[:-1]))
    BlastSequences(temp_1, temp_2, temp_out, blastdir,word_size, log_file)

    #Delete the temporary files
    os.remove(temp_1)
    os.remove(temp_2)

    parse_handle=open(temp_out, 'r')
    blast_parser=NCBIXML.parse(parse_handle)
    hsps=blast_parser.next()
    outhandle=open(outfile,'w')
    if len(hsps.alignments)==0: #No alignments found
        outhandle.close()
        return

    query_name=hsps.query
    query_seq=hsps.alignments[0].hsps[0].query
    outhandle.write('>{0}\n'.format(query_name))
    outhandle.write('{0}\n'.format(query_seq))

    subj_name=hsps.alignments[0].hit_id
    subj_seq=hsps.alignments[0].hsps[0].sbjct
    outhandle.write('>{0}\n'.format(subj_name))
    outhandle.write('{0}\n'.format(subj_seq))

    outhandle.close()





def BlastSequences(query,subject,outfile,blastdir, word_size, logfile):
    errhandle=open(logfile, 'a')
    p=subprocess.Popen([blastdir, '-query', query,'-subject', subject,'-out', outfile, '-outfmt', '5', '-max_hsps', '1', '-word_size', str(word_size) ],stderr=errhandle)
    p.communicate()
    errhandle.close()

def ClusterMSA(sequences, threshold=.8,):
    clusters=[]
##    print len(sequences)
    copy_seq=copy.copy(sequences)
    terminate=False
    count=0
    count_list=[]
    chr_list=[]
    while len(sequences)>0  :
        count+=1

        clusters.append([sequences[0]])

##        except:
##            break
        cluster_ind=[0]
        for s in range(1,len(sequences)):
            #Find the best alignment betwen the sequecnes using a heuristic
            perc_id=ComputeMSADistance(sequences[ s].seq, clusters[-1][0].seq)
##            print perc_id
            #Check whether the percent identity of the match exceeds a threshold
            if perc_id>=threshold:
                clusters[-1].append(sequences[s])
                cluster_ind.append(s)

        #Remove clusters from sequence
##        print count_list[-1],
        for s in sorted( cluster_ind, reverse=True):
            del sequences[s]

    return clusters


def ComputeMSADistance(seq1, seq2):
    seqarray_1=numpy.fromstring(seq1, '|S1')
    seqarray_2=numpy.fromstring(seq2, '|S1')
    matches=seqarray_1==seqarray_2
    return matches.sum()/float(len(matches))

def GetConsensusFromFasta(infile):
    seq=GetSeq(infile)


    consensus_dict={}
    if len(seq.keys())==0:
        return consensus_dict
    class FastaSeq():
        def __init__(self, name,sequence):
            self.name=name
            self.left=int( name.split('_')[-2])
            self.right=int( name.split('_')[-1])
            self.seq=sequence
            self.length=len(sequence)
    seq_dict=[]
    for key in seq.keys():
        seq_dict.append( FastaSeq(key, seq[key]))


    clusters=ClusterMSA(seq_dict)


    for i, cluster in enumerate( clusters):
        left_edges=[c.left for c in cluster]
        right_edges=[c.right for c in cluster]
        lengths=[c.length for c in cluster]
        sequences=[s.seq for s in cluster]

        consensus, information, diversity=GetConsensusFromSequences(sequences)
        mean_information=numpy.nanmean(information)
        information_string=ConvertFloatToAscii(information, 0, 2)
        complexity=ComputeKmerCompositionEntropy(consensus)
        name='{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(i,len(consensus), len(cluster), min(left_edges), max(right_edges), mean_information, diversity, complexity )
        consensus_dict[name]=(consensus,information_string)
    return consensus_dict


def GetConsensusFromSequences(seq):
    cons_array=numpy.ndarray((5, len(seq[0])))
    cons_array.fill(0.)
    nt_dict={'-':0, 'A':1, 'C':2,'T':3, 'G':4}
    NT_array=numpy.array(['-', 'A', 'C', 'T', 'G'])
    for s in seq:

        seq_array=numpy.fromstring(s, '|S1')
        for i,v in enumerate(NT_array):
            ind=seq_array==v

            cons_array[i,ind]+=1.

    cons_seq=NT_array[numpy.argmax(cons_array,0)]
    information=SequenceInformation(cons_array)
    div=SequenceDiversity(cons_array)

    return ''.join(cons_seq[cons_seq!='-']), information, div

def ConvertFloatToAscii(x,float_lower, float_upper, ascii_lower=33., ascii_upper=126.):
    conversion_factor=(ascii_upper-ascii_lower)/(float_upper-float_lower)
    ascii_number=numpy.round( (x-float_lower)*conversion_factor+ascii_lower).astype(int)
    char_string=''.join([chr(i) for i in ascii_number ])
    return char_string

def ConvertAsciiToFloat(char_string,float_lower, float_upper, ascii_lower=33., ascii_upper=126.):
    F=numpy.array( [float( ord(c)) for c in char_string])
    conversion_factor=(float_upper-float_lower)/(ascii_upper-ascii_lower)
    ascii_number=(F-ascii_lower)*conversion_factor+float_lower

    return ascii_number

def SequenceInformation(cons_array):
    norm_const=cons_array[0:,:].sum(0)

    #Identify where indels
    gap_ind=numpy.argmax(cons_array,0)==0
    f_ai=(cons_array[1:,:]/norm_const)[:,~gap_ind]


    H=-1*( f_ai*numpy.log2(f_ai))
##    H[numpy.where(f_ai==0)]=0

    H_i=numpy.nansum(H,0)
    e_n_i=(1/numpy.log(2))*(3/(2*norm_const[~gap_ind]))
##    e_n=0
##    print e_n

    R_i=2-(H_i+e_n_i)
    R_i[R_i<0]=0.
    R_i[R_i>2]=2.

    return R_i

def SequenceDiversity(cons_array):
##    print cons_array
    norm_const=cons_array[1:,:].sum(0)
    #Identify where indels
    gap_ind=numpy.argmax(cons_array,0)==0
    N_a=(cons_array[1:,:][:,~gap_ind])
    N_k=norm_const[~gap_ind]
##    print (N_a*(N_a-1))/(N_k*(N_k-1))

    SD=numpy.nansum ((N_a*(N_a-1))/(N_k*(N_k-1)) )/float(len(N_k))
##    print SD
    return SD
def ConsensusFromClustal(alignment):
    nt_list=[]
    for col in range( alignment.get_alignment_length()):
        column=alignment[:,col]

        nt_array=numpy.fromstring(column, '|S1')


        mode=scipy.stats.mode(nt_array)[0][0]

        nt_list.append(mode)
    nt_list= numpy.array(nt_list)
    return ''.join(nt_list[nt_list!='-'])

def ACF_summary(seq, cutoff=.4):
    acf=LaggedUngappedAlignment( seq,seq*2)[1:-1]
    pyplot.plot(acf)
    pyplot.show()
    hits=(acf>=cutoff)
    print hits.sum()
    print acf[hits].sum()
    print numpy.mean(acf[hits])
    print numpy.mean(acf[hits])*hits.sum()
    print numpy.median( numpy.diff(numpy.where(acf>cutoff)))

def Autocorrel(seq, max_size, complex_num=True,unbiased=False,  p=False):
##    seq_array=numpy.fromstring(seq.upper(), '|S1')
    if complex_num==True:
        char_array=numpy.fromstring(seq.upper(), '|S1')
        char_array[char_array!='N']
    ##    return LaggedIdentity(seq_array, max_size)
    ##    AT=numpy.array([0]*len(seq_array))
    ##    GC=numpy.array([0]*len(seq_array))
    ##    AT[(seq_array=='A')]=1
    ##    AT[(seq_array=='T')]=-1
    ##    GC[(seq_array=='G')]=1
    ##    GC[(seq_array=='C')]=-1

        seq_array=numpy.array([0.]*len(char_array), 'complex')
        #A/T is the real plane
        seq_array[char_array=='A']+=1
        seq_array[char_array=='T']-=1

        #G/C is the imaginary plane
        seq_array[char_array=='G']+=1j
        seq_array[char_array=='C']-=-1j
    else:
        seq_array=SeqToExanded(seq)
##    return seq_array
##    acf_AT, q_AT, p_AT=acf(AT,nlags=max_size,qstat=True, fft=True)
##    acf_GC, q_GC, p_GC=acf(AT,nlags=max_size,qstat=True, fft=True)
##    pyplot.plot(q_AT)
##    pyplot.show()
##    autocorr=numpy.sign(acf_AT* acf_GC)* ( acf_AT* acf_GC)**.5
##    autocorr=(acf_AT+acf_GC)/2.
    autocorr=stattools.acf(seq_array,unbiased=unbiased,nlags=max_size, fft=True) #demean=False, fft=True)# nlags=max_size, fft=True)
##    avf = stattools.acovf(seq_array, unbiased=False, demean=True, fft=True)
##    if complex_num==True:
##        autocorr = avf[:max_size + 1] / avf[0]
##    else:autocorr = avf[:max_size*4 + 1] / avf[0]
##    norm=

##    autocorr_2=acf(autocorr,nlags=max_size, fft=True)

##    autocorr[(p_AT>.01)*(p_GC>.01)]=0.
    if complex_num==True:
        return autocorr
    else:
        return autocorr[::4]

def SeqToExanded(seq):
    char_array=numpy.fromstring(seq.upper(), '|S1')
##    char_array[char_array!='N']
##    num_array=numpy.array([0.]*len(char_array)*4)
    num_array=numpy.ndarray(((len(char_array))*4))
    num_array.fill(0.)
##    num_array.fill(0.)
    nt=['A', 'T', 'G', 'C']
    for beg in range(4):
##        print len([char_array==nt[beg]])
##        print len( num_array[beg::4])
        num_array[beg::4][char_array==nt[beg]]=1.

    return num_array
def LaggedIdentity(array, max_lag):
    output=[]
    for i in range(1,max_lag):
        output.append(numpy.mean(array[:-i]==array[i:]))
    return numpy.array(output)

def GP_background_estimator(autocorr):
    """Estimate the background autocorrelation. Assert that """

    neg_ind=numpy.where( autocorr<0)[0]
    variance_neg=autocorr[neg_ind]**2
##    pyplot.scatter(neg_ind, variance_neg**.5
##)
    rk=sklearn.model_selection.GridSearchCV( RollingKernelAverage(), {'bandwidth': numpy.logspace(0,2.5, 20)})
    rk.fit(neg_ind[:,None], variance_neg)


##    kernel=kernels.RBF()*kernels.ConstantKernel()+kernels.WhiteKernel()
##    cv=sklearn.model_selection.GridSearchCV(sklearn.kernel_ridge.KernelRidge(kernel='rbf'), {'alpha': numpy.logspace(-10,2,12)})
####    gp=sklearn.kernel_ridge.KernelRidge(kernel='rbf')
##    cv.fit(neg_ind[:,None], variance_neg)
    exp_var=rk.predict(numpy.arange(len(autocorr))[:,None])
    pyplot.plot(numpy.arange(len(autocorr)), autocorr, c='dodgerblue', alpha=.5)
    pyplot.plot(numpy.arange(len(autocorr)), (exp_var**.5)*2, c='r', lw=2)
    pyplot.show()
    filtered_array=numpy.copy(autocorr)
    filter_ind=filtered_array<=(exp_var**.5)*2
    filtered_array[filter_ind]=0.
    return filtered_array

##
##def FindBestPhase(seq, repeat_size, identity_threshold):
##    char_array=numpy.fromstring(seq.upper(), '|S1')
##    char_array=char_array[char_array!='N']
##    char_mat=ArrayToMatrix(char_array, repeat_size)
##    #Convert this into complex numbers for ease of calculation
##    seq_array=numpy.array([0.]*len(char_array), 'complex')
##    #A/T is the real plane
##    seq_array[char_array=='A']+=1
##    seq_array[char_array=='T']-=1
##
##    #G/C is the imaginary plane
##    seq_array[char_array=='G']+=1j
##    seq_array[char_array=='T']-=1j
##
##    #divide the sequence into non-overlapping kmers of the repeat size
##    #storing these as an ndarray array where axis-0 represent kmers
##    #and axis-1 represents
##    shifted_mat=numpy.roll(char_array, repeat_size)
##    pyplot.scatter(numpy.arange(len(shifted_mat)), numpy.convolve( shifted_mat==char_array,[.01]*100, 'same') )
####    pyplot.show()
##    return
##
##    #Compute percent identity between neighboring k-mers
##    kmers,length=shifted_mat.shape
##    print char_mat.shape
##    print shifted_mat.shape
##    stored_kmers=[]
##    mode_list=[]
##    nonrpt=[]
##    for k in range( kmers-2):
##        perc_ident=numpy.mean( shifted_mat[k,:]==shifted_mat[k+1,:])
##        if perc_ident>identity_threshold:
##            stored_kmers.append(char_mat[k,:])
##        else:
##            if len(stored_kmers)>0:
##                stored_kmers.append(char_mat[k,:])
##                mode=''.join(scipy.stats.mode(stored_kmers).mode[0])
##                mode_list.append(mode)
##                stored_kmers=[]
##            else:
##                nonrpt.append(char_mat[k,:])
##
##    return mode_list, ''.join(numpy.array(nonrpt).flatten())
##
##
##    best_shift=shift_list[numpy.argmax(score_list)]
##    shifted_array=numpy.roll(seq_array, best_shift)
##    shifted_mat=ArrayToMatrix(shifted_array, repeat_size)
##    return shift_list, score_list, shifted_mat
##

def Autocorr_Nulldist(read_len,reps=10000):
    nt=numpy.array( ['A','T','C','G'])
    max_autocorr=[]
    best_period=[]
    p_dist=[]
    q_dist=[]
    for r in range(reps):
        nt_ind=numpy.random.randint(0,4, size=read_len)
        seq=''.join(nt[nt_ind])
        autocorr, p, q=Autocorrel(seq,int( read_len/2.))
        p_dist.append(p)
        q_dist.append(q)
        period=numpy.argmax(autocorr[1:])+1
        max_corr=autocorr[period]
        max_autocorr.append(max_corr)
        best_period.append(period)
    return max_autocorr, p_dist, q_dist #, best_period

def HammingDist(str1,str2):
    assert len(str1)==len(str2)
    dist=0

    for i in range(len(str1)):
        dist+= str1[i]==str2[i]
    return  float(dist)/len(str1)

def TestMethod(reps=20):
    null_dist=Autocorr_Nulldist(100,10000)
    nt=['A','T','G','C']
##    cutoff=scipy.stats.scoreatpercentile(null_dist, 90)
    cutoff=-1
    seq='ATGCTCGAGGTACTAAATCGCGTAGCT'
    mut_list=[]
    rpt_len=[]
    score_list=[]
    rpt_list=[]
    for l in range(20,21):
        repeat=BoothsAlgorithm( seq[:l])
        read=(repeat*10000)#[:10000]
        for i in range(0,200,5):
            mut_ind=numpy.array( [1]*i+[0]*(1000-i))
            print len(mut_ind)
            for r in range( reps):
                numpy.random.shuffle(mut_ind)
##                read=(repeat*1000)
                read_list=list(numpy.fromstring( read,'|S1'))
                mutations=numpy.random.choice(nt, size=i)
                indel_list=numpy.where(mut_ind==1)[0]
                count=0
                read_list=list(numpy.fromstring( read,'|S1'))
                for indel in indel_list:
                    count+1
                    if scipy.stats.bernoulli(.5)==1:
                        for k in range(numpy.random.randint(1,6)):
                            read_list.insert(indel, mutations[count])
                    else:
                        for k in range(numpy.random.randint(1,6)):
                            del read_list[indel]
                read_seq=''.join(read_list[:1000])
##                mutations=numpy.random.choice(nt, size=i)
##                seq_array=numpy.fromstring(read, '|S1')
##                seq_array[mut_ind==1]=mutations
##                read_seq=''.join(seq_array )
                autocorr=Autocorrel(read_seq,50)
                period=numpy.argmax(autocorr[1:])+1
                max_corr=autocorr[period]
                mut_list.append(i/1000.)
                if max_corr<cutoff:
                    rpt_len.append(0)
                    score_list.append(0)
                    continue

                modes, non_rpt, lengths=FindBestPhase(read_seq, period,.8)
                if len(modes)==0:
                    rpt_len.append(0)

                    score_list.append(0)
                    continue
                for rpt in modes:
                    rpt_list.append(BoothsAlgorithm( rpt))
                best_repeat=modes[numpy.argmax(lengths)]
                rpt_len.append(period)

                if len(best_repeat)==len(repeat):
                    score_list.append(HammingDist(BoothsAlgorithm( best_repeat), BoothsAlgorithm( repeat)))
                else:
                    score_list.append(0)
    return mut_list, rpt_len,score_list, read, rpt_list

def MutateSequence(seq, num_mutations):
    nt=['A','T','G','C']
    mutations=numpy.random.choice(nt, size=num_mutations)
    mut_ind=numpy.array( [1]*num_mutations+[0]*(len(seq)-num_mutations))
    numpy.random.shuffle(mut_ind)
    indel_list=numpy.where(mut_ind==1)[0]
    seq_array=numpy.fromstring(seq, '|S1')
    seq_array[mut_ind==1]=mutations
    return ''.join(seq_array)

def IndelSequence(seq, rpt_len, num_mutations):
    nt=['A','T','G','C']
    rpt_units=len(seq)/rpt_len
    mutations=numpy.random.choice(nt, size=num_mutations)
    mut_ind=numpy.array( [1]*num_mutations+[0]*(rpt_units-num_mutations))
    numpy.random.shuffle(mut_ind)
##    print mut_ind
    indel_list=sorted( numpy.where(mut_ind==1)[0]*rpt_len+numpy.random.randint(0,rpt_len-1, size=num_mutations), reverse=True)
##    print indel_list
    read_list=list(numpy.fromstring(seq,'|S1'))
    for indel in indel_list:
##        count+1
        if scipy.stats.bernoulli(.5)==1:
            for k in range(numpy.random.randint(1,2)):
                read_list.insert(indel, numpy.random.choice(nt))
        else:
            for k in range(numpy.random.randint(1,2)):
                del read_list[indel]
    return ''.join(read_list)



def GetACFDecay(seq, reps):
    num_mut=[]
    best_ACF=[]
    best_period=[]
    for mut_count in range( len(seq)):
        for rep in range(reps):
            mut_seq=MutateSequence(seq, mut_count)
            auto=Autocorrel(mut_seq, len(mut_seq))
            max_period=numpy.argmax(auto[1:])+1
            max_acf=numpy.max(auto[1:])
            num_mut.append(mut_count)
            best_ACF.append(max_acf)
            best_period.append(max_period)
    return num_mut, best_ACF, best_period

def GetIndelDecay(seq, reps):
    num_mut=[]
    best_ACF=[]
    best_period=[]
    for mut_count in range( 20):
        for rep in range(reps):

            mut_seq=IndelSequence(seq,5, mut_count)

            auto=Autocorrel(mut_seq, len(mut_seq))
            max_period=numpy.argmax(auto[1:])+1
            max_acf=numpy.max(auto[1:])
            num_mut.append(mut_count)
            best_ACF.append(max_acf)
            best_period.append(max_period)
    return num_mut, best_ACF, best_period


def ArrayToMatrix(seq_array, repeat_size):

    target_size=(len(seq_array)/repeat_size)*repeat_size
    seq_matrix=numpy.reshape(seq_array[:target_size], (target_size/repeat_size, repeat_size))

    return seq_matrix






def ScoreFunction(value, cluster):
    min_cluster=float( numpy.min(cluster))
    max_cluster=float( numpy.max(cluster))
    test_1=value/min_cluster
    test_2=max_cluster/value
    return max(test_1, test_2)

def CheckOverlap(group_1, group_2):
    x_=min(group_1)
    x_prime=max(group_1)
    y_=min(group_2)
    y_prime=max(group_2)
##    print (x_, x_prime), (y_, y_prime)
    return (x_prime>=y_)*(x_<=y_prime)

def GroupOverlappingClusters(clusters):
    clusters=copy.copy(clusters)
    groups=[[clusters[0]]]
    del clusters[0]
    count=0
    while len(clusters)>0 and count<100:
        count+=1
        add_list=[]
        for i, g in enumerate(clusters):
            left_edges=[min(c.perc_id_list) for c in groups[-1]]
            right_edges=[max(c.perc_id_list) for c in groups[-1]]
            edges=left_edges+right_edges
            if CheckOverlap([min(g.perc_id_list),max(g.perc_id_list)],edges )==1:
                add_list.append(i)
        if len(add_list)>0:
            for i in add_list[::-1]:
                groups[-1].append(clusters[i])
                del clusters[i]
        if len(clusters)>0:
            groups.append([clusters[0]])
            del clusters[0]
##    for i, g in enumerate( groups):
##        for c in g:
##            print i, min(c.perc_id_list), max(c.perc_id_list), numpy.mean(c.value_list), numpy.mean(c.count_list)
    return groups


def OrganizePeriods(clusters):
    #Group clusters with overlapping intervals
    groups=GroupOverlappingClusters(clusters)


    #Within groups sort by decreasing copy number
    sorted_by_CN=[]
    for g in groups:
        sorted_by_CN.append(sorted(g, key=lambda g:numpy.sum(numpy.array( g.count_list) ), reverse=True))

    #Sort groups by decreasing percent identity
    sorted_by_PI=[]
    for g in sorted_by_CN:
##        print g
##        print numpy.mean([ c.perc_id_list for c in g ])
##        [numpy.mean( c.perc_id_list) for c in g ]
        sort_ind=numpy.argsort( [numpy.mean( c.perc_id_list) for c in g ])[::-1]
        sorted_by_PI.append([g[i] for i in sort_ind])
    #Return a list of periods and percent identity cutoffs
    percent_ident=[]
    periods=[]
    counts=[]
    for g in sorted_by_PI:
        for c in g:
            percent_ident.append(numpy.mean(c.perc_id_list) )
            periods.append(numpy.mean(c.value_list))
            counts.append(len(c.value_list))
    return zip(periods, percent_ident, counts)

class ModePath():
    def __init__(self, mode, count, perc_id, empty=False):
        if empty==False:
            self.value_dict={perc_id:[mode]}
            self.count_dict={perc_id:[count]}
        else:
            self.value_list=[mode]
            self.count_list=[count]
            self.perc_id_list=[perc_id]
    def addmode(self, mode, count, perc_id, empty=False):
        if empty==False:
            if self.value_dict.has_key(perc_id)==False:
                self.value_dict[perc_id]=[]
                self.count_dict[perc_id]=[]
            self.value_dict[perc_id].append(mode)
            self.count_dict[perc_id].append(count)
        if empty==True:
            self.value_list.append(mode)
            self.count_list.append(count)
            self.perc_id_list.append(perc_id)
    def finalize(self):
        self.perc_id_list=sorted(self.value_dict.keys())
        self.value_list=[]
        self.count_list=[]
        for ident in self.perc_id_list:
            self.value_list.append(numpy.mean(self.value_dict[ident]))
            self.count_list.append(numpy.sum(self.count_dict[ident]))
    def checkmode(self, mode):
        min_values=[min(v) for v in self.value_dict.values()]
        max_values=[max(v) for v in self.value_dict.values()]
        counts=sum([sum(v) for v in self.count_dict.values()])
        min_value=numpy.min(min_values).astype(float)
        max_value=numpy.max(max_values).astype(float)
        return max( mode/min_value, max_value/mode), counts
    def checkcount(self, count):
        score=ScoreFunction(count, self.count_list)
##        min_value=numpy.min(self.count_list).astype(float)
        return score
    def clustercounts(self, threshold=.15, min_lifespan=.03):
        """Clusters the mode path to identify intervals of percent identity
        where the number of hits is stable."""
        clusters=[]
        for i, count in enumerate( self.count_list):
            if i==0:

                clusters.append(ModePath(self.value_list[i], self.count_list[i], self.perc_id_list[i], empty=True))
##                clusters[-1].finalize()
                continue
            score=clusters[-1].checkcount(self.count_list[i])
##            print score
            if  (1-threshold<=score) and score<=(1+threshold):
                clusters[-1].addmode(self.value_list[i], self.count_list[i], self.perc_id_list[i], empty=True)
            else:
                clusters.append(ModePath(self.value_list[i], self.count_list[i], self.perc_id_list[i], empty=True))
        #Remove shortlived modes:
        final_clusters=[]
        for m in clusters:
##            m.finalize()
            lifespan=numpy.max(m.perc_id_list)-numpy.min(m.perc_id_list)
            if lifespan>=min_lifespan:
                final_clusters.append(m)

        return final_clusters

def PlotDistances(ccf, min_identity=.8):
    for perc_id in numpy.arange(.995, min_identity,-.005):
        hits=numpy.where(ccf>=perc_id)[0]
        distances=numpy.diff(hits)
        pyplot.scatter(distances, [perc_id]*len(distances), alpha=.1, s=10,c='black')
    pyplot.ylabel('Percent identity')
    pyplot.xlabel('Periodicity')
    pyplot.tick_params(length=4, direction='in')
    pyplot.show()

def ConstructModeTree(sig, min_identity, cluster_threshold=.1, plot=False ):
    mode_list=[]
    terminate=False
    for perc_id in numpy.arange(.995, min_identity,-.005):
        print '{0}, '.format(perc_id),
        #Cluster the distances between hits at the %ident threshold to find
        #periods
        modes=IdentifyPeriods(sig, perc_id)
        if len(modes)==0: continue
        try:
            modes, counts=zip(*modes)
        except:
            print modes
            print perc_id
            print jabber
            continue

        for i, mode in enumerate( modes):
            new_cluster=True
            for m in mode_list:
                score, count= m.checkmode(mode)
                if count>=5000: terminate=True
                if score<1+cluster_threshold  and score>1-cluster_threshold :
                    new_cluster=False
                    m.addmode(mode, counts[i], perc_id)
            if new_cluster==True:
                mode_list.append(ModePath(mode, counts[i], perc_id ))
        if terminate==True: break
    for m in mode_list:
        m.finalize()
##        mode_list.append(modes)
##        print len(modes)
##        pyplot.scatter(modes, [perc_id]*len(modes), s=numpy.array( counts )/1.)
    if plot==True: PlotModeTree(mode_list)
    mode_list=ClusterModeTree(mode_list)
    if plot==True: PlotModeTree(mode_list)
    return mode_list


def IdentifyPeriods(signal, threshold=.95):
    hits=numpy.where(signal>=threshold)[0]
    if len(hits)==0:
        return []
    distance=numpy.diff(hits).astype(float)
    if len( distance)==0:
        return []
    clusters=HierarchicalCluster(distance,.1)

    periods=[( numpy.mean(c), len(c)) for c in clusters if len(c)>4]

    return periods

def HierarchicalCluster(lengths, threshold, plot=False):
    """This takes a list of lengths and clusters them as follows:
        For each adjacent pair, it computes the  """

    #Sort the lengths and sequences by length in ascending order
    sort_ind=numpy.argsort(lengths)
    lengths=numpy.array( lengths) [sort_ind].astype(float)

    distances=abs(numpy.diff(lengths))
##    if plot==True:
##        pyplot.plot(lengths)
##        pyplot.show()
##        print jabber
    midpoints=(lengths[1:]+lengths[:-1])/2.
    score=distances/(2*midpoints)
##    print score
##    print score*100
    clusters=[[lengths[0]]]

    for i,s in enumerate(lengths[1:]):
##        score=s/ min(clusters[-1])
        score=ScoreFunction(s, clusters[-1])
        if score<=(1.+threshold): clusters[-1].append(lengths[ i])
        else:
##            print score
            clusters.append([lengths[ i]])
    if plot==True:
        for c in clusters:
            pyplot.plot(c)
        pyplot.show()
##        print jabber

    return clusters


def ClusterModeTree(modetree):
    newtree=[]
    for modepath in modetree:
##        modepath.finalize()
        newtree+=modepath.clustercounts()
    return newtree

def PlotModeTree(modetree):
    colors=matplotlib.colors.cnames.keys()
    numpy.random.shuffle(colors)
    for i, branch in enumerate( modetree):
        pyplot.scatter(branch. value_list, branch.perc_id_list, s=numpy.array( branch.count_list).astype(float), c=colors[i])
    pyplot.ylabel('Percent Identity', size=10)
    pyplot.xlabel('Periodicity', size=10)
    pyplot.tick_params(length=4, direction='in')
    pyplot.show()

    for i, branch in enumerate( modetree):
        pyplot.plot(branch. perc_id_list, branch.count_list, c=colors[i])
    pyplot.show()
def ClusterNeighbors():
    pass

def AnalyzeRead(read):
    period_list=[]
    for j in range(1):

        if len(read)<10: break
        try:
            autocorr=Autocorrel(read,len(read), False)
        except:
            break
        pyplot.plot(autocorr)
        pyplot.show()
        pyplot.close()
        true_period=numpy.argmax(autocorr[1:])+1
        print true_period

        #Only interested in base periodicities>3
        #Remove repeats of this length from the sequence
        rpts=ExtractRepeatsOfSize(read, true_period,.6)
        for r,v in rpts:
            ccf=LaggedUngappedAlignment(r, read)
##            return ccf
            pyplot.plot(ccf)
            pyplot.show()

            break

    return repeats



def main(argv):
    print argv
    param={}
    for i in range(1, len(argv), 2):
        param[argv[i]]= argv[i+1]
    print param
    if param=={}: return
    TandemFinder(param['-i'], param['-o'], param['-m'] ,.8)

if __name__ == '__main__':
    main(sys.argv)
