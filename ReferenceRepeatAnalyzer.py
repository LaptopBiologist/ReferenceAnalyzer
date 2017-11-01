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
import matplotlib
##matplotlib.use('Agg')
from matplotlib import pyplot
from statsmodels.tsa import stattools
from statsmodels.tsa.stattools import acf
import sklearn
import sklearn.gaussian_process
import sklearn.kernel_ridge
import sklearn.metrics

import sklearn.gaussian_process.kernels as kernels
from sklearn.base import BaseEstimator
from Bio import Seq
from matplotlib import pyplot
import seaborn
import itertools

import gzip

import sys

import copy

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


def CountSequences(k=6, alphabet='ATGC'):
    """Enumerate all possible k-length products of the alphabet and then
    removes redundancies. Returns the cardinality."""
    sequences=itertools.product(alphabet, repeat=k)
    return list(sequences)

def ComputeKmerCompositionEntropy(sequence,k=5):
    """Summarize the complexity of a repeat unit by decomposing it into kmers
    and then modeling the expected counts of each observed kmer with a multinomial distribution.
    The multinomial probability is multiplied by the number of possible sets of kmers
    the same size as the observed set (the binomial coefficient), because we don't care at all about the
    identity of the kmers. This kept in log-space to avoid precision errors.
    The log-probability is divided by the number of kmers the sequence was decomposed into.

    Note: the binom coefficient becomes imprecise for k>5; 5-mers though provides
    a reasonable summary of sequece complexity."""

    #Decompose the sequence into kmer counts
    kmer_counts=CountKMERS(sequence, k)

    #Number of kmers possible
    num_poss_kmers=(4.**k)

    #Assume each kmer is equally likely
    pr=1./num_poss_kmers

    #The number of kmers contained in the sequence
    num_kmers=sum(kmer_counts.values())

    obs=kmer_counts.values()+[0]* (int(num_poss_kmers-len(kmer_counts.keys())))

    #Multiply the probablity by the binomial coefficient: Don't care which kmers are enriched!
    #In log space, this means add the logs
    logprob=numpy.log( scipy.special.binom(num_poss_kmers, len(kmer_counts.keys()))) + scipy.stats.multinomial.logpmf(obs, num_kmers, [pr]*int( num_poss_kmers))

    information=0

    if logprob!=0:
        information=-1*logprob
    else:
        information=numpy.inf

    #Output average information per kmer
    return information/(num_kmers)

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


def SplitFastaByEntropy(infile, outfile):
    low_outhandle=open(outfile+'_low.fa', 'w')
    high_outhandle=open(outfile+'_high.fa', 'w')
    seqs=GetSeq(infile, upper=True)
    for key in seqs.keys():
        entropy=ComputeKmerCompositionEntropy(seqs[key],5)
        if entropy<2:
            low_outhandle.write('>{0}_entropy={1}\n'.format(key, entropy))
            low_outhandle.write('{0}\n'.format(seqs[key]))
        else:
            high_outhandle.write('>{0}_entropy={1}\n'.format(key, entropy))
            high_outhandle.write('{0}\n'.format(seqs[key]))
    high_outhandle.close()
    low_outhandle.close()
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


def CountKMERS(sequence, k=10 ):
##    complement=str(Seq.Seq( sequence).reverse_complement())
    kmerSet={}
    for x in range(0,len(sequence)-k+1):
        kmer=str(sequence[x:x+k])
        if kmerSet.has_key(kmer)==False:
            kmerSet[kmer]=0.
        kmerSet[kmer]+=1
    return kmerSet
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
    print len(seq)

    period_list=[]
    window_list=[]
    corr_list=[]
    mode_list={}
    len_dict={}
    #Test over range of window sizes: 1000, 10000, 100000, 1000000
    for window_size in [10000,100000,500000]:
##        null_dist=Autocorr_Nulldist(window_size)
        step_size=window_size/2
        max_repeat=window_size/2

        window_edges=numpy.arange(0,len(seq)-window_size+step_size, step_size)
        for i in range(len(window_edges)-1) :
    ##        print i, window_edges[i],window_edges[i+1]
            seq_slice=seq[window_edges[i]:window_edges[i]+window_size]
            for j in range(5):
                if len(seq_slice)<100: break
                try:
                    autocorr=Autocorrel(seq_slice,max_repeat)
                except:
                    break
        ##        print len(autocorr)

                true_period=numpy.argmax(autocorr[1:])+1

                #Only interested in base periodicities>3
                if true_period>3:
                    #Only going to look for periodicities>20
                    period=numpy.argmax(autocorr[20:])+20

    ##                corr_list.append(autocorr[period-1])
                    modes, seq_slice,len_list= FindBestPhase(seq_slice, period, .85)

                    if len(modes)>0:
##                        print 'Minlen={0}'.format( min(len_list))
                        window_list.append(window_edges[i])
                        period_list.append(period)

                        mode_candidates=[]
                        for k in range(len(modes)):
                            m=BoothsAlgorithm( modes[k])
                            mode_candidates.append(m)
                            if len_dict.has_key(m)==False:
                                len_dict[m]={}
                            if len_dict[m].has_key(window_size)==False:
                                len_dict[m][window_size]=0.
                            len_dict[m][window_size]+=len_list[k]
                        for k in range(len(mode_candidates)):
                            m=mode_candidates[k]
                            if mode_list.has_key(true_period)==False:
                                mode_list[true_period]={}
                            if mode_list[true_period].has_key(len(m))==False:
                                mode_list[true_period][len(m)]= set()

                            mode_list[true_period][len(m)].add(m )

                            #Remember how many of instances of this repeat were identified


                else: break

    #Cluster modes
##    for length in mode_list.keys():
##        modes=list(mode_list[length])
##        jac_array=numpy.ndarray((len(modes), len(modes)))
##        jac_array.fill(0.)
##        for i in range(len(modes)):
##            for j in range(i):
##                kmer_i=set(GetKMERS(modes[i]))
##                kmer_j=set(GetKMERS(modes[j]))
##                jac_array[i,j]=float(len(kmer_i&kmer_j))/len(kmer_i|kmer_j)
##        seaborn.heatmap(jac_array)
##        pyplot.show()

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
        try:
            counts=float(key.split('_')[-1].split('=')[-1])
            chrom=key.split('_')[0]
        except:
            print key
            print jab
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
    print len(sequences)
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


def Autocorrel(seq, max_size, complex_num=True, p=False):
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
    autocorr=stattools.acf(seq_array,unbiased=False,nlags=max_size, fft=True) #demean=False, fft=True)# nlags=max_size, fft=True)
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

def ClusterNeighbors():
    pass


def main(argv):
    print argv
    param={}
    for i in range(1, len(argv), 2):
        param[argv[i]]= argv[i+1]
    print param
    if param=={}: return
    ExtractTandems(param['-i'], param['-o'])

if __name__ == '__main__':
    main(sys.argv)