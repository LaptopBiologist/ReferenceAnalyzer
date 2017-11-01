#-------------------------------------------------------------------------------
# Name:        InternalRptFinder
# Purpose:
#
# Author:      Michael mcGurk
#
# Created:     20/05/2014
# Copyright:   (c) Michael McGurk 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import subprocess
import shutil
import os
import sys
import csv
import networkx
import community
import numpy

import seaborn
import matplotlib
from matplotlib import pyplot

import Bio
from Bio import SearchIO
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord


class alignment():
    def __init__(self, row):
        self.query=row[0]
        self.subject=row[1]
        self.identity=float(row[2])
        self.length=int(row[3])
        self.covs=float(row[-2])
        self.pair=tuple(sorted((self.query, self.subject)))

def MakeDir(newdir):
    if os.path.exists(newdir)==False:
        os.mkdir(newdir)

def GetSeq(ref, parser='fasta', gzipped=False, clean_name=True):
    if gzipped==False:
        handle=open(ref, 'r')
    else:
        handle=gzip.open(ref)
    lib=SeqIO.parse(handle, parser)
    SeqLen={}
    for rec in lib:
        #if refName.count(rec.name)==0: continue
        if clean_name==True:
            SeqLen[CleanName(rec.name)]=rec
        else:
            SeqLen[rec.name]=rec
    handle.close()
    return SeqLen

def CleanName(name):
    illegal=['|', '!','>', '<', '?','/','*']
    for i in illegal:
        name='_'.join(name.split(i))
    return(name)


def ClusterIndexByLength(index, outDir, BlastDir, cutoff=85):
    #Split file by lengths
    inDir='/'.join(index.split('/')[:-1])
    tempDir=inDir+'/temp'
    SplitFileByLength(index, tempDir)
    file_list=os.listdir(tempDir)
    for f in file_list:
        infile=tempDir+'/'+f
        outfile=tempDir+'/'+f.split('.')[0]+'_out'
        ClusterIndex(infile, outfile, BlastDir, cutoff=80)


def ClusterIndex(index, outDir, BlastDir, cutoff=85):

    #Split the index into files containing 300 seq.
    inDir='/'.join(index.split('/')[:-1])
    tempDir=inDir+'/temp'
    splitFile(index, tempDir, 300)

    outfiles=[]

    fileList=os.listdir(tempDir)

    for j in range(len(fileList)):
        for i in range(j+1):
            print fileList[j], fileList[i]
            outname='alignments_{0}_{1}.tsv'.format(str(i), str(j))
            outfiles.append( BlastSeq_part(tempDir+'/'+ fileList[j],tempDir+'/'+ fileList[i], outDir, outname, BlastDir))
    AlignmentFile=outDir+'/alignment_output.tsv'
    JoinFiles(outfiles, outDir, AlignmentFile)


    aligned=ParseSequences(AlignmentFile)
    graph=BuildGraph(aligned, cutoff)
    clusters=GraphToCommunities(graph)
    groups, singletons=RemoveSingletons(clusters)


    WriteOutput(groups,singletons, outDir, index)
    shutil.rmtree(tempDir)
    return groups


def JoinFiles(infiles, inDir, outfile):
    outhandle=open(outfile, 'w')
    for f in infiles:
        inhandle=open(f,'r')

        for line in inhandle:
            outhandle.write(line)
        inhandle.close()
    outhandle.close()


def SplitFileByLength(infile, tempdir):
    MakeDir(tempdir)
    refSeq=GetSeq(infile, clean_name=False)
    #
    length_dict={}

    for read_name in refSeq.keys():
        length=float( read_name.split('_')[1].split('=')[1])
        if length_dict.has_key(length)==False:
            length_dict[length]={}
        length_dict[length][read_name]=refSeq[read_name]

    inDir='/'.join(infile.split('/')[:-1])
    root=infile.split('/')[-1]
    rootname='.'.join(root.split('.')[:-1])
    ext=root.split('.')[-1]
    for length in length_dict.keys():
        outfile="{0}/{1}_temp_{2}.{3}".format(tempdir, rootname,int( length), ext)
        outHandle=open(outfile, 'w')
        for read_name in length_dict[length].keys():
            outHandle.write('>'+read_name+'\n')
            outHandle.write(str(length_dict[length][read_name].seq)+'\n')

        outHandle.close()

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
        seq=str( refSeq[key].seq)
        name= CleanName( str(refSeq[key].name))

        outHandle.write('>'+name+'\n')
        outHandle.write(seq+'\n')
        count+=1

    outHandle.close()

def BlastSeq(Query, Subject, Out, BlastDir):
    """BLASTs two fastas against each other."""
    MakeDir(Out)
    OutPath='.'.join(Out.split('.'))
    OutFile=OutPath+'/output.csv'
    print (OutPath)
    errlog=open(OutPath+'/_err.log', 'a')
    column_spec='10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue btop'
    BLAST=subprocess.Popen([BlastDir+'/bin/blastn', '-query',Query, '-subject',Subject, '-outfmt', column_spec,  '-out', OutFile], stderr=errlog)
    BLAST.communicate()
    errlog.close()
    return OutFile

def BlastSeq_part(Query, Subject, OutPath, outname, BlastDir):
    """BLASTs two fastas against each other."""
    MakeDir(OutPath)
    OutFile=OutPath+'/'+outname
    print (OutPath)
    errlog=open(OutPath+'/_err.log', 'a')
    column_spec='10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qcovs btop'
    BLAST=subprocess.Popen([BlastDir+'/bin/blastn', '-query',Query, '-subject',Subject, '-outfmt', column_spec,  '-out', OutFile], stderr=errlog)
    BLAST.communicate()
    errlog.close()
    return OutFile


def BlastDir(insDir, consDir, outDir):

    insFile=os.listdir(insDir)
    consFile=os.listdir(consDir)
    MakeDir(outDir)
    for f in insFile:
        print "BLASTing {0}...".format(insFile)
        consPath=consDir+'/'+f
        insPath=insDir+'/'+f
        if os.path.exists(consPath)==False: continue
        BlastSeqII(insPath, consPath, outDir,  '.'.join(f.split('.')[:-1]) , 'c:/ncbi-blast')

def BlastDirII(insDir, consPath, outDir):

    insFile=os.listdir(insDir)
    MakeDir(outDir)
    for f in insFile:
        print "BLASTing {0}...".format(f)
        insPath=insDir+'/'+f
        BlastSeqII(insPath, consPath, outDir,  '.'.join(f.split('.')[:-1]) , 'c:/ncbi-blast')


def BlastSeqII(Query, Subject, Out, name, BlastDir):
    """BLASTs two fastas against each other."""
    OutPath='.'.join(Out.split('.'))
    OutFile=OutPath+'/'+name+'.csv'
    print (OutPath)
    errlog=open(OutPath+'/_err.log', 'a')
    column_spec='10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue btop'
    BLAST=subprocess.Popen([BlastDir+'/bin/blastn', '-query',Query, '-subject',Subject, '-outfmt', column_spec,  '-out', OutFile], stderr=errlog)
    BLAST.communicate()
    errlog.close()
    return OutFile


def ParseSequences(InFile):
    """Organizes the BLAST output and stores the information in python's working
     memory."""
    handle=open(InFile)
    table=csv.reader(handle )
    pairDict={}
    for row in table:
        aligned=alignment(row)
        if pairDict.has_key(aligned.pair)==False:
            pairDict[aligned.pair]=aligned
        if aligned.length>pairDict[aligned.pair]:pairDict[ aligned.pair]=aligned
    return(pairDict)


def BuildGraph(aligned, cutoff):
    """Takes the BLAST output and represents it as a graph using the following rules:
        Each sequence is represented as a node.
        Each alignment between sequences is represented as an unweighted,
        undirected edge."""
    #hold=set()
    #hold2=set()
    graph=networkx.Graph()
    for pair in aligned.keys():
        query,subject=pair
        if graph.has_node(subject)==False:
            graph.add_node(subject)
        if graph.has_node(query)==False:
            graph.add_node(query)
        if aligned[pair].identity>=cutoff and aligned[pair].covs>80: #aligned[pair].length>100:
            graph.add_edge(subject, query)

    return graph

def GraphToCommunities(Network):
    """Identifies communities in the graph using the Louvain method as implemented in community."""
    comm=community.best_partition(Network)
    clusters={}
    for k in comm.keys():
        if clusters.has_key(comm[k])==False:
            clusters[comm[k]]=[]
        clusters[comm[k]].append(k)
    return (clusters)

def RemoveSingletons(clusters):
    """Removes clusters containing a single node. That is, entries that align
    only to themselves."""
    groups=[]
    singletons=[]
    for key in clusters.keys():
        if len(clusters[key])>1: groups.append(clusters[key])
        else: singletons+=clusters[key]
    return groups, singletons

def WriteOutput(clusters, singletons, outDir, index):
    MakeDir(outDir)
    inhandle=open(index,'r')
    lib=SeqIO.parse(inhandle, 'fasta')
    hold={}

    for book in lib:
        name=book.id.strip()
        name=name.split('\t')[0]
##        if name[-1]=='.': name=name[:-1]
        hold[CleanName( name)]=book

    for h in hold.keys():
        print h

    outfile=outDir+'/communities.txt'
    outHandle=open(outfile, 'w')
    count=0

    for community in clusters:
##        if len(community)==1: continue
        dump=outHandle.write('Community {0}, Members = {1}:\n\n'.format(count, len(community)))
        outFasta=outDir+'/community_{0}.fa'.format(count)
        FASTAHandle=open(outFasta, 'w')
        sequences=[]
        for member in community:
            sequences.append(hold[member])
            dump=outHandle.write('\t{0}, '.format(member))
        dump=outHandle.write('\n\n')
        dump=SeqIO.write(sequences, FASTAHandle, 'fasta')
        FASTAHandle.close()
        count+=1

    #Write singletons deal with later
    dump=outHandle.write('Singletons, Members = {0}:\n\n'.format( len(singletons)))
    outFasta=outDir+'/singletons.fa'.format(count)
    FASTAHandle=open(outFasta, 'w')
    sequences=[]
    for member in singletons:
        sequences.append(hold[member])
        dump=outHandle.write('\t{0}, '.format(member))
    dump=outHandle.write('\n\n')
    dump=SeqIO.write(sequences, FASTAHandle, 'fasta')
    FASTAHandle.close()

    outHandle.close()
    inhandle.close()
    return(hold)


def ClusterByLengths(indir, outfile):
    file_list=os.listdir(indir)
    outhandle=open(outfile, 'w')
    for f in file_list:
        if f=='singletons.fa':
            sequences=GetSeq(indir+'/'+f, clean_name=False)
            for key in sequences.keys():
                outhandle.write('>{0}\n'.format(key))
                outhandle.write('{0}\n'.format(str(sequences[key].seq)))
            continue
        if f.split('_')[0]!='community':continue
        sequences=GetSeq(indir+'/'+f, clean_name=False)
        hc, clusters=HierarchicalCluster(sequences)
        for key in clusters:
            outhandle.write('>{0}\n'.format(key))
            outhandle.write('{0}\n'.format(str(sequences[key].seq)))
    outhandle.close()

def FilterByLengths(infile, outfile, cutoff=100):
    sequences=GetSeq(infile, clean_name=False)
    outhandle=open(outfile, 'w')
    for key in sequences:
        length=len(sequences[key])
        if length<cutoff: continue
        outhandle.write('>{0}\n'.format(key))
        outhandle.write('{0}\n'.format(str(sequences[key].seq)))
    outhandle.close()
def HierarchicalCluster(sequences):
    """This takes a list of lengths and clusters them as follows:
        For each adjacent pair, it computes the  """
    seq_keys=numpy.array(sequences.keys())
    lengths=numpy.array( [len(sequences[k]) for k in seq_keys], float)


    #Sort the lengths and sequences by length in ascending order
    sort_ind=numpy.argsort(lengths)
    lengths=lengths[sort_ind]
    seq_keys=seq_keys[sort_ind]




    distances=abs(numpy.diff(lengths))
    midpoints=(lengths[1:]+lengths[:-1])/2
    score=distances/(2*midpoints)
    clusters=[[seq_keys[0]]]
    for i,s in enumerate(score):
        if s<.02: clusters[-1].append(seq_keys[i+1])
        else:clusters.append([seq_keys[i+1]])

    retained_keys=[]
    for c in clusters:
        try:
            counts=numpy.array( [float(k.split('_')[-1].split('=')[-1]) for k in c])
        except:
            print k
            print jabber
        best_ind=numpy.argmax(counts)
        retained_keys.append(c[best_ind])

    return clusters, retained_keys

def PlotClusters(clusters):
    count=0
    colors=matplotlib.colors.cnames.keys()
    numpy.random.shuffle(colors)
    for c in clusters:

        counts=numpy.array( [float(k.split('_')[2].split('=')[-1]) for k in c])
        lengths=numpy.array( [float(k.split('_')[1].split('=')[-1]) for k in c])
        t=pyplot.scatter(lengths, numpy.log10( counts),c=colors[count%len(colors)] )

        count+=1
    pyplot.show()


def main(argv):
    param={}
    print argv
    for i in range(1, len(argv), 2):
        param[argv[i]]= argv[i+1]
    print param
    if param=={}: return()
    #inDir, outDir, AssemblyIndex, TEIndex, threads, length=70, ReadIDPosition=1, PositionRange=0, phred=33, shorten=True, Trim=True
    inDir=param['-i']
    outDir=param['-o']
    blast_path=param['-B']
    if param.has_key('-C')==True:
        cutoff=float( param['-C'])
    else:
        cutoff=.8
    ClusterIndex(inDir, outDir, blast_path, cutoff)
    pass

if __name__ == '__main__':
    main(sys.argv)

