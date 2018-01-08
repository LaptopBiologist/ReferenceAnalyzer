#-------------------------------------------------------------------------------
# Name:        module1
# Purpose: This contain functions for two purposes
# 1) Annotate tandem repeats with the features they span
# 2) Extract tandem repeats based on known coordinates
#
# Author:      I am
#
# Created:     06/12/2017
# Copyright:   (c) I am 2017
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy
import HTSeq
import csv
import pyliftover
import sys
##from ReferenceAnalyzer.ReferenceRepeatAnalyzer import *


class TandemLine():
    def __init__(self, row):
        self.id=row[0]
        self.method=row[1]
        self.chrom=row[2]
        self.rpt_length=int(row[3])
        self.copynumber=int(row[4])
        self.info=float(row[5])
        self.identity=float(row[6])
        try:self.repetition=float(row[7])
        except: self.repetition='err'
        self.start=int(row[8])
        self.end=int(row[9])
        self.starts=[int(s) for s in row[10].split(',') ]
        self.ends=[int(e) for e in row[11].split(',') ]
        self.info_string=row[12]
        self.seq=row[13]
        try:
            gene_parts= [item.split(':') for item in row[14].split(',')]
            self.gene=dict(gene_parts)
        except:
            self.gene={}
        try:
            rpt_parts= [item.split(':') for item in row[15].split(',')]
            self.repeat=dict(rpt_parts)
        except:

            self.repeat={}


class RefseqLine():
    def __init__(self, row):
##        pass
        self.bin=int(row[0])
        self.name=row[1]
        self.chrom=row[2]
        self.strand=row[3]
        self.txStart=int(row[4])
        self.txEnd=int( row[5])
        self.cdsStart=int( row[6])
        self.cdsEnd=int( row[7])
        self.exonCount=int( row[8])
        self.exonStarts=[int(i) for i in  row[9].split(',')[:-1]]
        self.exonEnds=[int(i) for i in  row[10].split(',')[:-1]]
        self.score=row[11]
        self.name2=row[12]
        self.cdsStartStat=row[13]
        self.cdsEndStat=row[14]
        self.exonFrames=[int(i) for i in  row[15].split(',')[:-1]]

class WarburtonLine():
    def __init__(self, row):
##        pass
        self.bin=int(row[0])
        self.name=row[1]
        self.chrom=row[2]
        self.strand=row[3]
        self.txStart=int(row[4])
        self.txEnd=int( row[5])
        self.cdsStart=int( row[6])
        self.cdsEnd=int( row[7])
        self.exonCount=int( row[8])
        self.exonStarts=[int(i) for i in  row[9].split(',')[:-1]]
        self.exonEnds=[int(i) for i in  row[10].split(',')[:-1]]
        self.score=row[11]
        self.name2=row[12]
        self.cdsStartStat=row[13]
        self.cdsEndStat=row[14]
        self.exonFrames=[int(i) for i in  row[15].split(',')[:-1]]

class RepeatLine():
    def __init__(self, row):
        pass

        self.chrom=row[5]
        self.start= int(row[6])
        self.end= int(row[7])
        self.name=row[10]



def LoadRefseqAnnotations(infile):
    gen_array=HTSeq.GenomicArrayOfSets('auto', False)
    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter='\t')
    header=intable.next()
    for row in intable:

        line=RefseqLine(row)
        for i in range(len( line.exonStarts)):
            gene_iv=HTSeq.GenomicInterval(line.chrom, line.exonStarts[i], line.exonEnds[i])
            gen_array[gene_iv]+=line.name2


    inhandle.close()
    return gen_array

def LoadWarburtonAnnotations(infile):
    gen_array=HTSeq.GenomicArrayOfSets('auto', False)
    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter='\t')
    header=intable.next()
    warburton_set=set()
    for row in intable:

        chrom,interval=row[-1].split(':')
##        chrom='Cchrom[1:]
        left, right=interval.split('-')
        try:
            gene_iv=HTSeq.GenomicInterval(chrom,int( left),int( right))
        except: continue


        warburton_set.add(row[-3])
        gen_array[gene_iv]+=row[-3]


    inhandle.close()
    return gen_array, warburton_set

def LoadRepeatAnnotations(infile):
    gen_array=HTSeq.GenomicArrayOfSets('auto', False)
    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter='\t')
    header=intable.next()
    for row in intable:

        line=RepeatLine(row)

        gene_iv=HTSeq.GenomicInterval(line.chrom, line.start, line.end)
        gen_array[gene_iv]+=line.name


    inhandle.close()
    return gen_array

def AnnotateWithRefSeq(infile, outfile, annot_file, rpt_file):
    print 'Loading gene annotations...'
    gen_array=LoadRefseqAnnotations(annot_file)
    print 'Loading repeat annotations...'
    rpt_array=LoadRepeatAnnotations(rpt_file)
    print 'Annotations loaded'
    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter='\t')
    outhandle=open(outfile, 'w')
    outtable=csv.writer(outhandle, delimiter='\t')
    header=intable.next()
    outtable.writerow(header+['Features'])
    counter=0
    for row in intable:
        if counter%100==0:
            print '.',
        counter+=1
        line=TandemLine(row)
        gene_annotations={}
        rpt_annotations={}
        for i in range(len( line.starts)):
            s,e=line.starts[i], line.ends[i]

            tandem_iv=HTSeq.GenomicInterval(line.chrom, s, e)
            for gen_annot in set.union(* [ seq[1] for seq in gen_array[tandem_iv].steps() ]):
##                for gen_annot in set( sets):
##                    print gen_annot
                if gene_annotations.has_key(gen_annot)==False:
                    gene_annotations[gen_annot]=0
                gene_annotations[gen_annot]+=1
            for gen_annot in set.union(* [ seq[1] for seq in rpt_array[tandem_iv].steps() ]):

                if rpt_annotations.has_key(gen_annot)==False:
                    rpt_annotations[gen_annot]=0
                rpt_annotations[gen_annot]+=1

##        gene_annotations=','.join(list(gene_annotations))
##        rpt_annotations=','.join(list(rpt_annotations))
        gene_annotations=','.join(['{0}:{1}'.format(key, gene_annotations[key]) for key in gene_annotations ])
        rpt_annotations=','.join(['{0}:{1}'.format(key, rpt_annotations[key]) for key in rpt_annotations ])
##        rpt_annotations=''
        outtable.writerow(row+[gene_annotations, rpt_annotations])
    outhandle.close()
    inhandle.close()

def LoadRepeatMaskerAnnotations():
    pass



def UpdateWarburtonTable1(infile, ref_file,outfile):
    inhandle=open(infile,'r')
    intable=csv.reader(inhandle, delimiter='\t')

    outhandle=open(outfile, 'w')

    updated_file='{0}_hg38.tsv'.format( '.'.join( infile.split('.')[:-1]))
    updated_handle=open(updated_file,'w')
    updated_table=csv.writer(updated_handle, delimiter='\t')
    header=intable.next()
    updated_table.writerow(header)
    lo=pyliftover.LiftOver('hg18', 'hg38')


##    seq=GetSeq(ref_file)"
    for row in intable:
        chrom, interval=row[-1].split(':')
        left,right=interval.split('-')
        left=int(''.join(left.split(',')))
        right=int(''.join(right.split(',')))

        coord_left=lo.convert_coordinate(chrom,left)[0][1]
        chromosome,coord_right=lo.convert_coordinate(chrom,right)[0][:2]
        print chromosome, left, coord_left
        print chromosome, right, coord_right

        new_line=row[:-1]+['{0}:{1}-{2}'.format(chromosome, coord_left, coord_right)]
##        seq_name='>{0}_{1}_{2}_Up{3}_{4}_{5}\n',format(row[7],row[2],row[3], chromosome, coord_left, coord_right,)
##        outfile.write(seq_name)
##        outfile.write ( '{0}\n'.format( seq[chromosome][coord_left:coord_right].upper() ) )
        updated_table.writerow(new_line)
    inhandle.close()
    updated_handle.close()
    outfile.close()

def AnnotateWarburton(infile,annot_file):
    print 'Loading gene annotations...'
    gen_array, gen_set=LoadWarburtonAnnotations(annot_file)
##    print gen_set
    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter='\t')
##    outhandle=open(outfile, 'w')
##    outtable=csv.writer(outhandle, delimiter='\t')
    header=intable.next()
##    outtable.writerow(header+['Features'])
    counter=0
    found_set=set()
    for row in intable:
##        if counter%100==0:
##            print '.',
        counter+=1
        line=TandemLine(row)
        gene_annotations=set()
        rpt_annotations=set()
        for i in range(len( line.starts)):
            s,e = line.starts[i], line.ends[i]

            tandem_iv = HTSeq.GenomicInterval(line.chrom, s, e)
            gene_annotations |= set.union(* [ seq[1] for seq in gen_array[tandem_iv].steps()])
##            if len( gene_annotations) > 0 : print '.'
            found_set |= gene_annotations
##        gene'_annotations=','.join(list(gene_annotations))

##        rpt_annotations=''

##    outhandle.close()
    print gen_set-found_set
    print '{0}/{1}: {2}'.format(len(found_set), len(gen_set), float(len(found_set))/ len(gen_set))
    inhandle.close()

def ConcatenateAnnotations(file1, file2, outfile):
    handle1=open(file1, 'r')
    handle2=open(file2, 'r')
    table1=csv.reader(handle1, delimiter='\t')
    table2=csv.reader(handle2, delimiter='\t')
    outhandle=open(outfile, 'w+')
    outtable=csv.writer(outhandle, delimiter='\t')
    header=['Repeat ID','Settings']+ table1.next()
    outtable.writerow(header)
    count=0
    for line in table1:
        outtable.writerow([count, '80']+line)
        count+=1
    header=table2.next()
    for line in table2:
        outtable.writerow([count, 'auto']+line)
        count+=1
    handle1.close()
    handle2.close()
    outhandle.close()

def ExtractFASTAs(annot_file, outfile):
    #
    inhandle=open(annot_file, 'r')
    intable=csv.reader(inhandle, delimiter='\t')

    gene_file='{0}_genes.fa'.format(outfile)
    alpha_file='{0}_alpha.fa'.format(outfile)
    rpt_file='{0}_rpt.fa'.format(outfile)

    gene_handle=open(gene_file, 'w')
    rpt_handle=open(rpt_file, 'w')
    alpha_handle=open(alpha_file, 'w')
    intable.next()
    count=0
    for row in intable:
        line=TandemLine(row)
        if len(line.repeat.keys())>0: count+=1
        if count%20==0:
            print line.repeat
            count+=1
        if len( line.gene.keys()) !=0:
            keys,values=line.gene.keys(),numpy.array(  line.gene.values(), int)
            name=keys[numpy.argmax(values)]
            seq_name='>{0}_Name={1}_CN={2}_Length={3}_Identity={4}_Rpt={5}_method={6}_Chr={7}:{8}-{9}\n'.format(line.id,\
            name, line.copynumber, line.rpt_length, line.identity, line.repetition, line.method, line.chrom,\
            '{:,}'.format( line.start), '{:,}'.format( line.start))
            gene_handle.write(seq_name)
            gene_handle.write('{0}\n'.format(line.seq))
        elif line.repeat.has_key('ALR/Alpha'):
            name='ALR/Alpha'
            seq_name='>{0}_CN={2}_Length={3}_Identity={4}_Rpt={5}_method={6}_Name={1}_Chr={7}:{8}-{9}\n'.format(line.id,\
            name, line.copynumber, line.rpt_length, line.identity, line.repetition, line.method, line.chrom,\
            '{:,}'.format( line.start), '{:,}'.format( line.start))
            alpha_handle.write(seq_name)
            alpha_handle.write('{0}\n'.format(line.seq))
        else:
            if len( line.repeat.keys()) ==0: name='Unannot'
            elif len(line.repeat.keys())==1: name=line.repeat.keys()[0]
            else: name='Complex'
            seq_name='>{0}_Name={1}_CN={2}_Length={3}_Identity={4}_Rpt={5}_method={6}_Chr={7}:{8}-{9}\n'.format(line.id,\
            name, line.copynumber, line.rpt_length, line.identity, line.repetition, line.method, line.chrom,\
            '{:,}'.format( line.start), '{:,}'.format( line.start))
            rpt_handle.write(seq_name)
            rpt_handle.write('{0}\n'.format(line.seq))

    inhandle.close()
    gene_handle.close()
    rpt_handle.close()
    alpha_handle.close()
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
    gene_path=param['-gene']
    rpt_path=param['-rpt']
    if param.has_key('-C')==True:
        cutoff=float( param['-C'])
    else:
        cutoff=.8
    AnnotateWithRefSeq(inDir, outDir, gene_path, rpt_path)
    pass

if __name__ == '__main__':
    main(sys.argv)

