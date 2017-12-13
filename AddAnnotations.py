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
from ReferenceAnalyzer.ReferenceRepeatAnalyzer import *


class TandemLine():
    def __init__(self, row):
        self.chrom=row[0]
        self.rpt_length=int(row[1])
        self.copynumber=int(row[2])
        self.info=float(row[3])
        self.identity=float(row[4])
        try:self.repetition=float(row[5])
        except: self.repetition='err'
        self.start=int(row[6])
        self.end=int(row[7])
        self.starts=[int(s) for s in row[8].split(',')[:-1] ]
        self.ends=[int(e) for e in row[9].split(',')[:-1] ]
        self.info_string=row[10]
        self.seq=row[11]



class RefseqLine():
    def __init__(self, row):
        pass
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
        self.name=row[12]



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
    gen_array=LoadRefseqAnnotations(annot_file)
    rpt_array=LoadRepeatAnnotations(rpt_file)
    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter='\t')
    outhandle=open(outfile, 'w')
    outtable=csv.writer(outhandle, delimiter='\t')
    header=intable.next()
    outtable.writerow(header+['Features'])
    for row in intable:
        line=TandemLine(row)
        tandem_iv=HTSeq.GenomicInterval(line.chrom, line.start, line.end)
        annotations=','.join( list(set.union(* [ s[1] for s in gen_array[tandem_iv].steps()])))
        rpt_annotations=','.join( list(set.union(* [ s[1] for s in rpt_array[tandem_iv].steps()])))
        outtable.writerow(row+[annotations, rpt_annotations])
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


    seq=GetSeq(ref_file)
    for row in intable:
        chrom, interval=row[-1].split(':')
        left,right=interval.split('-')
        left=int(''.join(left.split(',')))
        right=int(''.join(right.split(',')))
        try:
            coord_left=lo.convert_coordinate(chrom,left)[0][1]
        except:
            print row
            print left
            print lo.convert_coordinate(chrom,left)
        chromosome,coord_right=lo.convert_coordinate(chrom,right)[0][:2]
        new_line=row[:-1]+['{0}:{1}-{2}'.format(chromosome, coord_left, coord_right)]
        seq_name='>{0}_{1}_{2}_{3}_{4}_{5}\n',format(row[7],row[2],row[3], chromosome, coord_left, coord_right,)
        outfile.write(seq_name)
        outfile.write ( '{0}\n'.format( seq[chromosome][coord_left:coord_right].upper() ) )
        updated_table.writerow(new_line)
    inhandle.close()
    updated_handle.close()
    outfile.close()


def main():
    pass

if __name__ == '__main__':
    main()
