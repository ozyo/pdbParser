#Things to do
#obtain files from previous processes and create a sequence alignment between CA atoms
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SeqUtils import seq1 as letter
import numpy as np

#tags: a aligned, s sequence, map residue name-nr map

matrix = matlist.blosum62

def getseq(ca):
    seq=letter("".join(ca['resname'].tolist()))
    resmap=zip(list(seq),ca['resnr'].tolist())
    return seq, resmap

def align(seq1,seq2):
    if isinstance(seq1, basestring):
        alignments=pairwise2.align.globaldd(seq1,seq2,matrix,-5,-5,-5,-.5)
        return alignments
    elif isinstance(seq1, list):
        alignments=pairwise2.align.globaldd(seq1,seq2,matrix,-5,-5,-5,-.5,gap_char=['-'])
        return alignments

def findgap(aca):
    global start
    global end
    if aca[0] == '-' or aca[-1] == '-':
        start=next(let for let,v in enumerate(aca) if v != '-')
        end=(aca[start:-1].index('-')+start)-len(aca)-1
    else:
        start=None
        end=None
    return start,end

def getaligned(ca1,ca2):
    sca1,mapca1=getseq(ca1)
    sca2,mapca2=getseq(ca2)
    #Grabs the first alignment
    aligned=align(sca1,sca2)[0]
    aca1=list(aligned[0])
    aca2=list(aligned[1])
    shift1=findgap(aca1)
    print shift1
    shift2=findgap(aca2)
    if start or end is None:
        pass
    else:
        for i in range(0,shift1[0]):
            mapca1.insert(i,('-','-'))

