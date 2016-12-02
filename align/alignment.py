#See COPYING for license 
#Things to do
#obtain files from previous processes and create a sequence alignment between CA atoms
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SeqUtils import seq1 as letter
import numpy as np
import logging

#tags: a aligned, s sequence, map residue name-nr map

matrix = matlist.blosum62

def getseq(ca):
    seq=letter("".join(ca['resname'].tolist()))
    resmap=zip(list(seq),ca['resname'].tolist(),ca['resnr'].tolist())
    return seq, resmap

def align(seq1,seq2):
    if isinstance(seq1, basestring):
        alignments=pairwise2.align.globalms(seq1,seq2,2,0,-5,-5)
        return alignments
    elif isinstance(seq1, list):
        alignments=pairwise2.align.globalms(seq1,seq2,2,0,-5,-5,gap_char=['-'])
        return alignments

def findgap(aca):
    global start
    global end
    if aca[0] == '-' and aca[-1] != '-':
        start=next(ind for ind,gap in enumerate(aca) if gap != '-')
        end=-1
    elif aca[-1] == '-' and aca[0] != '-':
        end=aca.index('-')-len(aca)-1
        start=0
    elif aca[0] == '-' and aca[-1] == '-':
        start=next(ind for ind,gap in enumerate(aca) if gap != '-')
        end=aca[start:-1].index('-')+start
    else:
        start=0
        end=-1
    return start,end

def getaligned(ca1,ca2):
    sca1,mapca1=getseq(ca1)
    sca2,mapca2=getseq(ca2)
    #Grabs the first alignment
    aligned=align(sca1,sca2)[0]
    aca1=list(aligned[0])
    aca2=list(aligned[1])
    shift1=findgap(aca1)
    shift2=findgap(aca2)
    #map the gaps to the opposite structure, this is important!!!
    nter1=mapca1[shift2[0]:shift2[1]][0]
    cter1=mapca1[shift2[0]:shift2[1]][-1]
    nter2=mapca2[shift1[0]:shift1[1]][0]
    cter2=mapca2[shift1[0]:shift1[1]][-1]
    #since we mapped it to the opposite in the steps above, we can return back to normal
    core1=ca1[(ca1['resnr']>nter1[2]) & (ca1['resnr']<cter1[2])]
    core2=ca2[(ca2['resnr']>nter2[2]) & (ca2['resnr']<cter2[2])]
    if len(core1) == len(core2):
        return core1, core2
    else:
        print len(core1), len(core2)
        logging.critical('This was not supposed to happen')
        logging.critical('I am not extracting the same region for both sequences, returning the structures for debugging.')
        return core1, core2

def multialigned(ca1,ca2,mer):
    chains1=np.array_split(ca1,mer)
    chains2=np.array_split(ca2,mer)
    whole1=None
    whole2=None
    for a,b in zip(chains1,chains2):
        cores=getaligned(a,b)
        if whole1 is None:
            whole1=cores[0]
            whole2=cores[1]
        elif whole1 is not None:
            whole1=np.append(whole1,cores[0],axis=0)
            whole2=np.append(whole2,cores[1],axis=0)

    return whole1, whole2
