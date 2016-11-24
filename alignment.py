#Things to do
#obtain files from previous processes and create a sequence alignment between CA atoms
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SeqUtils import seq1 as letter
import numpy as np

matrix = matlist.blosum62

def getseq(ca):
    seq=letter("".join(ca['resname'].tolist()))
    resmap=zip(list(seq),ca['resnr'].tolist())
    return seq, resmap

def align(seq1,seq2):
    if isinstance(seq1, basestring):
        alignments=pairwise2.align.globaldd(seq1,seq2,matrix,-5,-4,-3,-.1)
        return alignments
    elif isinstance(seq1, list):
        alignments=pairwise2.align.globaldd(seq1,seq2,matrix,-5,-4,-3,-.1,gap_char=['-'])
        return alignments
