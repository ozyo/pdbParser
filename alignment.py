#Things to do
#obtain files from previous processes and create a sequence alignment between CA atoms
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

matrix = matlist.blosum62

def align(seq1,seq2):
    if isinstance(seq1, basestring):
        alignments=pairwise2.align.globaldd(seq1,seq2,matrix,-5,-4,-3,-.1)
        return alignments
    elif isinstance(seq1, list):
        alignments=pairwise2.align.globaldd(seq1,seq2,matrix,-5,-4,-3,-.1,gap_char=['-'])
        return alignments
