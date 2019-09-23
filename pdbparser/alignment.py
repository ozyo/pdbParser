#See COPYING for license 
#Things to do
#obtain files from previous processes and create a sequence alignment between CA atoms
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SeqUtils import seq1 as letter
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
import numpy as np
import logging,os

#tags: a aligned, s sequence, map residue name-nr map

matrix = matlist.blosum62

def getseq(ca):
    seq=letter("".join(ca['resname'].tolist()))
    resmap=zip(list(seq),ca['resname'].tolist(),ca['resnr'].tolist())
    return seq, resmap

def align(seq1,seq2):
    if isinstance(seq1, basestring):
        alignments=pairwise2.align.globalds(seq1,seq2,matrix,-11,-11,penalize_end_gaps=True,)
        return alignments
    elif isinstance(seq1, list):
        alignments=pairwise2.align.globalds(seq1,seq2,matrix,-11,-11,gap_char=['-'],penalize_end_gaps=True)
        return alignments

def findgap(aca):
    #This function finds the first patch without gap. For an MSA we need something that removes the structure if it is broken.
    if aca[0] == '-' and aca[-1] != '-':
        start=next(ind for ind,gap in enumerate(aca) if gap != '-')
        end=len(aca)+1
    elif aca[-1] == '-' and aca[0] != '-':
        end=aca.index('-')-len(aca)#-1
        start=0
    elif aca[0] == '-' and aca[-1] == '-':
        start=next(ind for ind,gap in enumerate(aca) if gap != '-')
        try:
            end=aca[start:-1].index('-')+start-1
        except ValueError: #If only the last is - and the ones before are a sequence
            end=aca[start:].index('-')+start-1
    else:
        start=0
        end=len(aca)+1
    return start,end

def findinsertions(aca):
    #Use this function with a reference sequence alignment result to find insertions
    for i in range(0,len(aca)):
        if aca[i] != '-':
            start=i
            break
    subaca=aca[i:]
    lastelem=len(subaca) - 1
    while lastelem>=0:
        if subaca[lastelem]!='-':
            end=lastelem+start
            break
        lastelem=lastelem-1
    insertions=[]
    for i in range(start,end+1):
        if aca[i] == '-':
            insertions.append(i)
    return insertions
        

def findtergap(aca,insertions):
    for i in range(0,len(aca)):
        if aca[i] != '-':
            start=i
            break
    subaca=aca[i:]
    lastelem=len(subaca) - 1
    while lastelem>=0:
        if subaca[lastelem]!='-':
            end=lastelem+start
            break
        lastelem=lastelem-1
    broken=False
    insertion=False
    for i in range(start,end+1):
        if i in insertions and aca[i]=='-':
            continue
        elif i in insertions and aca[i]!='-':
            insertion=True
            break
        else:
            if aca[i]=='-':
                broken=True
                break
    if broken is True:
        return None,None
    if insertion is True:
        return None,None
    else:
        return start,end

def getaligned(ca1,ca2):
    sca1,mapca1=getseq(ca1)
    sca2,mapca2=getseq(ca2)
    #Grabs the first alignment
    aligned=align(sca1,sca2)[0]
    print(aligned)
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
    core1=ca1[(ca1['resnr']>=nter1[2]) & (ca1['resnr']<=cter1[2])]
    core2=ca2[(ca2['resnr']>=nter2[2]) & (ca2['resnr']<=cter2[2])]
    if len(core1) == len(core2):
        logging.info('Run successful proceeding with eBDIMS calculation')
        logging.critical('SUCCESS')
        return core1, core2, True
    else:
        logging.critical('Different number of atoms.')
        logging.critical('I am not extracting the same region for these structures')
        logging.critical('Please upload your own structures to continue. ')
        logging.critical('FAIL')
        return core1, core2, False
        exit()
        

def multialigned(ca1,ca2,mer):
    whole1=None
    whole2=None
    correct=[]
    if len(np.unique(ca1['ch'])) == len(np.unique(ca2['ch'])) and len(np.unique(ca1['ch'])) == mer:        
        for a,b in zip(np.unique(ca1['ch']),np.unique(ca2['ch'])):
            cores=getaligned(ca1[ca1['ch']==a],ca2[ca2['ch']==b])
            if whole1 is None:
                whole1=cores[0]
                whole2=cores[1]
                correct.append(cores[2])
            elif whole1 is not None:
                whole1=np.append(whole1,cores[0],axis=0)
                whole2=np.append(whole2,cores[1],axis=0)
                correct.append(cores[2])
    if False in correct:
        logging.critical('Different number of atoms.')
        logging.critical('I am not extracting the same region for these structures')
        logging.critical('Please upload your own structures to continue. ')
        logging.critical('FAIL')        
        return whole1, whole2, False
    if False not in correct:
        return whole1, whole2, True


def findresid(shifts,nter,cter,resmap):
    resmap=open(resmap,'r').read().split('>')
    resid={}
    for pdb in resmap:
        if len(pdb)>2:
            id,nr=pdb.split('\n')[0:2]
            if id in shifts.keys():
                if len(shifts[id]) == 0:
                    resid[id]=[None,None]
                    continue
                noffset=nter-shifts[id][0]
                coffset=cter-nter+noffset
                nrs=nr.split('-')[noffset:coffset+1]
                try:
                    resid[id]=[int(nrs[0]),int(nrs[-1])]
                except IndexError:
                    logging.critical('PDB ID %s contains too little sequence or the alignment is problematic. Remove this structure and try again.' %id)
                    exit()
    return(resid)

def msa_clustal(infile,resmap,outfile,clustalopath,cwd,merinfo,query,totmer,alnf=None):
    resmap=cwd+'/'+resmap
    if alnf is None:
        clustalomega_cline = ClustalOmegaCommandline(infile=infile, outfile=outfile, verbose=False, auto=True, force=True)
        clustalomega_cline()
        msa=AlignIO.read(outfile, "fasta")
    else:
        msa=AlignIO.read(cwd+'/'+alnf, "fasta")
    broken=[]
    nter=0
    cter=msa.get_alignment_length()
    shifts={}
    for record in msa:
        if record.id == 'refseq':
            insertions=findinsertions(str(record.seq))
            break
    for record in msa:
        if record.id == 'refseq':
            continue
        aln=str(record.seq)
        start,end=findtergap(aln,insertions)
        if start is None and end is None:
            broken.append(record.id.split('|')[0])
            shifts[record.id]=[]
            logging.critical('%s contains insertions or missing residues, skipping this chain and the assembly it belongs to.' %record.id)
            continue
        else:
            shifts[record.id]=[start,end]
            if start > nter:
                nter=start
            if end < cter:
                cter=end
    core=msa[:,nter:cter+1]
    resid=findresid(shifts,nter,cter,resmap)
    completemers={}
    fullids=list(set([key.split('|')[0] for key in resid.keys()]))
    for pdb in fullids:
        pdbid,mer=pdb.split('.')[0].split('_')
        amer=[]
        merinfo[pdbid]
        for ch in merinfo[pdbid][1][int(mer)-1]:
            if pdb+'|'+ch+'|' in broken:
                continue
            else:
                amer.append(ch)
        completemers[pdb]=amer
    AlignIO.write(core,cwd+'/'+query+'_core.fasta',format='fasta')
    return(completemers,resid,list(set(broken)))
