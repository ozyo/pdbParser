#We need to read quite a bit of information from the REMARKS lines. So it is better to keep them in a seperate module.
import numpy as np
from urllib2 import urlopen
import logging

def getpdb(pdbid):
    url='http://files.rcsb.org/view/%s.pdb' %pdbid
    return urlopen(url).readlines()

def readcompnd(pdb):
    compnd=None
    read=False
    for line in pdb:
        if 'COMPND' in line:
            if 'MOL_ID: 1' in line:
                read=True
            elif 'MOL_ID: 2' in line:
                read=False
            if read is True and 'CHAIN:' in line:
                compnd=[i.strip().strip(';') for i in line.split(':')[1].split(',')]
    return compnd

def readremark(pdb,compnd):
    remark=[]
    read=False
    for line in pdb:
        if 'REMARK 465' in line:
            remark.append(line)
    r465=np.genfromtxt(remark[7:-1],names=['REMARK','465','rname','ch','rid'],dtype=['S6',int,'S3','S1',int])
    filt=r465[np.in1d(r465['ch'],compnd)]
    return filt

def readseq(pdb,compnd):
    seqres={}
    for ch in compnd:
        seqres[ch]=[]
    read=False
    for line in pdb:
        if 'SEQRES' in line:
            ch=line[11]
            if ch in compnd:
                seq=line[19:-1].split()
                seqres[ch]=seqres[ch]+seq
    return seqres

def readatom(pdb):
    atoms=[]
    read=False
    for line in pdb:
        if line.startswith('ATOM'):
            atoms.append(line)
    return atoms

def coord(atomlines,compnd):
    coords=[]
    for atom in atomlines:
        atnr=int(atom[6:11].strip())
        atname=str(atom[12:16].strip())
        altloc=str(atom[16].strip())
        resname=str(atom[17:20].strip())
        ch=str(atom[21])
        resnr=int(atom[22:26])
        icode=str(atom[26].strip())
        x=float(atom[30:38].strip())
        y=float(atom[38:46].strip())
        z=float(atom[46:54].strip())
        occu=float(atom[54:60].strip())
        tfact=float(atom[60:66].strip())
        element=str(atom[76:78].strip())
        charge=str(atom[78:80].strip())
        coords.append((atnr,atname,altloc,resname,ch,resnr,icode,x,y,z,occu,tfact,element,charge))
    coords=np.array(coords,dtype=('i,S4,S4,S4,S4,i,S4,f,f,f,f,f,S4,S4'))
    coords.dtype.names=('atnr','atname','altloc','resname','ch','resnr','icode','x','y','z','occu','tfact','element','charge')
    filt=coords[np.in1d(coords['ch'],compnd)]
    return filt

#pdb=open('3RIF.pdb').readlines()
#read=pdbread()
#filt=read.readRemark(pdb)
#print filt
