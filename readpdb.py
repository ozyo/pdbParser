#We need to read quite a bit of information from the REMARKS lines. So it is better to keep them in a seperate module.
import numpy as np
from urllib2 import urlopen
import logging

def getpdb(pdbid):
    url='http://files.rcsb.org/view/%s.pdb' %pdbid
    return urlopen(url).readlines()

def readcompnd(pdb,pdbid,mer):
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
    print 'Detected chains for %s are ' % (pdbid)+' '.join(i for i in compnd) 
    if len(compnd) is not 1:
        print 'If multimeric option not chosen continuing with the chain %s' %(compnd[0])
    if mer is not False:
        logging.info('Multimeric option chosen...')
        if mer < len(compnd):
            logging.warning('%s structure contains more biological assemblies than one') % (pdbid)
            if mer % len(compnd) != 0:
                warnings.warn('Cannot determine the biological assembly, please provide your own input files')
                exit()
            else:
                logging.warning('Attempting to seperate them')
        elif mer == 1 and len(compnd) != 1:
            logging.warning('Assuming you know what you are doing.')
            logging.warning('You chose only 1 chain for the analysis but there are more chains in the biological assembly')
        elif mer > len(compnd):
            logging.warning("Structure contains less chains then the supplied multimeric option. Assuming not complete structure, please provide your own input files")
            exit()
    return compnd

def readmissing(pdb):
    compnd=readcompnd(pdb)
    remark=[]
    read=False
    for line in pdb:
        if 'REMARK 465' in line:
            remark.append(line)
    r465=np.genfromtxt(remark[7:-1],names=['REMARK','465','rname','ch','rid'],dtype=['S6',int,'S3','S1',int])
    filt=r465[np.in1d(r465['ch'],compnd)]
    return filt

def readatom(pdb):
    atoms=[]
    read=False
    for line in pdb:
        if line.startswith('ATOM'):
            atoms.append(line)
    return atoms

def coord (pdb,pdbid,mer):
    coords=[]
    atoms=readatom(pdb)
    compnd=readcompnd(pdb,pdbid,mer)
    for atom in atoms:
        atnr=int(atom[6:11].strip())
        atname=atom[12:16].strip()
        altloc=atom[16]
        resname=atom[17:20].strip()
        ch=atom[21]
        resnr=int(atom[22:26])
        icode=atom[26]
        x=float(atom[30:38].strip())
        y=float(atom[38:46].strip())
        z=float(atom[46:54].strip())
        occu=float(atom[54:60].strip())
        tfact=float(atom[60:66].strip())
        element=atom[76:78]
        charge=atom[78:80]
        coords.append((atnr,atname,altloc,resname,ch,resnr,icode,x,y,z,occu,tfact,element,charge))
    coords=np.array(coords,dtype=('int,str,str,str,str,int,str,float,float,float,float,float,str,str'))
    coords.dtype.names=('atnr','atname','altloc','resname','ch','resnr','icode','x','y','z','occu','tfact','element','charge')
    filt=coords[np.in1d(coords['ch'],compnd)]
    return filt

#pdb=open('3RIF.pdb').readlines()
#read=pdbread()
#filt=read.readRemark(pdb)
#print filt
