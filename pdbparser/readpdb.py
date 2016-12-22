#See COPYING for license 

#We need to read quite a bit of information from the REMARKS lines. So it is better to keep them in a seperate module.
import numpy as np
import urllib2 
import logging

def getpdb(pdbid):
    try:
        url='http://files.rcsb.org/view/%s.pdb' %pdbid
        return urllib2.urlopen(url).readlines()
    except urllib2.HTTPError as err:
        if err.code == 404:
            logging.critical('PDB code not found')
        else:
            logging.critical('Cannot reach PDB database, unknown error')
        exit()

def readcompnd(pdb):
    compnd=None
    compndlist=[]
    read=False
    for line in pdb:
        if 'COMPND' in line:
            compndlist.append(line)

    for line in compndlist:
        try:
            if 'MOLECULE:' and '5\'' in line:
                compndlist.remove('COMPND MOL_ID: 1')
        except ValueError:
            pass

    for line in compndlist:
        if 'MOL_ID:' in line:
            read=True
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
    if compnd is None:
        compnd=[np.unique(coords['ch']).tolist()]
    filt=coords[np.in1d(coords['ch'],compnd)]
    return filt
