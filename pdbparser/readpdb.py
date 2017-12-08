
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
def checkmulti(pdb):
    for line in pdb:
        if 'NMR' in line:
            logging.critical('Current version cannot handle multipdb files.')
            logging.critical('Please upload your own structures.')
            exit()
        else:
            pass
def readcompnd(pdb):
    compnd=None
    compndlist=[]
    molin=[]
    read=False
    for line in pdb:
        if 'COMPND' in line:
            compndlist.append(line[10:-1].strip())
    for line in compndlist:
        if line.startswith('MOL_ID:'):
            molin.append(compndlist.index(line))
    if len(molin) >2:
        logging.critical('Current version cannot handle this pdb file.')
        logging.critical('Please provide your input files')
    elif len(molin) == 2:
        molid=1
        molid1=compndlist[molin[0]:molin[1]] 
        molid2=compndlist[molin[1]:-1]
        checks=['DNA','RNA','5\'','3\'']
        for line in molid1:
            if any(s in line for s in checks):
                molid=2
        if molid == 1:
            for line in molid1:
                if 'CHAIN:' in line:
                    compnd=[i.strip().strip(';') for i in line.split(':')[1].split(',')]
        else:
            for line in molid2:
                if 'CHAIN:' in line:
                    compnd=[i.strip().strip(';') for i in line.split(':')[1].split(',')]
    else:
        for line in compndlist:
            if 'CHAIN:' in line:
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

def readhetatm(pdb):
    hetatoms=[]
    read=False
    for line in pdb:
        if line.startswith('HETATM'):
            hetatoms.append(line)
    return hetatoms

def coord(atomlines,coortype):
    # This will read the coordinates with no icode
    # Instead it will use the filed to write residue numbers above 9999
    # If icode is required read it with decreasing resnr to 26 and adding
    # an icode field
    coords=[]
    for atom in atomlines:
        atnr=str(atom[6:11].strip())
        atname=str(atom[12:16].strip())
        altloc=str(atom[16].strip())
        resname=str(atom[17:21].strip())
        ch=str(atom[21])
        resnr=str(atom[22:27])
        x=str(atom[30:38].strip())
        y=str(atom[38:46].strip())
        z=str(atom[46:54].strip())
        try:
            occu=str(atom[54:60].strip())
        except ValueError:
            occu=0.0
        try:
            tfact=str(atom[60:66].strip())
        except ValueError:
            tfact=0.0
        if coortype == 'charmm':
            segid=str(atom[71:77].strip())
            coords.append((atnr,atname,altloc,resname,ch,resnr,x,y,z,occu,tfact,segid))
        else:
            element=str(atom[76:78].strip())
            charge=str(atom[78:80].strip())
            coords.append((atnr,atname,altloc,resname,ch,resnr,x,y,z,occu,tfact,element,charge))
    if coortype == 'charmm':
        coords=np.array(coords,dtype=('S5,S4,S1,S4,S1,S5,S8,S8,S8,S6,S6,S5'))
        coords.dtype.names=('atnr','atname','altloc','resname','ch','resnr','x','y','z','occu','tfact','segid')
    else:
        coords=np.array(coords,dtype=('S5,S4,S1,S4,S1,S5,S8,S8,S8,S6,S6,S2,S2'))
        coords.dtype.names=('atnr','atname','altloc','resname','ch','resnr','x','y','z','occu','tfact','element','charge')
    return coords
