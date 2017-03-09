
#See COPYING for license 

import numpy as np
import logging
from writepdb import writecharmm
import numpy.lib.recfunctions
import os
from readpdb import readall
from sep_seg import remove_field_name
def read_charmm(atomlines):
    atomlines=readall(atomlines)
    coords=[]
    for atom in atomlines:
        atnr=int(atom[6:11].strip())
        atname=str(atom[12:16].strip())
        altloc=str(atom[16].strip())
        resname=str(atom[17:21].strip())
        ch=str(atom[21])
        resnr=int(atom[22:26])
        icode=str(atom[26].strip())
        x=float(atom[30:38].strip())
        y=float(atom[38:46].strip())
        z=float(atom[46:54].strip())
        try:
            occu=float(atom[54:60].strip())
        except ValueError:
            occu=0.0
        try:
            tfact=float(atom[60:66].strip())
        except ValueError:
            tfact=0.0
        segid=str(atom[72:76].strip())
        coords.append((atnr,atname,altloc,resname,ch,resnr,icode,x,y,z,occu,tfact,segid))
    coords=np.array(coords,dtype=('i,S4,S4,S4,S4,i,S4,f,f,f,f,f,S4'))
    coords.dtype.names=('atnr','atname','altloc','resname','ch','resnr','icode','x','y','z','occu','tfact','segid')
    return coords

def renumber(coord,seg):
    totnr=len(coord[coord['segid']==seg])/1 #3 for full atomistic
    print totnr
    i=np.where(coord['segid']==seg)
    for res in range(0,totnr):
        coord[i[0][1*res]]['resnr'] = res+1 # 3* res for full atomistic
        #coord[i[0][3*res+1]]['resnr'] = res+1 # same above
        #coord[i[0][3*res+2]]['resnr'] = res+1 #same above
    return coord

def renumberatom(coord):
    totnr=len(coord)
    for atnr in range(0,totnr):
        print coord[atnr]['atnr']
        coord[atnr]['atnr'] = atnr+1
    return coord
    
def rnrseg_charmm(coord,seg,cwd):
    atoms=read_charmm(coord)
    #renumbered=remove_field_name(renumber(atoms,seg),'icode')
    renumbered=renumber(atoms,seg)
    writecharmm(renumbered,cwd+'/renumbered.pdb')
    exit

def rnratnr_charmm(coord,cwd):
    atoms=read_charmm(coord)
    #renumbered=remove_field_name(renumberatom(atoms),'icode')
    renumbered=renumberatom(atoms)
    writecharmm(renumbered,cwd+'/renumbered.pdb')
    exit
