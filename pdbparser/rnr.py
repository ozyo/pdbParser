
#See COPYING for license 

import numpy as np
import logging
from writepdb import writecharmm_noicode
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
        resnr=int(atom[22:27])
        #icode=str(atom[27].strip())
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
        coords.append((atnr,atname,altloc,resname,ch,resnr,x,y,z,occu,tfact,segid))
    coords=np.array(coords,dtype=('i,S4,S4,S4,S4,i,f,f,f,f,f,S4'))
    coords.dtype.names=('atnr','atname','altloc','resname','ch','resnr','x','y','z','occu','tfact','segid')
    return coords

def renumber(coord,seg):
    changes=np.copy(coord)
    curres=np.unique(coord[coord['segid']==seg]['resnr'])
    for res in range(0,len(curres)):
        resnr=curres[res]
        new=res+1
        location=np.where((changes['segid']==seg)&(changes['resnr']==resnr))
        #changes['resnr'][location[0]]=new
        #print location
        changes['resnr'][location[0]]=new
    return changes

def renumberatom(coord):
    totnr=len(coord)
    for atnr in range(0,totnr):
        coord[atnr]['atnr'] = atnr+1
    return coord
    
def rnrseg_charmm(coord,seg,cwd):
    atoms=read_charmm(coord)
    #renumbered=remove_field_name(renumber(atoms,seg),'icode')
    renumbered=renumber(atoms,seg)
    writecharmm_noicode(renumbered,cwd+'/renumbered.pdb')
    exit

def rnratnr_charmm(coord,cwd):
    atoms=read_charmm(coord)
    #renumbered=remove_field_name(renumberatom(atoms),'icode')
    renumbered=renumberatom(atoms)
    writecharmm_noicode(renumbered,cwd+'/renumbered.pdb')
    exit
