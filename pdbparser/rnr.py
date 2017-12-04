
#See COPYING for license 

import numpy as np
import logging
from writepdb import writecharmm_noicode
import numpy.lib.recfunctions
import os
from readpdb import readall
from sep_seg import remove_field_name
from itertools import groupby
from operator import itemgetter

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

def check_resnr_dup(location):
    loc_list=[]
    for k, g in groupby(enumerate(location), lambda (i, x): i-x):
        loc_list.append(map(itemgetter(1), g))
    if len(loc_list) == 1:
        return loc_list[0]
    elif len(loc_list) == 2:
        return loc_list[1]
    else:
        'Hey I detected multiple locations, I cannot renumber these set of residues'
        exit

def renumber(coord,seg,resstart):
    changes=np.copy(coord)
    #curres=np.unique(coord[coord['segid']==seg]['resnr'])
    curres=coord[coord['segid']==seg]['resnr']
    unires=[k for k,g in groupby(curres) if k!=0]
    for res in range(0,len(unires)):
        resnr=unires[res]
        new=resstart+1+res
        location=np.where((changes['segid']==seg)&(changes['resnr']==resnr))
        checkedloc=check_resnr_dup(location[0])
        changes['resnr'][checkedloc]=new
    return changes

def renumberatom(coord):
    totnr=len(coord)
    for atnr in range(0,totnr):
        coord[atnr]['atnr'] = atnr+1
    return coord
    
def rnrseg_charmm(coord,seg,cwd):
    atoms=read_charmm(coord)
    if len(seg)==1:
        print seg
        #renumbered=remove_field_name(renumber(atoms,seg),'icode')
        renumbered=renumber(atoms,seg,0)
        writecharmm_noicode(renumbered,cwd+'/renumbered.pdb')
        exit
    else:
        renumbered=np.copy(atoms)
        segid=0
        while segid <len(seg):
            if segid==0:
                resstart=0
            else:
                resstart=renumbered[renumbered['segid']==seg[segid-1]]['resnr'][-1]
            renumbered=renumber(renumbered,seg[segid],resstart)
            segid=segid+1
        writecharmm_noicode(renumbered,cwd+'/renumbered.pdb')
        exit

def rnratnr_charmm(coord,cwd):
    atoms=read_charmm(coord)
    #renumbered=remove_field_name(renumberatom(atoms),'icode')
    renumbered=renumberatom(atoms)
    writecharmm_noicode(renumbered,cwd+'/renumbered.pdb')
    exit
