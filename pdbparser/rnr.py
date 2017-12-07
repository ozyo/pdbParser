
#See COPYING for license 

import numpy as np
import logging
import numpy.lib.recfunctions
import os
from itertools import groupby
from operator import itemgetter

def check_resnr_dup(location,mode):
    if mode == 2:
        # This only ensures that the renumbered residues will not overlap with the remainings
        # If there are multiple occurences from the beginning, may the force help you. 
        # I don't have a check function for that yet, implementing
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
    elif mode == 1:
        loc_list=[]
        for k, g in groupby(enumerate(location), lambda (i, x): i-x):
            loc_list.append(map(itemgetter(1), g))
            return loc_list

def renumberseg(coord,seg,resstart):
    #Renumbering is only done based on segid here
    changes=np.copy(coord)
    #Group by is a better way in case there residues with the same number.
    #curres=np.unique(coord[coord['segid']==seg]['resnr'])
    curres=coord[coord['segid']==seg]['resnr']
    unires=[k for k,g in groupby(curres) if k!=0]
    alreadypassed=[]
    for res in range(0,len(unires)):
        resnr=unires[res]
        new=resstart+1+res
        location=np.where((changes['segid']==seg)&(changes['resnr']==resnr))
        alloc=check_resnr_dup(location[0],1)
        if alloc[0] in alreadypassed:
            checkedloc=alloc[1]
        else:
            alreadypassed.append(alloc[0])
            checkedloc=alloc[0]
            #checkedloc=check_resnr_dup(location[0],2)
        changes['resnr'][checkedloc]=new
    return changes

def renumberatom(coord):
    totnr=len(coord)
    for atnr in range(0,totnr):
        coord[atnr]['atnr'] = atnr+1
    return coord
    
def rnrseg_charmm(coord,seg):
    if len(seg)==1:
        renumbered=renumberseg(coord,seg,0)
        return renumbered
    else:
        renumbered=np.copy(coord)
        segid=0
        while segid <len(seg):
            if segid==0:
                resstart=0
            else:
                resstart=renumbered[renumbered['segid']==seg[segid-1]]['resnr'][-1]
            renumbered=renumberseg(renumbered,seg[segid],resstart)
            segid=segid+1
        return renumbered

def replace(coor,old,new):                             
    for x,y in zip(old,new):                        
        ind=coor['resname'] == x                       
        coor['resname'][ind] = y                       
    return coor                                        
