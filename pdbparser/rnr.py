
#See COPYING for license 

import numpy as np
import numpy.lib.recfunctions
from itertools import groupby, repeat
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
    changelock=np.where(coord['segid']==seg)
    changes=np.copy(coord[coord['segid']==seg])
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
    coord[changelock]=changes
    return coord

def renumberatom(coord):
    totnr=len(coord)
    for atnr in range(0,totnr):
        coord[atnr]['atnr'] = atnr+1
    return coord
    
def rnrseg_charmm(coord,seg,merge):
    if len(seg)==1:
        renumbered=renumberseg(coord,seg,0)
        return renumbered
    elif len(seg) > 1 and merge is True:
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
    else:
        renumbered=np.copy(coord)
        segid=0
        while segid <len(seg):
            resstart=0
            renumbered=renumberseg(renumbered,seg[segid],resstart)
            segid=segid+1
        return renumbered

def replace(coor,old,new):                             
    for x,y in zip(old,new):                        
        ind=coor['resname'] == x                       
        coor['resname'][ind] = y                       
    return coor                                        

def addbfact(coor,val):
    #($resn =~ /CLA/) is excluded
    reslist=['SOL','HOH','TIP3','SWM4','POT','SOD','CAL','NA','ETOH']
    restorest=coor['resname'].tolist()
    for res in reslist:
        if res in restorest:
            restorest.remove(res)
    for res in restorest:
        location=np.where(coor['resname']==res)
        coor['tfact'][location]= val
    # Now put 0 back for drude and lone pair and hydrogen
    for at in ['H','D','L']:
        print 'Here'
        location=np.where(coor['atname'][0:1]==at)
        print coor['atname']
        coor[location]['tfact']=0.0
        # Didn't work
        #[np.flatnonzero(np.char.find(coor['atname'][0:][0][0:1],at)!=-1)]=0.0
    return coor
