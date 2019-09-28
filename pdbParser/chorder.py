import numpy as np
from copy import deepcopy
from pdbParser import parser as pP
from pdbParser import writepdb as wp

def com(ca,chain=None):
    #This is not really a center of mass calculation, no mass weight is considered
    if chain:
        ca=ca[ca['ch']==chain]
    com=[np.mean(np.array(ca[['x']],dtype="float64")),np.mean(np.array(ca[['y']],dtype="float64")),np.mean(np.array(ca[['z']],dtype="float64"))]
    return(com)

def direction(ca,chains):
    coms=[]
    for ch in range(len(chains)):
        if ch < 2:
            coms.append(com(ca,chain=chains[ch]))
    coms.append(com(ca))
    coms=np.array(coms)
    if np.linalg.det(coms) > 0:
        #If matrix determinant is positive, anticlockwise. Else clockwise
        return "ac"
    else:
        return "c"

def change_direction(chlist):
    rev=list(reversed(chlist))
    fixed=chlist[0:1]+rev[:-1]
    return(fixed)

def direction_check(uni_list,chorderlist):
    directions={}
    for ids,uni in uni_list.items():
        directions[ids]=direction(uni,chorderlist[ids])
    acs=list(directions.values()).count('ac')
    cs=list(directions.values()).count('c')
    reorder=[]
    if acs > cs:
        for ids in directions:
            if directions[ids] == "c":
                reorder.append(ids)
    else:
        for ids in directions:
            if directions[ids] == "ac":
                reorder.append(ids)
    neworders={}
    for ids in reorder:
        neworders[ids]=change_direction(chorderlist[ids])
    return(neworders)

def sepchains(ca,chlist):
    sep={}
    for ch in chlist:
        sep[ch]=deepcopy(ca[ca['ch']==ch])
    return sep

def reorder(id_ch_dict,mer,cwd,altloc='A',tag=None):
    uni_list={}
    for ids in id_ch_dict.keys():
        if tag:
            name=tag+ids
        else:
            name=ids
        uni_list[ids]=pP.parse_ca(open(cwd+'/'+name,'r'),'NA',mer,altloc,list(sorted(id_ch_dict[ids])))
    neworders=direction_check(uni_list,id_ch_dict)
    for ids in uni_list:
        if ids in neworders.keys():
            newchains=[]
            sep_chains=sepchains(uni_list[ids],id_ch_dict[ids])
            oldchs=id_ch_dict[ids]
            for ind,ch in enumerate(neworders[ids]):
                oldch=oldchs[ind]
                tmp=deepcopy(sep_chains[ch])
                changelock=np.where(tmp['ch']==ch)
                tmp['ch'][changelock]=oldch
                newchains.append(tmp)
            merged=np.concatenate(newchains)
            wp.writeca(merged,cwd+'/'+'reordered_'+ids)




