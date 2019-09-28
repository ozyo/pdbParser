import math
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
    chains=list(sorted(chains))
    com_all=np.array(com(ca))
    xyz=np.array(ca[['x','y','z']].tolist())
    xyz=xyz-com_all
    # For ordering https://github.com/pierrepo/principal_axes/blob/master/principal_axes.py
    # and mdanalysis routines merged
    inertia = np.dot(xyz.transpose(), xyz)
    e_values, e_vectors = np.linalg.eig(inertia)
    order = np.argsort(e_values)
    eval3, eval2, eval1 = e_values[order]
    axis3, axis2, axis1 = e_vectors[:, order].transpose()
    rotxyz=np.array([axis1,axis2,axis3])
    for i in range(xyz.shape[0]):
        t=np.dot(rotxyz,xyz[i])
        ca['x'][i]=round(t[0],3)
        ca['y'][i]=round(t[1],3)
        ca['z'][i]=round(t[2],3)
    print("XYZ axis limits:\n   " , xyz[:,0].min(), "to", xyz[:,0].max() )
    print("   " , xyz[:,1].min(), "to", xyz[:,1].max() )
    print("   " , xyz[:,2].min(), "to", xyz[:,2].max())
    #ca['x'][np.where(ca['x'])]=xyz[:,0]
    #ca['y'][np.where(ca['y'])]=xyz[:,1]
    #ca['z'][np.where(ca['z'])]=xyz[:,2]

    coms=[]
    #for ch in range(len(chains)):
    #    if ch < 2:
    #        coms.append(com(ca,chain=chains[ch]))
    coms.append(com(ca))
    coms.append(com(ca,chain=chains[0]))
    coms.append(com(ca,chain=chains[1]))
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
        print(ids)
        directions[ids]=direction(uni,chorderlist[ids])
    print(directions)
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
        #We could simply loop over neworders improve this part.
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
    return list(neworders.keys())



