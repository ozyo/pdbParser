import math
import numpy as np
from copy import deepcopy
from pdbParser import parser as pP
from pdbParser import writepdb as wp
from scipy.spatial.transform import Rotation as R

def com(xyz,weight,tweight):
    coms=np.sum(xyz,axis=0)
    coms=(coms*weight)/tweight
    print(coms)
    return(coms)
def inertia(xyz,com,weight):
    xyz=xyz-com
    Ixx=np.sum(weight*(xyz[:,1]*xyz[:,1]*xyz[:,2]*xyz[:,2]))
    Iyy=np.sum(weight*(xyz[:,0]*xyz[:,0]*xyz[:,2]*xyz[:,2]))
    Izz=np.sum(weight*(xyz[:,1]*xyz[:,1]*xyz[:,2]*xyz[:,2]))
    Ixy=np.sum(-(weight*(xyz[:,0]*xyz[:,1])))
    Ixz=np.sum(-(weight*(xyz[:,0]*xyz[:,2])))
    Iyz=np.sum(-(weight*(xyz[:,1]*xyz[:,2])))
    return(np.array([[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]]))

def get_rotation_matrix(i_v, unit=None):
    #https://stackoverflow.com/questions/43507491/imprecision-with-rotation-matrix-to-align-a-vector-to-an-axis
    # From http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q38
    if unit is None:
        unit = [1.0, 0.0, 0.0]
    # Normalize vector length
    i_v /= np.linalg.norm(i_v)

    # Get axis
    uvw = np.cross( unit, i_v)

    # compute trig values - no need to go through arccos and back
    rcos = np.dot(i_v, unit)
    rsin = np.linalg.norm(uvw)

    #normalize and unpack axis
    if not np.isclose(rsin, 0):
        uvw /= rsin
    u, v, w = uvw

    # Compute rotation matrix - re-expressed to show structure
    return (
        rcos * np.eye(3) +
        rsin * np.array([
            [ 0, -w,  v],
            [ w,  0, -u],
            [-v,  u,  0]
        ]) +
        (1.0 - rcos) * uvw[:,None] * uvw[None,:]
    )

def direction(ca,chains,ids):
    chains=chains
    xyz=np.array(ca[['x','y','z']].tolist(),dtype="float64")
    m=np.float(12.0107)
    totm=m*xyz.shape[0]
    print(totm)
    coms=com(xyz,m,totm)
    I=inertia(xyz,coms,m)
    #inertia = np.dot( xyz.transpose(),xyz)
    #print(inertia)
    xyz=xyz-coms
    e_values, e_vectors = np.linalg.eig(I)#inertia)
    print(I)
    order=np.argsort(e_values)
    axis3,axis2,axis1 =  e_vectors[:,order].transpose()
    #rotxyz=np.array([axis1,axis2,axis3])
    #rotxyz=np.dot(np.array([[1,0,0],[0,1,0],[0,0,1]]),rotxyz)
    #axess=np.array([1,0,0])
    print(I.shape)
    rotxyz=get_rotation_matrix(axis2,unit=[0.0,0.0,1.0])
    #print(rotxyz)
    for i in range(ca.shape[0]):
        xyz[i]=np.dot(xyz[i].T, rotxyz.T)
        #xyz[i]=np.dot(xyz[i],axess)
    #    xyz[i]=np.dot(rotxyz,xyz[i])
    print(xyz[0])
    #xyz=xyz-xyz.min(axis=0)
    print(xyz[0])
    for i in range(ca.shape[0]):
        ca['x'][i]=round(xyz[i][0],3)
        ca['y'][i]=round(xyz[i][1],3)
        ca['z'][i]=round(xyz[i][2],3)
    wp.writeca(ca,"aln_%s" %ids)

    coms=[]
    for ch in sorted(chains):
        loc=np.where(ca['ch']==ch)
        chcom=np.mean(xyz[loc],0)
        coms.append(chcom)
        #np.mean(np.array(ca[ca['ch']==ch][['x','y','z']].tolist()),0))
    
    edges=[]
    for i in range(len(coms)):
        if i == len(coms)-1:
            nextc=coms[0]
        else:
            nextc=coms[i+1]
        c=coms[i]
        edges.append((nextc[0]-c[0])*(nextc[1]+c[1]))
    mult=1
    if xyz[0][2] > np.mean(xyz,0)[2]:
        mult=-1
    else:
        mult=1

    print(sum(edges))
    s=sum(edges)*mult
    print(s)
    if s < 0:
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
        directions[ids]=direction(uni,chorderlist[ids],ids)
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
        uni_list[ids]=pP.parse_ca(open(cwd+'/'+name,'r'),'NA',mer,altloc,id_ch_dict[ids])
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



