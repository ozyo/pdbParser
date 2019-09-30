import math
import numpy as np
from copy import deepcopy
from pdbParser import parser as pP
from pdbParser import writepdb as wp

def com(ca,chain=None):
    #This is not really a center of mass calculation, no mass weight is considered
    if chain:
        tmp=ca[ca['ch']==chain]
    else:
        tmp=ca
    com=[np.mean(np.array(tmp[['x']],dtype="float64")),np.mean(np.array(tmp[['y']],dtype="float64")),np.mean(np.array(tmp[['z']],dtype="float64"))]
    return(com)

def orient(tca,ids,coms_only=True):
    tca=deepcopy(tca)
    x=deepcopy(tca['x'])
    y=deepcopy(tca['y'])
    z=deepcopy(tca['z'])
    xx=deepcopy(tca['x'])
    yy=deepcopy(tca['y'])
    zz=deepcopy(tca['z'])
    lens=tca.shape[0]
    m=np.array([12.0107]*lens) #carbon
    comx=0
    comy=0
    comz=0
    totalm=0
    for i in range(lens):
        # use the abs of the weights                                                                                                                                                                            
        mm=abs(m[i])
        totalm+=mm
        comx+=x[i]*mm
        comy+=y[i]*mm
        comz+=z[i]*mm
    comx=comx/totalm
    comy=comy/totalm
    comz=comz/totalm
    coms=[comx,comy,comz]
    if coms_only:
        return(coms)
    #get inertia
    Ixx=0
    Ixy=0
    Ixz=0
    Iyy=0
    Iyz=0
    Izz=0
    for i in range(lens):
        mm=abs(m[i])
        totalm+=mm
        xx[i]=xx[i]-coms[0]
        yy[i]=yy[i]-coms[1]
        zz[i]=zz[i]-coms[2]
        rr=xx[i]+yy[i]+zz[i]
        Ixx=Ixx+mm*(yy[i]*yy[i]*zz[i]*zz[i])
        Ixy=Ixy-mm*(xx[i]*yy[i])
        Ixz=Ixz-mm*(xx[i]*zz[i])
        Iyy=Iyy+mm*(xx[i]*xx[i]*zz[i]*zz[i])
        Iyz=Iyz-mm*(yy[i]*zz[i])
        Izz=Izz+mm*(xx[i]*xx[i]*yy[i]*yy[i])
    I=np.array([[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]],dtype="float64")
    #get vectors
    e_values, e_vectors=np.linalg.eig(I)
    #order = np.argsort(e_values)
    #e_values = e_values[order]
    #print(e_vectors[:, order])
    #e_vectors = e_vectors[:, order].transpose()
    #print(e_vectors)
    for i in range(lens):
        xyz=np.array([xx[i],yy[i],zz[i]])
        t=np.dot(e_vectors,xyz)
        xx[i]=t[0]
        yy[i]=t[1]
        zz[i]=t[2]
    if max(xx)-min(xx) > max(yy)-min(yy):
        tmpx=xx
        xx=yy
        yy=tmpx
    elif max(yy)-min(yy) > max(zz)-min(zz):
        tmpy=yy
        yy=zz
        zz=tmpy
    #print(ca['x'][0])
    for i in range(lens):
        tca['x'][i]=round(xx[i],3)
        tca['y'][i]=round(yy[i],3)
        tca['z'][i]=round(zz[i],3)
    wp.writeca(tca,"aln_%s" %ids)
    return(tca)

def direction(ca,chains,ids):
    chains=chains
    com_all=orient(ca,ids,coms_only=True)
    #print(com_all)
    #if ca['z'][0] > com_all[2] or ca['y'][0] > com_all[1] or ca['x'][0] > com_all[0]:
    #    for i in range(ca.shape[0]):
    #        ca['x'][i]=ca['x'][i]*-1
    #        ca['y'][i]=ca['y'][i]*-1
    #        ca['z'][i]=ca['z'][i]*-1
    #Making sure all of them looks at the same direction in the coordinate axes, otherwise assignments are wrong.
    xyz=np.array(ca[['x','y','z']].tolist(),dtype="float64")
    xyz=xyz-com_all
    #np.save("test_%s.npy" %ids,ca[['atnr','ch','x','y','z']])
    # For ordering https://github.com/pierrepo/principal_axes/blob/master/principal_axes.py
    inertia = np.dot(xyz.transpose(), xyz)
    #print(inertia)
    e_values, e_vectors = np.linalg.eig(inertia)
    order = np.argsort(e_values)
    eval3, eval2, eval1 = e_values[order]
    #print(e_vectors[:,order])
    rotxyz = e_vectors[:, order].transpose()
    #print(e_values[order])
    #print(rotxyz)
    #rotxyz=np.array([axis1,axis2,axis3])
    for i in range(ca.shape[0]):
        t=np.dot(rotxyz,xyz[i])
        xyz[i][0]=t[0]
        xyz[i][1]=t[1]
        xyz[i][2]=t[2]
    #print(np.mean(xyz,0))
    if max(xyz[:,0])-min(xyz[:,0]) > max(xyz[:,1])-min(xyz[:,1]):
        tmpx=xyz[:,0]
        xyz[:,0]=xyz[:,1]
        xyz[:,1]=tmpx
    elif max(xyz[:,1])-min(xyz[:,1]) > max(xyz[:,2])-min(xyz[:,2]):
        tmpy=xyz[:,1]
        xyz[:,1]=xyz[:,2]
        xyz[:,2]=tmpy
    for i in range(ca.shape[0]):
        ca['x'][i]=round(xyz[i][0],3)
        ca['y'][i]=round(xyz[i][1],3)
        ca['z'][i]=round(xyz[i][2],3)
    #wp.writeca(ca,"aln_%s" %ids)
    #com_all=np.array(com(ca))
    #if ca['z'][0] > com_all[2]:
    #    for i in range(ca.shape[0]):
    #        ca['x'][i]=ca['x'][i]*-1
    #        ca['y'][i]=ca['y'][i]*-1
    #        ca['z'][i]=ca['z'][i]*-1
    #    wp.writeca(ca,"rot_%s" %ids)
        # xyz=np.array(ca[['x','y','z']].tolist())
        # angle=180 * math.pi / 180
        # zn=xyz[:,2]
        # yn=xyz[:,0]*math.cos(angle) - xyz[:,1]*math.sin(angle)
        # xn=xyz[:,0]*math.sin(angle) + xyz[:,1]*math.cos(angle)
        # for i in range(xyz.shape[0]):
        #     ca['x'][i]=round(xn[i],3)
        #     ca['y'][i]=round(yn[i],3)
        #     ca['z'][i]=round(zn[i],3)
        #wp.writeca(ca,"rot_%s" %ids)
    #ca=orient(ca,ids)
    coms=[]
    for ch in chains:
        coms.append(orient(ca[ca['ch']==ch],ids,coms_only=True))
    edges=[]
    for ind,c in enumerate(coms):
        if ind == len(chains)-1:
            nextc=coms[0]
        else:
            nextc=coms[ind+1]
        edges.append((nextc[2]-c[2])*(nextc[1]+c[1]))
    if sum(edges) > 0:
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



