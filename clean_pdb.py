#Things to do
#func 1 identify the multiple rotamers
#clean anything that is not CA
#return CA coordinates.
import numpy as np


aa=['ALA','HIS','ILE','LEU','LYS','MET','PHE','PRO','THR','TYR','VAL','GLU','ASN','ARG','GLY','CYS','SER','TRP','ASP','GLN']
daa=['A'+i for i in aa]

def getca(pdb):
    calines=[]
    seenresid=[]
    for line in pdb:
            if line.startswith('ATOM') and 'CA' in line:
                resid=line[22:26].strip()
                ch=line[20:21].strip()
                if [ch,resid] not in seenresid:
                    seenresid.append([ch,resid])
                    calines.append(line)
                else:
                    print 'Identifiied multiple confomers for residue %s . Continuing with the first confomer' %(resid)
    for ca in calines:
        for rep in daa:
            if rep in ca:
                calines[calines.index(ca)]=ca.replace(rep,' '+aa[daa.index(rep)])
    return calines

