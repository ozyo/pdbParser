#Things to do
#func 1 identify the multiple rotamers
#clean anything that is not CA
#return CA coordinates.
import numpy as np
from readpdb import coord

#left from old list method. Might be useful in somewhere just leave it 
#aa=['ALA','HIS','ILE','LEU','LYS','MET','PHE','PRO','THR','TYR','VAL','GLU','ASN','ARG','GLY','CYS','SER','TRP','ASP','GLN']
#daa=['A'+i for i in aa]

def getca(pdb,pdbid,mer):
    crdarray=coord(pdb,pdbid,mer)
    delalter=crdarray[np.in1d(crdarray['altloc'],'B',invert=True)]
    print 'Cleaning alternative location B if present'
    print 'Currently no support is provided for chosing a different alternative location'
    return delalter[np.in1d(delalter['atname'],'CA')]
