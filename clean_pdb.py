#Things to do
#func 1 identify the multiple rotamers
#clean anything that is not CA
#return CA coordinates.
import numpy as np
import logging
#left from old list method. Might be useful in somewhere just leave it 
#aa=['ALA','HIS','ILE','LEU','LYS','MET','PHE','PRO','THR','TYR','VAL','GLU','ASN','ARG','GLY','CYS','SER','TRP','ASP','GLN']
#daa=['A'+i for i in aa]

def getca(coord,compnd):
    delalter=coord[coord['altloc'] != 'B']
    logging.warning('Cleaning alternative location B if present')
    logging.warning('Currently no support is provided for chosing a different alternative location')
    ca=delalter[np.in1d(delalter['atname'],'CA')]
    return ca
