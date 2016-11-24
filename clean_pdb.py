#Things to do
#func 1 identify the multiple rotamers
#clean anything that is not CA
#return CA coordinates.
import numpy as np
import logging

def getca(coord,compnd):
    delalter=coord[coord['altloc'] != 'B']
    logging.warning('Cleaning alternative location B if present')
    logging.warning('Currently no support is provided for chosing a different alternative location')
    ca=delalter[np.in1d(delalter['atname'],'CA')]
    return ca
