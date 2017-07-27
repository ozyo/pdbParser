#See COPYING for license 

import numpy as np
import logging

def getch(coord,chlist):
    print chlist
    delch=coord[np.in1d(coord['ch'],chlist)]
    return delch

def getca(coord,altloc,chlist):
    subcoord=getch(coord,chlist)
    delalter=subcoord[(subcoord['altloc'] == altloc) | (subcoord['altloc'] == '')]
    logging.warning('Cleaning alternative locations if present')
    logging.warning('Default alternative location is A')
    ca=delalter[np.in1d(delalter['atname'],'CA')]
    return ca
