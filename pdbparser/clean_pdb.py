#See COPYING for license 

import numpy as np
import logging

def getca(coord,altloc):
    delalter=coord[(coord['altloc'] == altloc) | (coord['altloc'] == '')]
    logging.warning('Cleaning alternative locations if present')
    logging.warning('Default location is A')
    ca=delalter[np.in1d(delalter['atname'],'CA')]
    return ca
