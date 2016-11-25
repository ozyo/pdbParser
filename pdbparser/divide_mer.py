#See COPYING for license 
import numpy as np
import logging

def divide_mer(ca,compnd,mer,missinginfo):
    slice=len(compnd) / mer
    nrba={}
    for i in range(slice):
        nrba[i]=compnd[i*mer:i*mer+mer]
    batouse=[]
    for ba in nrba:
        bainfo=0
        for ch in nrba[ba]:
          bainfo=bainfo+missinginfo[ch]
        if bainfo==mer:
            batouse.append(ba)
            continue
    if len(batouse) == 0:
        logging.warning('All the assemblies are broken. Provide your own input files')
        exit()
    else:
        divide_mer=np.array_split(ca,slice)[batouse[0]]
        logging.info('I will continue with the %i biological assembly' %(batouse[0]+1))
        return divide_mer
