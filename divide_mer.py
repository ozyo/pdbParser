#We need to seperate the biological assemblies
#
import numpy as np

def divide_mer(ca,compnd,mer,missinginfo):
    slice=len(compnd) / mer
    divide_mer=np.array_split(ca,slice)
    if all(missinginfo.values()==1):
        divide_mer=np.array_split(ca,slice)[0]
        logging.info('I will continue with the first biological assembly')
    else:
        

    return divide_mer[0]
