#See COPYING for license 

import numpy as np

def write(coord,file,coortype,bysegid):
    if bysegid is False:
        if coortype == "charmm":
            np.savetxt(file,coord,fmt='ATOM  %5s %-4s%1s%-4s%1s%5s   %8s%8s%8s%6s%6s     %-5s',delimiter='')
        else:
            np.savetxt(file,coord,fmt='ATOM  %5s %-4s%1s%-4s%1s%5s   %8s%8s%8s%6s%6s          %-2s%2s',delimiter='')
    else:
        segids=np.unique(coord['segid'])
        for seg in segids:
            to_write=coord[coord['segid']==seg]
            np.savetxt(seg.lower()+'.pdb',to_write,fmt='ATOM  %5s %-4s%1s%-4s%1s%5s   %8s%8s%8s%6s%6s     %-5s',delimiter='')
