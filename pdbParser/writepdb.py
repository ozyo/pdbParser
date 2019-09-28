#See COPYING for license 

import numpy as np

def writeca(div,file):
    np.savetxt(file,div,fmt='ATOM  %5s  %2s %1s%3s %1s%4s%1s   %8s%8s%8s%6s%6s          %-2s%2s',delimiter='')
