#We need to seperate the biological assemblies
#
import numpy as np

def divide_mer(coord,compnd,mer):
    slice=len(compnd) / mer
    divide_mer=np.array_split(coord,slice)
    print divide_mer[-1]['ch']
    return divide_mer[0]

