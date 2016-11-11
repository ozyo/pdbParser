#We need to seperate the biological assemblies
#
import numpy as np

def divide_mer(ca,compnd,mer):
    slice=len(compnd) / mer
    divide_mer=np.array_split(ca,slice)
    return divide_mer[0]
