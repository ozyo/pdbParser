import itertools
import numpy as np

def renumber(pdb:np.recarray)->np.recarray:
    numbers = list(pdb[["icode","resnr","ch"]])
    counts = [ sum( 1 for _ in group ) for key, group in itertools.groupby( numbers ) if key ]
    switch_chains = [ sum( 1 for _ in group ) for key, group in itertools.groupby( list(pdb["ch"]) ) if key ]
    new_numbers = []
    i = 1
    for count in counts:
        for x in range(count):
            new_numbers.append(i)
            if i == switch_chains[0]:
                switch_chains.pop(0)
                i = 1
            else:
                i+=1
    pdb["resnr"]=new_numbers
    return pdb