import numpy as np

def direction(ca):
    chains=set(ca['ch'].flatten())
    coms=[]
    for ch in chains:
        coms.append(np.mean(ca[ca['ch']==ch][['x','y','z']]))
    coms=np.array(coms)
