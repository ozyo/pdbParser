import numpy as np

def direction(ca):
    chains=set(ca['ch'].flatten())
    coms=[]
    for ch in sorted(chains):
        coms.append=list(np.mean(ca[ca['ch']==ch][['x','y','z']]).flatten())
    coms=np.array(coms)
    center=np.mean(coms)
    np.linalg.det(np.stack(coms[0:2],center))