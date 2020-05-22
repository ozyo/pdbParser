#See COPYING for license 

import numpy as np
import logging

def getch(coord,chlist,returnch=False):
    _, idx = np.unique(coord['ch'], return_index=True)
    chains=coord['ch'][np.sort(idx)]
    flatchlist=[]
    for mol in chlist:
        for ch in mol:
            flatchlist.append(ch)
    ordch=[]
    for ch in chains:
        if ch in flatchlist:
            ordch.append(ch)
    delch=coord[np.in1d(coord['ch'],ordch)]
    if returnch:
        return delch, ordch
    return delch

def getca(coord,altloc,chlist=[],returnch=False):
    delalter=coord[(coord['altloc'] == altloc) | (coord['altloc'] == '')]
    ca=delalter[np.in1d(delalter['atname'],'CA')]
    if len(chlist) > 0:
        if returnch:
            ca,ordch=getch(ca,chlist,returnch=True)
            return ca,ordch
        else:
            ca=getch(ca,chlist,returnch=False)
            return ca

def clean_altloc(coord,altloc):
    return coord[(coord['altloc'] == altloc) | (coord['altloc'] == '')]
