#See COPYING for license 

from pdbParser import readpdb as rp
from pdbParser import clean_pdb as cp
from importlib import reload
import numpy as np
reload(rp)
reload(cp)
import logging

def pdbTitle(pdb):
    titles=[]
    for line in pdb:
        if 'TITLE' in line:
            titles.append(line)
    for i in titles:
        if any(True for x in ['fused','chimeric','chimera','chimaeric'] if x in i.lower()):
            return True 
        else:
            return False

def parse_ca(pdb,pdbid,mer,altloc,chains=[],caonly=True,returnch=False):
    logging.info('Retriving CA coordinates')
    logging.info('Checking for missing residues')
    atomlines=rp.readatom(pdb)
    coords=rp.coord(atomlines)
    if caonly:
        if returnch:
            coords,ordch=cp.getca(coords,altloc,chlist=chains,returnch=True)
            return coords,ordch
        else:
            coords=cp.getca(coords,altloc,chlist=chains,returnch=False)
            return coords
    else:
        if returnch and len(chains) > 0:
            coords,ordch=cp.getch(coords,chains)
            return coords,ordch
        else:
            return coords
            
def parse_chlist(pdb):
    logging.info('Parsing the chain information only')
    atomlines=rp.readatom(pdb)
    coords=rp.coord(atomlines)
    _, idx = np.unique(coords['ch'], return_index=True)
    chains=coords['ch'][np.sort(idx)]
    return chains