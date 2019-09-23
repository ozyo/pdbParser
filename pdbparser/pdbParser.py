#See COPYING for license 

from pdbParser import readpdb, clean_pdb
import logging

def pdbTitle(pdb):
    titles=[]
    for line in pdb:
        if 'TITLE' in line:
            titles.append(line)
    for i in titles:
        if 'CHIMERA' in i or 'FUSED' in i:
            return True 
        else:
            return False

def pdbParser(pdb,pdbid,mer,altloc,chlist):
    logging.info('Retriving CA coordinates')
    logging.info('Checking for missing residues')
    atomlines=readpdb.readatom(pdb)
    coords=readpdb.coord(atomlines)
    print(coords)
    ca=clean_pdb.getca(coords,altloc,chlist)
    print(ca)
    return ca

    #Above loop is the part where I check the missing residues. I have a feeling we should do this before. But if we want a sequence alignment between two structure and retrive a region automatically for eBDIMS then doing it after is better so that we can also use the remarks and seqres to chop and align the sequences.
