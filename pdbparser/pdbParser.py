#See COPYING for license 

from urllib2 import urlopen
import readpdb, clean_pdb, divide_mer
import missing, writepdb
import logging

def pdbHeader(pdb,pdbid,mer,compnd):
    for line in pdb:
        if 'REMARK' in line:
            return True

def pdbParser(pdb,pdbid,mer,altloc,chlist):
    print pdbid
    logging.info('Retriving CA coordinates')
    logging.info('Checking for missing residues')
    atomlines=readpdb.readatom(pdb)
    coords=readpdb.coord(atomlines)
    ca=clean_pdb.getca(coords,altloc,chlist)
    return ca

    #Above loop is the part where I check the missing residues. I have a feeling we should do this before. But if we want a sequence alignment between two structure and retrive a region automatically for eBDIMS then doing it after is better so that we can also use the remarks and seqres to chop and align the sequences.
