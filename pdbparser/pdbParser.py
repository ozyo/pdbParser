#See COPYING for license 

from urllib2 import urlopen
import readpdb, clean_pdb, divide_mer
import missing, writepdb
import logging
from sep_seg import Segsep

def pdbParser(pdb,pdbid,mer,cwd,altloc,rname,bysegid):
    print pdbid
    logging.info('Checking for missing residues')
    logging.info('Retriving CA coordinates')
    atomlines=readpdb.readall(pdb)
    coords=readpdb.coord(atomlines)
    ca=clean_pdb.getall(coords,altloc)
    segs=Segsep(ca)
    final_segs=segs.sep_segs(ca,cwd,rname)

    #This is the part where we check the missing residues. I have a feeling we should do this before. But if we want a sequence alignment between two structure and retrive a region automatically for eBDIMS then doing it after is better so that we can also use the remarks and seqres to chop and align the sequences.
    return ca
    
def pdbParselocal(pdb,cwd,charmm,altloc,rname,bysegid):
    if charmm is False:
        compnd=None
        atomlines=readpdb.readatom(pdb)
        coords=readpdb.coord(atomlines)
        ca=clean_pdb.getall(coords,altloc)
        return ca
    else:
        compnd=None
        atomlines=readpdb.readall(pdb)
        coords=readpdb.coordcharm(atomlines)
        ca=clean_pdb.getall(coords,altloc)
        segs=Segsep(ca)
        if bysegid is True:
            final_segs=segs.sep_bysegids(ca,cwd)
        else:
            final_segs=segs.sep_segs(ca,cwd,rname)
        return final_segs
