#See COPYING for license 

from urllib2 import urlopen
import readpdb, clean_pdb, divide_mer
import missing, writepdb
import logging
from sep_seg import Segsep
from rnr import replace, rnrseg_charmm, renumberatom, addbfact

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
    
def pdbParselocal(pdb,cwd,coortype,altloc,rname,bysegid,renumber,segids,bfact):
    atomlines=readpdb.readatom(pdb)
    coords=readpdb.coord(atomlines,coortype)
    #Cleans alternative location
    ca=clean_pdb.getall(coords,altloc)
    if bysegid is True:
        segs=Segsep(ca)
        ca=segs.sep_segs(ca)
    if rname is True:
        old=['A','C','G','T']
        new=['ADE','CYT','GUA','THY']
        ca=replace(ca,old,new)
        ca=replace(ca,['HOH'],['TIP3'])
    if renumber == 'seg':
        ca=rnrseg_charmm(ca,segids)
        print 'Here'
    elif renumber == 'atnr':
        ca=renumberatom(ca)
    else:
        ca=ca
    if bfact is True:
        ca=addbfact(ca,3.0)
    return ca
