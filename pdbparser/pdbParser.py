#See COPYING for license 

from urllib2 import urlopen
import readpdb, clean_pdb, divide_mer
import missing, writepdb
import logging
from sep_seg import Segsep
from rnr import replace, rnrseg_charmm, renumberatom, addbfact , rnrch, removedrude, addchid

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
    
def pdbParselocal(pdb,cwd,coortype,altloc,rname,bysegid,renumber,segids,bfact,merge,chains,drude,bfact_nr,bfact_type,namena,chid):
    atomlines=readpdb.readatom(pdb)
    coords=readpdb.coord(atomlines,coortype)
    #Cleans alternative location
    ca=clean_pdb.getall(coords,altloc)
    if drude is True:
        ca=removedrude(ca)
    if bysegid is True and coortype == 'pdb':
        #if len(ca['segid']) > 1:
        #    pass
        #else:
        #    segs=Segsep(ca)
        #    ca=segs.sep_segs(ca)
        segs=Segsep(ca)
        ca=segs.sep_segs(ca)
    if rname is True:
        old=['A','C','G','T']
        new=['ADE','CYT','GUA','THY']
        ca=replace(ca,old,new)
        ca=replace(ca,['HOH'],['TIP3'])
        old=['DA','DC','DG','DT']
        new=['ADE','CYT','GUA','THY']
        ca=replace(ca,old,new)
    if chid is True:
        ca=addchid(ca)
    if renumber == 'seg':
        ca=rnrseg_charmm(ca,segids,merge)
    elif renumber == 'atnr':
        ca=renumberatom(ca)
    elif renumber == 'ch':
        ca=rnrch(ca,chains,merge)
    else:
        ca=ca
    if bfact is True:
        ca=addbfact(ca,bfact_nr,bfact_type)

    if namena[0] == '3l':
        old=['A','C','G','T']
        new=['ADE','CYT','GUA','THY']
        ca=replace(ca,old,new)
        old=['DA','DC','DG','DT']
        new=['ADE','CYT','GUA','THY']
        ca=replace(ca,old,new)
    elif namena[0] == '2l':
        old=['A','C','G','T']
        new=['DA','DC','DG','DT']
        ca=replace(ca,old,new)
        old=['ADE','CYT','GUA','THY']
        new=['DA','DC','DG','DT']
        ca=replace(ca,old,new)
    elif namena[0] == '1l':
        new=['A','C','G','T']
        old=['DA','DC','DG','DT']
        ca=replace(ca,old,new)
        old=['ADE','CYT','GUA','THY']
        new=['A','C','G','T']
        ca=replace(ca,old,new)

    return ca
