#See COPYING for license 

from urllib2 import urlopen
import readpdb, clean_pdb, divide_mer
import missing, writepdb
import logging

def pdbParser(pdb,pdbid,mer):
    print pdbid
    compnd=readpdb.readcompnd(pdb)
    logging.info('Detected chains for %s are ' % (pdbid)+' '.join(i for i in compnd))
    if len(compnd) is not 1:
        logging.info('If multimeric option not chosen continuing with the chain %s' %(compnd[0]))
    if mer != 1:
        logging.info('Multimeric option chosen...')
        if mer < len(compnd):
            logging.warning('%s structure contains more biological assemblies than one' % (pdbid))
            if len(compnd) % mer != 0:
                logging.critical('Cannot determine the biological assembly, please provide your own input files')
                exit()
            else:
                logging.warning('Will attempt to seperate them')
        elif mer == 1 and len(compnd) > 1:
            logging.warning('You chose only 1 chain for the analysis but there are more chains in the biological assembly.')
            logging.warning('Assuming monomeric assembly.')
        elif mer > len(compnd):
            logging.warning("Structure contains less chains then the supplied number. Assuming not complete structure, please provide your own input files")
            exit()
        elif mer == len(compnd):
            pass
    logging.info('Checking for missing residues')
    logging.info('Retriving CA coordinates')
    atomlines=readpdb.readatom(pdb)
    coords=readpdb.coord(atomlines,compnd)
    ca=clean_pdb.getca(coords,compnd)
    #This is the part where we check the missing residues. I have a feeling we should do this before. But if we want a sequence alignment between two structure and retrive a region automatically for eBDIMS then doing it after is better so that we can also use the remarks and seqres to chop and align the sequences.
    r465=readpdb.readremark(pdb,compnd)
    misinfo=missing.missinginfo(r465,compnd,ca)
    global div
    if misinfo is not None:
        pass
    else:
        logging.critical('All chains are broken, provide your input files.')
        exit()
    if len(compnd) > mer:
        div=divide_mer.divide_mer(ca,compnd,mer,misinfo)
    else:
        div=ca
    return div
