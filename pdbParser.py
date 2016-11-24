###DEAL with the license options laters!!!!


#This is the main script where we call
from urllib2 import urlopen
from readpdb import *
from clean_pdb import *
from divide_mer import *
from missing import *
from writepdb import *
from alignment import getaligned
import argparse
import logging

logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser(description='Identification of Missing residues')
parser.add_argument('--start', metavar='PDB Code', nargs=1 , help='Crystal structure from PDB database with the header')
parser.add_argument('--end', metavar='PDB Code', nargs=1 , help='Crystal structure from PDB database with the header')
parser.add_argument('--local',help='User PDB files',dest='local', action='store_true')
parser.add_argument('--no-local',help='Files are fetched from RCSB database',dest='local',action='store_false')
parser.set_defaults(local=False)
parser.add_argument('--multimeric', dest='mer', help='If the protein is multimeric provide the number of mers. Default is 1', type=int)
parser.set_defaults(mer=1)

args=parser.parse_args()
sid=args.start[0]
eid=args.end[0]

def pdbParser(pdb,pdbid,mer,file):
    print pdbid
    compnd=readcompnd(pdb)
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
            logging.warning("Structure contains less chains then the supplied multimeric option. Assuming not complete structure, please provide your own input files")
            exit()
        elif mer == len(compnd):
            pass
    logging.info('Checking for missing residues')
    logging.info('Retriving CA coordinates')
    atomlines=readatom(pdb)
    coords=coord(atomlines,compnd)
    ca=getca(coords,compnd)
    #This is the part where we check the missing residues. I have a feeling we should do this before. But if we want a sequence alignment between two structure and retrive a region automatically for eBDIMS then doing it after is better so that we can also use the remarks and seqres to chop and align the sequences.
    r465=readremark(pdb,compnd)
    missing=missinginfo(r465,compnd,ca)
    global div
    if missing is not None:
        pass
    else:
        logging.critical('All chains are broken, provide your input files.')
        exit()
    if len(compnd) > mer:
        div=divide_mer(ca,compnd,mer,missing)
        #writeca(div,file)
    else:
        div=ca
        #writeca(div,file)
    return div
global toAlign
if args.local is True:
    #We want to check the integrigty of the files anyways. It is possible that the user will provide anything else than CA atoms or missmatching files.
    #Deal with this later since those files will not have the header info
    logging.info('You have provided PDB files, assuming you have fixed the missing residues and made sure that the structures have the same number of residues')
    sca=getca(open(args.start[0]).readlines())
    eca=getca(open(args.end[0]).readlines())
    toAlign=False
else:
    logging.info('Fetching PDB files from RCSB database')
    start=getpdb(sid)
    end=getpdb(eid)
    logging.info('Processing PDB files')
    sca=pdbParser(start,sid,args.mer,'start.pdb')
    eca=pdbParser(end,eid,args.mer,'end.pdb')
    logging.info('Retriving CA coordinates successful')
    toAlign=True

if toAlign is True:
    getaligned(sca,eca)
