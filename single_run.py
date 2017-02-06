#!/usr/bin/env python

from pdbparser.pdbParser import pdbParser
from pdbparser.readpdb import getpdb
from pdbparser.readpdb import checkmulti
from pdbparser.writepdb import writeca
#from align.alignment import getaligned, multialigned
from pdbparser.pdbParser import pdbParselocal
import argparse
import logging
from os import getcwd

parser = argparse.ArgumentParser(description='Identification of Missing residues')
parser.add_argument('--start', metavar='PDB Code', nargs=1 , help='Crystal structure from PDB database with the header')
#parser.add_argument('--target', metavar='PDB Code', nargs=1 , help='Crystal structure from PDB database with the header')
parser.add_argument('--local',help='User PDB files',dest='local', action='store_true')
parser.set_defaults(local=False)
parser.add_argument('--multimeric', dest='mer', help='If the protein is multimeric provide the number of mers. Default is 1', type=int)
parser.set_defaults(mer=1)
parser.add_argument('--dir', dest='cwd', help='The directory to save the output files. Default is current work directory.')
parser.set_defaults(cwd=getcwd())
parser.add_argument('--charmm', dest='charmm', help='Write charmm segments',action='store_true')
parser.set_defaults(charmm=False)

args=parser.parse_args()
try:
    sid=args.start[0]
#     eid=args.target[0]
except TypeError:
    parser.print_help()
    exit()
print args.cwd

outf=args.cwd+'/error.dat'
logging.basicConfig(level=logging.CRITICAL,filename=outf)

if args.local is True:
    logging.warning('You have provided PDB files, assuming you have fixed the missing residues.')
    logging.info('Reading PDB files, extracting the core region...')
    if args.charmm is True:
        sca=pdbParselocal(open(args.start[0]).readlines(),args.cwd,True)
    else:
        sca=pdbParselocal(open(args.start[0]).readlines())
else:
    logging.info('Fetching PDB files from RCSB database')
    start=getpdb(sid)
    checkmulti(start)
    logging.info('Processing PDB files')
    sca=pdbParser(start,sid,args.mer,args.cwd)
    logging.info('Retriving CA coordinates successful')
    if args.charmm[0] is True:
        logging.info('Finished')
    else:
        print args.charmm
        writeca(sca,args.cwd+'/'+sid+'_clean.pdb')
