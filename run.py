#!/usr/bin/env python

from pdbparser.pdbParser import pdbParser
from pdbparser.readpdb import getpdb
from pdbparser.readpdb import checkmulti
from pdbparser.writepdb import writeca
from align.alignment import getaligned, multialigned
import argparse
import logging
from os import getcwd

parser = argparse.ArgumentParser(description='Identification of Missing residues')
parser.add_argument('--start', metavar='PDB File', nargs=1 , help='Starting structure')
parser.add_argument('--target', metavar='PDB File', nargs=1 , help='Target structure')
parser.set_defaults(local=False)
parser.add_argument('--multimeric', dest='mer', help='If the protein is multimeric provide the number of chains. Default is 1', type=int)
parser.set_defaults(mer=1)
parser.add_argument('--dir', dest='cwd', help='The directory to save the output files. Default is current work directory.')
parser.set_defaults(cwd=getcwd())
parser.add_argument('--altloc', dest='altloc', help='Alternative location to be extracted, default is A')
parser.set_defaults(altloc='A')

args=parser.parse_args()
try:
    sid=args.start[0]
    eid=args.target[0]
except TypeError:
    parser.print_help()
    exit()
print args.cwd

outf=args.cwd+'/error.dat'
logging.basicConfig(level=logging.INFO,filename=outf)

logging.info('Fetching PDB files')
start=getpdb(sid)
checkmulti(start)
logging.info('Processing PDB files')
sca=pdbParser(start,sid,args.mer,args.altloc)
logging.info('Retriving CA coordinates for starting structure successful')

logging.info('Fetching PDB files')
end=getpdb(eid)
checkmulti(end)
logging.info('Processing PDB files')
eca=pdbParser(end,eid,args.mer,args.altloc)
logging.info('Retriving CA coordinates for the target structure successful')

toAlign=True

if toAlign is True and args.mer == 1:
    logging.info('Extracting the core region.')
    score,ecore=getaligned(sca,eca)
    print 'here'
    writeca(score,args.cwd+'/start.pdb')
    writeca(ecore,args.cwd+'/target.pdb')
elif toAlign is True and args.mer !=1:
    logging.info('Extracting the core region.')
    score,ecore=multialigned(sca,eca,args.mer)
    writeca(score,args.cwd+'/start.pdb')
    writeca(ecore,args.cwd+'/target.pdb')
