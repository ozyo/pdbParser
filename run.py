#!/usr/bin/env python

from pdbparser.pdbParser import pdbParser
from pdbparser.readpdb import getpdb
from pdbparser.writepdb import writeca
from align.alignment import getaligned
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

if args.local is True:
    #We want to check the integrigty of the files anyways. It is possible that the user will provide anything else than CA atoms or missmatching files.
    #Deal with this later since those files will not have the header info
    logging.warning('You have provided PDB files, assuming you have fixed the missing residues.')
    logging.info('Reading PDB files, extracting the core region...')
    sca=clean_pdb.getca(open(args.start[0]).readlines())
    eca=clean_pdb.getca(open(args.end[0]).readlines())
    toAlign=True
else:
    logging.info('Fetching PDB files from RCSB database')
    start=getpdb(sid)
    end=getpdb(eid)
    logging.info('Processing PDB files')
    sca=pdbParser(start,sid,args.mer)
    eca=pdbParser(end,eid,args.mer)
    logging.info('Retriving CA coordinates successful')
    toAlign=True

if toAlign is True:
    logging.info('Extracting the core region.')
    score,ecore=getaligned(sca,eca)
    writeca(score,'start.pdb')
    writeca(ecore,'end.pdb')
