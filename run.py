#!/usr/bin/env python

from pdbparser.pdbParser import pdbParser
from pdbparser.readpdb import getpdb
from pdbparser.writepdb import writeca
from align.alignment import getaligned, multialigned
from pdbparser.pdbParser import pdbParselocal
import argparse
import logging
from os import getcwd

logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser(description='Identification of Missing residues')
parser.add_argument('--start', metavar='PDB Code', nargs=1 , help='Crystal structure from PDB database with the header')
parser.add_argument('--target', metavar='PDB Code', nargs=1 , help='Crystal structure from PDB database with the header')
parser.add_argument('--local',help='User PDB files',dest='local', action='store_true')
parser.set_defaults(local=False)
parser.add_argument('--multimeric', dest='mer', help='If the protein is multimeric provide the number of mers. Default is 1', type=int)
parser.set_defaults(mer=1)
parser.add_argument('--dir', dest='cwd', help='The directory to save the output files. Default is current work directory.')
parser.set_defaults(cwd=getcwd())

args=parser.parse_args()
try:
    sid=args.start[0]
    eid=args.target[0]
except TypeError:
    parser.print_help()
    exit()

if args.local is True:
    #We want to check the integrigty of the files anyways. It is possible that the user will provide anything else than CA atoms or missmatching files.
    #Deal with this later since those files will not have the header info
    logging.warning('You have provided PDB files, assuming you have fixed the missing residues.')
    logging.info('Reading PDB files, extracting the core region...')
    sca=pdbParselocal(open(args.start[0]).readlines())
    eca=pdbParselocal(open(args.target[0]).readlines())
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

if toAlign is True and args.mer == 1:
    logging.info('Extracting the core region.')
    score,ecore=getaligned(sca,eca)
    writeca(score,'start.pdb')
    writeca(ecore,'target.pdb')
elif toAlign is True and args.mer !=1:
    logging.info('Extracting the core region.')
    score,ecore=multialigned(sca,eca,args.mer)
    writeca(score,'start.pdb')
    writeca(ecore,'target.pdb')
