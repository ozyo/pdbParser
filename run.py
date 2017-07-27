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
parser.add_argument('--schains',dest='schains',help='Supply the chain ids of the starting structure: A-E or A,B,C')
parser.set_defaults(schains='A')
parser.add_argument('--echains',dest='echains',help='Supply the chain ids of the target structure: A-E or A,B,C')
parser.set_defaults(echains='A')
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

def chlistsplit(chlist):
    finlist=[]
    tmplist=chlist.split(",")
    for sub in tmplist:
        if len(sub)==3 and sub[1]=="-":
            sublist=sub.split("-")
            firstch=sublist[0]
            lastch=sublist[1]
            finsublist=[chr(i) for i in range(ord(firstch),ord(lastch)+1)]
            for ch in finsublist:finlist.append(ch)
        elif len(sub)==1 and isinstance(sub, str):
            finlist.append(sub)
        else:
            logging.critical("Please supply the chain ids in these formats: A-E or A-C,E or A,B")
            exit()
    return finlist

schains=chlistsplit(args.schains)
echains=chlistsplit(args.echains)

toAlign=True

outf=args.cwd+'/error.dat'
logging.basicConfig(level=logging.INFO,filename=outf)
if args.mer == 1:
    for ch1,ch2 in zip(schains,echains):
        start=getpdb(sid)
        checkmulti(start)
        logging.info('Processing PDB files')
        sca=pdbParser(start,sid,args.mer,args.altloc,ch1)
        end=getpdb(eid)
        checkmulti(end)
        logging.info('Processing PDB files')
        eca=pdbParser(end,eid,args.mer,args.altloc,ch2)
        logging.info('Extracting the core region.')
        score,ecore,correct=getaligned(sca,eca)
        if correct is True:
            writeca(score,args.cwd+'/start.pdb')
            writeca(ecore,args.cwd+'/target.pdb')
            break
        else:
            continue

elif args.mer !=1:
    start=getpdb(sid)
    checkmulti(start)
    logging.info('Processing PDB files')
    sca=pdbParser(start,sid,args.mer,args.altloc,schains)
    end=getpdb(eid)
    checkmulti(end)
    logging.info('Processing PDB files')
    eca=pdbParser(end,eid,args.mer,args.altloc,echains)
    logging.info('Extracting the core region.')
    score,ecore,correct=multialigned(sca,eca,args.mer)
    if correct is True:
        writeca(score,args.cwd+'/start.pdb')
        writeca(ecore,args.cwd+'/target.pdb')

