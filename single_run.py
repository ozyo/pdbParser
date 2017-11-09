#!/usr/bin/env python

from pdbparser.pdbParser import pdbParser
from pdbparser.readpdb import getpdb
from pdbparser.readpdb import checkmulti
from pdbparser.writepdb import writeca, writecrysol
#from align.alignment import getaligned, multialigned
from pdbparser.pdbParser import pdbParselocal, pdbParsecrysol
import argparse
import logging
from os import getcwd
from pdbparser.rnr import rnrseg_charmm
from pdbparser.rnr import rnratnr_charmm

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
parser.add_argument('--altloc', dest='altloc', help='Alternative location to be extracted, default is A')
parser.set_defaults(altloc='A')
parser.add_argument('--renumber',dest='renumber',help='renumber Water segment residue numbers or renumber atoms; seg or atnr')
parser.set_defaults(renumber=False)
parser.add_argument('--segid',dest='segid',help='Segment identifier to renumber')
parser.add_argument('--rname', dest='rname',help='Renaming of the water segnames, default is False',action='store_true')
parser.set_defaults(rname=False)
parser.add_argument('--bysegid', dest='bysegid',help='Write each segment without chain identifier to a seperate file',action='store_true')
parser.set_defaults(segid=False)
parser.add_argument('--crysol', dest='crysol', help='Write the pdb in correct format with 3 character atomnames',action='store_true')
parser.set_defaults(segid=False)
parser.add_argument('--out',dest='out',help='outputname')
args=parser.parse_args()
try:
    sid=args.start[0]
    parser.set_defaults(out=sid)
#     eid=args.target[0]
except TypeError:
    parser.print_help()
    exit()
print args.cwd

outf=args.cwd+'/error.dat'
logging.basicConfig(level=logging.CRITICAL,filename=outf)

print args.renumber

if args.local is True:
    logging.warning('You have provided PDB files, assuming you have fixed the missing residues.')
    logging.info('Reading PDB files, extracting the core region...')
    if args.charmm is True:
        sca=pdbParselocal(open(args.start[0]).readlines(),args.cwd,True,args.altloc,args.rname,args.bysegid)
    elif str(args.renumber) == 'seg':
        rnrseg_charmm(open(args.start[0]).readlines(),args.segid,args.cwd)
    elif str(args.renumber) == 'atnr':
        rnratnr_charmm(open(args.start[0]).readlines(),args.cwd)
    elif args.crysol is True:
        sca=pdbParsecrysol(open(args.start[0]).readlines(),args.cwd,args.altloc)
        writecrysol(sca,args.cwd+'/'+str(args.out))
    else:
        sca=pdbParselocal(open(args.start[0]).readlines(),args.cwd,False,args.altloc,args.rname,args.bysegid)
        writeca(sca,args.cwd+'/'+sid+'_clean.pdb')
else:
    logging.info('Fetching PDB files from RCSB database')
    start=getpdb(sid)
    checkmulti(start)
    logging.info('Processing PDB files')
    sca=pdbParser(start,sid,args.mer,args.cwd,args.altloc,args.rname)
    logging.info('Retriving CA coordinates successful')
    if args.charmm is True:
        logging.info('Finished')
    else:
        print args.charmm
        writeca(sca,args.cwd+'/'+sid+'_clean.pdb')
