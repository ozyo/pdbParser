#!/usr/bin/env python

from pdbparser.pdbParser import pdbParser
from pdbparser.readpdb import getpdb
from pdbparser.readpdb import checkmulti
from pdbparser.writepdb import write
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
parser.add_argument('--coortype', dest='coortype', help='Read and write formats, charmm or charmm crysol, charmm read format creates segment ids',nargs='+')
parser.set_defaults(coortype=['None','None'])
parser.add_argument('--altloc', dest='altloc', help='Alternative location to be extracted, default is A')
parser.set_defaults(altloc='A')
parser.add_argument('--renumber',dest='renumber',help='renumber Water segment residue numbers or renumber atoms; seg or atnr or chain ids')
parser.set_defaults(renumber=False)
parser.add_argument('--segid',dest='segid',help='Segment identifier to renumber',nargs='+')
parser.add_argument('--chains',dest='chains',help='Chain identifier to renumber',nargs='+')
parser.add_argument('--rname', dest='rname',help='Renaming of the water segnames, default is False',action='store_true')
parser.set_defaults(rname=False)
parser.add_argument('--bysegid', dest='bysegid',help='Write each segment without chain identifier to a seperate file',action='store_true')
parser.set_defaults(bysegid=False)
parser.add_argument('--bfact', dest='bfact', help='Add Bfactors, i.e. for restraints in NAMD',action='store_true')
parser.set_defaults(bfact=False)
parser.add_argument('--out',dest='out',help='outputname',type=str)
parser.add_argument('--merge',dest='merge',help='Merging the segment numbering together, default is false',action='store_true')
parser.set_defaults(merge=False)
parser.add_argument('--drude',dest='drude',help='Remove drudes and lonepairs from the output, no renumbering is done!',action='store_true')
parser.set_defaults(drude=False)
args=parser.parse_args()

try:
    sid=args.start[0]
except TypeError:
    parser.print_help()
    exit()
print args.cwd

if args.out is None:
    out=sid+'_processed.pdb'
else:
    out=args.out

if len(args.coortype) == 1:
    args.coortype.append(args.coortype[0])

outf=args.cwd+'/error.dat'
logging.basicConfig(level=logging.CRITICAL,filename=outf)

if args.local is True:
    #Read the coor
    rcoortype=args.coortype[0]
    wcoortype=args.coortype[1]
    sca=pdbParselocal(open(args.start[0]).readlines(),args.cwd,rcoortype,args.altloc,args.rname,args.bysegid,args.renumber,args.segid,args.bfact,args.merge,args.chains,args.drude)
    write(sca,args.cwd+'/'+out,wcoortype,args.bysegid)
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
