#This is the main script where we call
from urllib2 import urlopen
from readpdb import getpdb, readcompnd
from clean_pdb import *
#from alignment import *
import argparse
import logging
parser = argparse.ArgumentParser(description='Identification of Missing residues')

parser.add_argument('--start', metavar='PDB Code', nargs=1 , help='Crystal structure from PDB database with the header')
parser.add_argument('--end', metavar='PDB Code', nargs=1 , help='Crystal structure from PDB database with the header')
parser.add_argument('--local',help='User PDB files',dest='local', action='store_true')
parser.add_argument('--no-local',help='Files are fetched from RCSB database',dest='local',action='store_false')
parser.set_defaults(local=False)
parser.add_argument('--multimeric', dest='mer', help='If the protein is multimeric provide the number of mers.', type=int)
parser.set_defaults(mer=1)

args=parser.parse_args()
sid=args.start[0]
eid=args.end[0]

if args.local is True:
    sca=getca(open(args.start[0]).readlines())
    eca=getca(open(args.end[0]).readlines())
    logging.info('You have provided PDB files, assuming you have fixed the missing residues and made sure that the structures have the same number of residues')
    logging.info('Moving to eBDIMS calculations')
else:
    logging.info('Fetching PDB files from RCSB database')
    start=getpdb(sid)
    end=getpdb(eid)
    logging.info('Processing PDB files')
    

#identify missing residues
#Needs two routines, one from PDB Header and the other one from sequence when Header is not possible 

if args.local is True:
    logging.info('Moving to eBDIMS calculations')
else:
    sca=getca(start,sid,args.mer)
    eca=getca(end,eid,args.mer)
    logging.info('Checking for missing residues')
#create sequences from CA 
#We still have a problem of multiple confomers. They are not present in the file but we have not renamed them either. So the shell script below requires more lines to replace them.
#     call(['./check_seq.sh end.pdb'],shell=True)
#     call(['./check_seq.sh start.pdb'],shell=True)
#     eseq=open('end.seq','r').readlines()[0].strip()
#     sseq=open('start.seq','r').readlines()[0].strip()
#Grab CA and create an alignment
#     align(eseq,sseq)
#Check for misaligned regions, report any residue mismatch.

#Pass the inputs as start.pdb and end.pdb to the eBDIMS
