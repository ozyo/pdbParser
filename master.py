#This is the main script where we call
from subprocess import call
from get_missing import *
from clean_pdb import *
from alignment import *
import argparse
import logging
parser = argparse.ArgumentParser(description='Identification of Missing residues')

parser.add_argument('--start', metavar='PDB Code', nargs=1 , help='Crystal structure from PDB database with the header')
parser.add_argument('--end', metavar='PDB Code', nargs=1 , help='Crystal structure from PDB database with the header')
parser.add_argument('--local',help='User PDB files',dest='local', action='store_true')
parser.add_argument('--no-local',help='Files are fetched from RCSB database',dest='local',action='store_false')
parser.set_defaults(local=False)

parser.add_argument('--multimeric', metavar='Provide how many chains are present in this protein', nargs=1, help='If the protein is multimeric provide the number of mers.')

args=parser.parse_args()
if args.local is True:
    start=open(args.start[0]).readlines()
    end=open(args.end[0]).readlines()
    logging.warning('You have provided PDB files, assuming you have fixed the missing residues and made sure that the structures have the same number of residues')
else:
    logging.warning('Fetching PDB files from RCSB database')
    call(["./get-fix_pdb.sh", args.start[0]])
    call(["./get-fix_pdb.sh", args.end[0]])
    start=open(args.start[0]+'.pdb').readlines()
    end=open(args.end[0]+'.pdb').readlines()


#identify missing residues
#Needs two routines, one from PDB Header and the other one from sequence when Header is not possible 

if args.local is True:
    logging.warning('Moving to eBDIMS calculations')
else:
    logging.warning('Checking for missing residues')
    getmissing=getmissing()

    sdata=getmissing.getmissing(start)
    edata=getmissing.getmissing(end)

    sca=getca(start)
    eca=getca(end)

    outf=open('end.pdb', 'w')
    for line in eca:
        outf.write(line)
    outf.close()
    outf=open('start.pdb','w')
    for line in sca:
        outf.write(line)
    outf.close()

#create sequences from CA 
#We still have a problem of multiple confomers. They are not present in the file but we have not renamed them either. So the shell script below requires more lines to replace them.
    call(['./check_seq.sh end.pdb'],shell=True)
    call(['./check_seq.sh start.pdb'],shell=True)
    eseq=open('end.seq','r').readlines()[0].strip()
    sseq=open('start.seq','r').readlines()[0].strip()
#Grab CA and create an alignment
    align(eseq,sseq)
#Check for misaligned regions, report any residue mismatch.

#Pass the inputs as start.pdb and end.pdb to the eBDIMS
