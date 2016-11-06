#This is the main script where we call
from subprocess import call
from get_missing import *
import argparse
parser = argparse.ArgumentParser(description='Identification of Missing residues')
#This should be a redirect from getpdb.sh command, but for now this will do
parser.add_argument('--start', metavar='Path to pdb file', nargs=1 , help='Crystal structure from PDB database with the header')
parser.add_argument('--end', metavar='Path to pdb file', nargs=1 , help='Crystal structure from PDB database with the header')
#maybe it is a good idea to have a fixed name for this. Users can rename it themselves later if they want. Tag should be PDBID_fix

parser.add_argument('--multimeric', metavar='Provide how many chains are present in this protein', nargs=1, help='If the protein is multimeric provide the number of mers.')

args=parser.parse_args()
call(["./get-fix_pdb.sh", args.start[0]])
call(["./get-fix_pdb.sh", args.end[0]])
start=open(args.start[0]+'.pdb').readlines()
end=open(args.end[0]+'.pdb').readlines()


#identify missing residues
#Needs two routines, one from PDB Header and the other one from sequence when Header is not possible 
getmissing=getmissing()
sdata=getmissing.getmissing(start)
edata=getmissing.getmissing(end)

print spdb[1]
print epdb[1]

#clean pdb files and return CA

#Grab CA and create an alignment

#Check for misaligned regions, report any residue mismatch.

#Pass the inputs as start.pdb and end.pdb to the eBDIMS
