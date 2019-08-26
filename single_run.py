#!/usr/bin/env python

from pdbParser.pdbParser import pdbParser
from pdbParser.readpdb import getpdb
from pdbParser.readpdb import checkmulti
from pdbParser.writepdb import write
#from align.alignment import getaligned, multialigned
from pdbParser.pdbParser import pdbParselocal
import argparse
import logging
from os import getcwd
from pdbParser.sep_seg import Segsep


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
parser.add_argument('--rname', dest='rname',help='Renaming of the water segnames, and bases to Charmm values default is False',action='store_true')
parser.set_defaults(rname=False)
parser.add_argument('--namena', dest='namena',help='Renaming of bases with 3, 2 or 1 letter code. 3l is ADE, 2l is DA, 1l is A.',nargs=1)
parser.set_defaults(namena='0l')
parser.add_argument('--bysegid', dest='bysegid',help='Write each segment without chain identifier to a seperate file',action='store_true')
parser.set_defaults(bysegid=False)
parser.add_argument('--bfact', dest='bfact', help='Add Bfactors, i.e. for restraints in NAMD',action='store_true')
parser.set_defaults(bfact=False)
parser.add_argument('--bfact_nr', dest='bfact_nr',help='Force constant',nargs=1)
parser.set_defaults(bfact_nr=3.0)
parser.add_argument('--bfact_type', dest='bfact_type',help='all,heavy',nargs=1)
parser.set_defaults(bfact_type='heavy')
parser.add_argument('--out',dest='out',help='outputname',type=str)
parser.add_argument('--merge',dest='merge',help='Merging the segment numbering together, default is false',action='store_true')
parser.set_defaults(merge=False)
parser.add_argument('--drude',dest='drude',help='Remove drudes and lonepairs from the output, no renumbering is done!',action='store_true')
parser.set_defaults(drude=False)
parser.add_argument('--addchid',dest='chid',help='Adding chainids to each segment, starts from A',action='store_true')
parser.set_defaults(chid=False)
parser.add_argument('--reord',dest='reord',help='Reorder segments',action='store_true')
parser.set_defaults(reord=False)
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
    sca=pdbParselocal(open(args.start[0]).readlines(),args.cwd,rcoortype,args.altloc,args.rname,args.bysegid,args.renumber,args.segid,args.bfact,args.merge,args.chains,args.drude,args.bfact_nr,args.bfact_type,args.namena,args.chid,args.reord)
    if rcoortype == wcoortype:
        write(sca,args.cwd+'/'+out,wcoortype,args.bysegid)
    else:
        #For now this is broken but there should be functions to convert the array between coord types.
        #Attempt in fixing this. One might want to add the segids before processing file, support for that will come later
        if wcoortype == 'charmm' and args.bysegid is False:
            segs=Segsep(sca)
            sca=segs.add_segid(sca)
            write(sca,args.cwd+'/'+out,wcoortype,args.bysegid)
        else:
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
