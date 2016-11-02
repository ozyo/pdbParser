import os, argparse, itertools, operator
parser = argparse.ArgumentParser(description='Identification of Missing residues')
#This should be a redirect from getpdb.sh command, but for now this will do
parser.add_argument('--pdb1', metavar='Path to pdb file', nargs=1 , help='Crystal structure from PDB database with the header')
parser.add_argument('--pdb2', metavar='Path to pdb file', nargs=1 , help='Crystal structure from PDB database with the header')
#maybe it is a good idea to have a fixed name for this. Users can rename it themselves later if they want. Tag should be PDBID_fix
args=parser.parse_args()
pdb1=open(args.pdb[0]).readlines()
pdb2=open(args.pdb[0]).readlines()

def get_missing(pdb):
    misres=[]
    misatom=[]
    missing_res=False
    missing_atm=False
    for line in pdb:
        if "REMARK" in line:
            if "465" in line:
                missing_res=True
            else:
                missing_res=False
            if missing_res is True:
                try:
                    ch=line.split(' ',8)[7].strip()
                    res=line.split(' ',8)[8].strip()
                    if len(res) == 3:
                        tofix=(ch,res)
                        misres.append(tofix)
                except IndexError:
                    pass
    misres_sort_ch=sorted((list(group) for key,group in itertools.groupby(misres,operator.itemgetter(0))),key=operator.itemgetter(2))
    for chain in misres_sort_ch:
        res_list=[]
        for pair in range(0,len(chain)):
            res_list.append(chain[pair][1])
    ch_nr=len(misres_sort_ch)

    return misres_sort_ch, res_list, ch_nr

#A complex generator. The first argument argument in sorted groups the misssing residues list by its chain, and creates a sublist. This sublist is then passed to sort to be sorted based on the resid. This is to ensure we create residue ranges for modelling. 
# In theory it should be sorted but you cannot be too careful when it comes to pdb files.
data_pdb1=get_missing(pdb1)
data_pdb2=get_missing(pdb2)


misres=get_missing(pdb)[0]
res_list=get_missing(pdb)[1]
ch_nr=get_missing(pdb)[2]

print 'There are %s chains in this protein. If not multimeric we will use the chains that are not broken' % (ch_nr)
print 'Missing residues halting... Please provide a pdb file that is fixed, an example modeller script to fix broken molecules can be found :...'

# I am getting tired trying to figure this out!!!
