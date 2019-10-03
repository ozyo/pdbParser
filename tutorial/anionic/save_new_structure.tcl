mol load pdb [lindex $argv 0]

#Residue numbers obtained from 6D6U  with python code below
#for tag,block in info.core_blocks.items():
#    print(tag)
#    print("chain A/B2 ",info.coreresids.loc['6D6U_1.pdb|A|'][block[0]],'to',info.coreresids.loc['6D6U_1.pdb|A|'][block[-1]])
#    print("chain B/A1 ",info.coreresids.loc['6D6U_1.pdb|B|'][block[0]],'to',info.coreresids.loc['6D6U_1.pdb|B|'][block[-1]])
#  print("chain E/G2 ",info.coreresids.loc['6D6U_1.pdb|E|'][block[0]],'to',info.coreresids.loc['6D6U_1.pdb|E|'][block[-1]])

    

set sel [atomselect top "name CA and ((chain A C and resid 12 to 30 33 to 81 88 to 174 181 to 191 192 to 306 313 to 339 ) or (chain B D and resid 15 to 33 36 to 84 91 to 177 186 to 196 197 to 311 318 to 344) or (chain E and resid 27 to 45 48 to 96 103 to 189 196 to 206 207 to 321 328 to 354))"]

$sel writepdb "correct_newGABAaR_bicul.pdb"

quit