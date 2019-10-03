from modeller import *
from modeller.automodel import *
log.verbose()
env = environ()
env.io.atom_files_directory = ['.']

mdl = model(env, file='4TNV_fix_BA2.B99990001.pdb')
mdl.rename_segments(segment_ids=('P', 'Q','R','S','T'),renumber_residues=[1,1,1,1,1])
mdl.write(file='4TNV_fix_BA2_rnr.pdb')

mdl = model(env, file='4TNV_fix_BA1.B99990001.pdb')
mdl.rename_segments(segment_ids=('A', 'B','C','D','E'),renumber_residues=[1,1,1,1,1])
mdl.write(file='4TNV_fix_BA1_rnr.pdb')
