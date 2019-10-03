from modeller import *
from modeller.automodel import *
log.verbose()
env = environ()
env.io.atom_files_directory = ['.']

mdl = model(env, file='6D6U_chE_repair.pdb')
mdl.rename_segments(segment_ids=('E'),renumber_residues=[25])
mdl.write(file='6D6U_chE_rnr.pdb')

