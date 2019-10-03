from modeller import *
from modeller.automodel import *
from glob import glob
log.verbose()
env = environ()
env.io.atom_files_directory = ['.']

cpdb=glob('correct_*pdb')
for p in cpdb:
    name='_'.join(['rnr']+p.split('_')[1:])
    mdl = model(env, file=p)
    mdl.rename_segments(segment_ids=('A', 'B','C','D','E'),renumber_residues=[1,1,1,1,1])
    mdl.write(file=name)

