from modeller import *
from modeller.automodel import *    # Load the automodel class


from myloop import MyLoop #MODIFY THIS CLASS TOO!!!!!
#Example taken from :
#https://salilab.org/modeller/manual/node35.html

#log.verbose()
env = Environ()

env.io.atom_files_directory = ['.', '../atom_files']
a = MyLoop(env,
           alnfile  = 'seg.ali',      # alignment filename
           knowns   = '4NPQ_4',               # codes of the templates
           sequence = '4NPQ',               # code of the target
           loop_assess_methods=assess.DOPE) # assess each loop with DOPE

a.starting_model= 1                 # index of the first model 
a.ending_model  = 1                 # index of the last model

a.loop.starting_model = 1           # First loop model
a.loop.ending_model   = 1           # Last loop model

a.make()                            # do modeling and loop refinement

mdl = model(env, file='4NPQ.BL00010001.pdb')
# Assign new segment names and write out the new model:
mdl.rename_segments(segment_ids=('P','Q','R','S','T'),renumber_residues=[5,5,5,5,5])
s=selection(mdl)
sel=s.only_atom_types('CA')
sel.write(file='4NPQ_4_fixed.pdb', model_format='PDB', no_ter=False)
quit()
