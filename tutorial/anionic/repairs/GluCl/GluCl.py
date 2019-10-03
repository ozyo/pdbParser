from modeller import *
from modeller.automodel import *

log.verbose()
env = environ()
env.io.hetatm = True
env.io.atom_files_directory = ['./']
class MyModel(automodel):
    def select_atoms(self):
        return selection(self.residue_range('102:A', '106:A'),
                         self.residue_range('441:B', '445:B'),
                         self.residue_range('780:C', '784:C'),
                         self.residue_range('1119:D', '1123:D'),
                         self.residue_range('1458:E', '1462:E'))

a = MyModel (env, 
               alnfile='result.pir_aln',
               knowns='4TNV_assembly2_onlyprotein',
               sequence='4TNV_fix_BA2',
               assess_methods=(assess.DOPE,
                              #soap_protein_od.Scorer(),
                              assess.GA341))

a.starting_model = 1
a.ending_model = 1
a.make()
