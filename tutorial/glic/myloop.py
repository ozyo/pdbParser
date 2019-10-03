#Example taken from
#https://salilab.org/modeller/manual/node35.html
from modeller import *
from modeller.automodel import *
class MyLoop(loopmodel):
    # This routine picks the residues to be refined by loop modeling
    # Note that the residue numbering is based on internal MODELLER numbering
    def select_loop_atoms(self):
        # Two residue ranges (both will be refined simultaneously)
        # 180 to 181
        return selection(self.residue_range('1121:D', '1122:D')
                         #self.residue_range('146:A', '146:A'),
                         #self.residue_range('164:A', '179:A')
                         )

    #Uncomment the lines below to apply a secondary structure restrain
    #def special_restraints(self, aln):
        #rsr = self.restraints
        #at = self.atoms
        #rsr.add(secondary_structure.alpha(self.residue_range('37:A', '38:A')))
        #rsr.add(secondary_structure.alpha(self.residue_range('53:A', '56:A')))
        #rsr.add(secondary_structure.alpha(self.residue_range('126:A', '128:A')))
        #rsr.add(secondary_structure.alpha(self.residue_range('130:A', '130:A')))
        #rsr.add(secondary_structure.alpha(self.residue_range('134:A', '134:A')))
