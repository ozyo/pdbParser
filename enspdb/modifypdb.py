import itertools
from typing import Union
from biopython.Bio.PDB.Structure import Structure
from biopython.Bio.PDB.Model import Model
from biopython.Bio.PDB.Chain import Chain
from biopython.Bio.PDB.Residue import Residue

def change_residue_number(res:Residue,res_nr:int):
    """
    Change residue number
    """
    hetero,cur_res_nr,icode= res.id
    res.id=(hetero,res_nr,icode)

def renumber_polychains(pdb:Union[Structure,Model,Chain]):
    """
    Renumber residues of an entity continously.
    """
    if pdb.get_level() not in ["S","M","C"]:
        raise ValueError("This function is for renumbering the S,M or C level entities.")
    new_nr = 1
    prev_model=0
    prev_chain="A"
    for res in pdb.get_residues():
        _,model,chain,_ = res.get_full_id()
        if model != prev_model or chain != prev_chain:
            new_nr=1
        change_residue_number(res,new_nr)
        new_nr=+1