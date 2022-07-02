import itertools
from typing import Union
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

def change_residue_number(res:Residue,res_nr:int):
    """
    Change residue number
    """
    hetero,_,icode= res.id
    res.id=(hetero,res_nr,icode)

def change_icode(res:Residue,icode:str):
    """
    Change icode character
    """
    hetero,res_nr,_= res.id
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