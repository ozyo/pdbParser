from copy import copy
import itertools
from pathlib import Path
from typing import Union
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Entity import Entity
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.PDBIO import Select, PDBIO
from Bio.PDB.PDBExceptions import PDBException

ENTITIY_LEVELS=["A", "R", "C", "M", "S"]

class ParserError(Exception):
    """
    Raised when a known error is caught.
    """

    pass


class ChainSelect(Select):
    def __init__(self, chains):
        self.chains = chains

    def accept_chain(self, chain):
        return chain.id in self.chains


class CAlphaSelect(Select):
    def accept_atom(self, atom):
        return atom.id == "CA"

class SelectResidues(Select):
    def __init__(self, chain_res_ids):
        self.chain_res_ids = chain_res_ids

    def accept_chain(self, chain):
        return chain.id in self.chain_res_ids.keys()

    def accept_residue(self, residue):
        return residue.id[1] in self.chain_res_ids[residue.get_full_id()[2]]


def get_selected(struct: Union[Structure,Model,Chain,Residue], selection: Select, return_level:str="S", fold:bool=True) -> list:
    selected = []

    for atom in struct.get_atoms():
        residue = atom.get_parent()
        chain = residue.get_parent()
        model = chain.get_parent()
        if (
            selection.accept_atom(atom)
            and selection.accept_residue(residue)
            and selection.accept_chain(chain)
            and selection.accept_model(model)
        ):
            selected.append(atom)
    if len(selected) == 0:
        raise ValueError(f"Selection did not return any atoms for {struct.get_full_id()[0]}")
    if fold:
        return fold_entities(selected, return_level)
    return unfold_entities(selected,return_level)


def write_pdb(struct: Structure, outname: Path):
    writer = PDBIO()
    writer.set_structure(struct)
    writer.save(outname.as_posix())


def fold_entities(entity_list, target_level):
    """
    Return the desired target level by removing every other entity that are not
    in the entity list.

    This method doesn't create a new structure entry, it instead modifies
    a copy of the entities' parents.
    """
    if target_level not in ENTITIY_LEVELS:
        raise PDBException("%s: Not an entity level." % target_level)
    if entity_list == []:
        return []
    if isinstance(entity_list, (Entity, Atom)):
        entity_list = [entity_list]

    level = entity_list[0].get_level()
    if not all(entity.get_level() == level for entity in entity_list):
        raise PDBException("Entity list is not homogeneous.")

    target_index = ENTITIY_LEVELS.index(target_level)
    level_index = ENTITIY_LEVELS.index(level)

    if level_index == target_index:  # already right level
        return entity_list
    if level_index != 0:
        full_ids = [atom.get_full_id() for ent in entity_list for atom in list(ent.get_atoms())]
    else:
        full_ids = [atom.get_full_id() for atom in entity_list]
    for i in range(level_index, target_index):
        if i != target_index-1:
            entity_list = list(dict.fromkeys(entity.get_parent() for entity in entity_list))
        else:
            # parent is detached when you copy with bioppython so we do it at the last entry
            entity_list = list(dict.fromkeys(entity.get_parent().copy() for entity in entity_list))
    for parent in entity_list:
        for atom in list(parent.get_atoms()):
            if atom.get_full_id() not in full_ids:
                atom.get_parent().detach_child(atom.id)

    for parent in entity_list:
        for i in range(1,target_index):
            if i == 1:
                for res in list(parent.get_residues()):
                    if len(list(res)) == 0:
                        res.get_parent().detach_child(res.id)
            if i == 2:
                for ch in list(parent.get_chains()):
                    if len(list(ch)) == 0:
                        ch.get_parent().detach_child(ch.id)
            if i == 3:
                for model in list(parent.get_models()):
                    if len(list(model)) == 0:
                        model.get_parent().detach_child(model.id)
    return list(entity_list)