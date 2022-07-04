# See LICENSE for license

from Bio.PDB.Structure import Structure
from Bio.PDB.Selection import unfold_entities
from enspdb.utils import CAlphaSelect, ChainSelect, fold_entities, get_selected
from enspdb import utils as u
from enspdb.modifypdb import change_icode

def check_chains_exists(struct: Structure, chains:list):
    for chain in struct.get_chains():
        if chain.id not in chains:
            raise ValueError(
                f"{struct.id} chain {chain.id} is not found in {''.join(chains)}."
            )
    return True


def getca_forchains(
    struct: Structure, chains: list, order: bool = False
) -> Structure:
    """
    Filter CA coordinates that belongs to certain chain ids.
    """
    struct_chains = get_selected(struct, ChainSelect(chains))[0]
    check_chains_exists(struct_chains, chains)
    struct_ca = get_selected(struct_chains, CAlphaSelect())[0]
    if order:
        return fold_entities(
            [
                chain
                for ord_chain in chains
                for chain in struct_ca.get_chains()
                if chain.id == ord_chain
            ],
            "S",
        )[0]
    return struct_ca


def clean_altloc(struct: Structure, altloc: str) -> Structure:
    """
    Filter altloc lines and clean altloc character.
    """
    atoms=[]
    for atom in struct.get_atoms():
        if atom.disordered_flag == 1 and len(atom.child_dict) > 1:
            for child_label in list(atom.disordered_get_id_list()):
                if child_label == altloc:
                    atoms.append(atom.child_dict[child_label])
                else:
                    atom.disordered_remove(child_label)
        elif atom.disordered_flag == 1 and len(atom.child_dict) == 1:
            atoms.append(list(atom.child_dict.values())[0])
        else:
            atoms.append(atom)
    return u.fold_entities(atoms, "S")[0]


def clean_icode(struct: Structure) -> None:
    """
    Clean icode character. If this function is called before renumbering it might raise
    an error due to overlapping residue numbers
    """
    for res in struct.get_residues():
        change_icode(res,"")
