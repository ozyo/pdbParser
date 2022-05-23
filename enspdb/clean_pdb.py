# See LICENSE for license

from typing import List, Union
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.PDBIO import Select


def get_chain_order(coord: np.recarray) -> List[str]:
    """
    Return order of the chains in the coordinates.
    """
    _, idx = np.unique(coord["ch"], return_index=True)
    return list(coord["ch"][np.sort(idx)])


def order_chainids(order: List[str], chlist: List[str]) -> List[str]:
    """
    Order chain list.
    """
    return [ch for ch in order if ch in chlist]


class ChainSelect(Select):
    def __init__(self, chains):
        self.chains = chains

    def accept_chain(self, chain):
        return chain in self.chains


class CAlphaSelect(Select):
    def accept_atom(self, atom):
        return atom.id == "CA"


class NotDisordered(Select):
    def __init__(self, keep_altloc):
        self.altloc = keep_altloc

    def accept_atom(self, atom):
        if not atom.is_disordered() or atom.get_altloc() == self.altloc:
            atom.set_altloc(" ")
            return True
        return False


def get_selected(struct: Structure, selection: Select) -> Structure:
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

    return unfold_entities(selected, "S")[0]


def check_chains_exists(struct: Structure, chains: list[str]):
    for chain in struct.get_chains():
        if chain not in chains:
            raise ValueError(
                f"{struct.id} chain {chain} is not found in {''.join(chains)}."
            )
    return True


def getca_forchains(
    struct: Structure, chains: list[str], order: bool = False
) -> Structure:
    """
    Filter CA coordinates that belongs to certain chain ids.
    """
    struct_chains = get_selected(struct, ChainSelect(chains))
    check_chains_exists(struct_chains, chains)
    struct_ca = get_selected(struct_chains, CAlphaSelect()).copy()
    if order:
        return unfold_entities(
            [
                chain
                for ord_chain in chains
                for chain in struct_ca.get_chains()
                if chain.id == ord_chain
            ],
            "S",
        )[0]
    return struct_ca


def clean_altloc(struct: Structure, altloc: str)->Structure:
    """
    Filter altloc lines and clean altloc character.
    """
    return unfold_entities(get_selected(struct, NotDisordered(altloc)),"S")[0]


def clean_icode(struct: Structure)->None:
    """
    Clean icode character.
    """
    for atom in struct.get_atoms():
        # TODO: Implement support for cleaning icode in BioPython
        pass
