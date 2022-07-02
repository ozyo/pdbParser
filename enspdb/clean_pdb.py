# See LICENSE for license

from Bio.PDB.Structure import Structure
from Bio.PDB.Selection import unfold_entities
from enspdb.utils import CAlphaSelect, ChainSelect, CleanDisorderedChilds, get_selected
from enspdb.modifypdb import change_icode

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


def clean_altloc(struct: Structure, altloc: str) -> Structure:
    """
    Filter altloc lines and clean altloc character.
    """
    return unfold_entities(get_selected(struct, CleanDisorderedChilds(altloc)), "S")[0]


def clean_icode(struct: Structure) -> None:
    """
    Clean icode character. If this function is called before renumbering it might raise
    an error due to overlapping residue numbers
    """
    for res in struct.get_residues():
        change_icode(res,"")
