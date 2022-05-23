from pathlib import Path
from Bio.PDB.Structure import Structure
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.PDBIO import Select, PDBIO


class ParserError(Exception):
    """
    Raised when a known error is caught.
    """

    pass


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


class SelectResidues(Select):
    def __init__(self, chain, res_ids):
        self.chain = chain
        self.res_ids = res_ids

    def accept_residue(self, residue):
        return all(
            [residue.parent.id == self.chain, residue.full_id()[1] in self.res_ids]
        )


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


def write_pdb(struct: Structure, outname: Path):
    writer = PDBIO()
    writer.set_structure(struct)
    writer.save(outname.as_posix())
