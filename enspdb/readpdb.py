# See LICENSE for license

from os import getcwd
from pathlib import Path
from enspdb.utils import ParserError, fold_entities
from typing import List
import urllib.request, urllib.error, urllib.parse
import logging

from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.Polypeptide import PPBuilder

FILTER_MOLTYPES = ["DNA", "RNA", "5'", "3'", "FAB", "ScFv"]
FILTER_TITLES = ["fused", "chimeric", "chimera", "chimaeric"]

def getpdb(
    pdbid: Path, download: bool = False, clean: bool = False, cwd: Path = Path(getcwd())
) -> Structure:
    """
    Download or read in PDB files.
    """
    if pdbid.suffix != ".pdb":
        pdbid = pdbid.with_suffix(".pdb")
    if download:
        logging.info(f"Downloading {pdbid}")
        try:
            urllib.request.urlretrieve(
                f"http://files.rcsb.org/download/{pdbid}", f"{cwd}/{pdbid}"
            )
            pdb_data = return_proteins_from_file(cwd/pdbid)
            if clean:
                (cwd / pdbid).unlink()
            return pdb_data
        except urllib.error.HTTPError as err:
            logging.error("FAIL", exc_info=err)
            raise ParserError
    else:
        if (cwd / pdbid).is_file():
            logging.info(f"Opening file {pdbid}")
            return return_proteins_from_file(cwd / pdbid)
        else:
            try:
                raise FileNotFoundError(f"{pdbid} not found in {cwd.absolute()}")
            except Exception as err:
                logging.error("FAIL", exc_info=err)
                raise ParserError


def check_moltype(compndlist: List[str]) -> bool:
    """
    Check if mol can be filtered.
    """
    for line in compndlist:
        if any(s in line for s in FILTER_MOLTYPES):
            return True
    return False


def extract_relevant_chains(struct: Structure) -> List[str]:
    """
    Extract relevant chains using the compound information from the header.
    After filtering there must be a single MOLID remaining, otherwise an error is raised.
    This limitation is to prevent any errors, the code is yet to be tested with complexes.
    """
    mols = struct.header["compound"]
    passing_mol = dict()

    for mol in range(mols):
        if not check_moltype(mols[str(mol)]["molecule"]):
            passing_mol[str(mol)] = mols[str(mol)]
    if len(passing_mol) != 1:
        try:
            raise ValueError(
                "Current version cannot handle this pdb file, both interaction partners are proteins."
            )
        except Exception as err:
            logging.error("FAIL", exc_info=err)
            raise ParserError
    return list(map(str.upper, list(passing_mol.values())[0]["chain"]))


def extract_missing_residues(struct: Structure, chains: List[str]):
    """
    Extracts missing residues.
    """
    missing = []
    for res in struct.header["missing_residues"]:
        if res["chain"] in chains:
            missing.append(res)
    return missing


def flatten_polypeptide_to_residue(polypep: list) -> list:
    """
    Turn polypeptide entry to residues
    """
    residues = []
    for pep in polypep:
        for res in pep:
            residues.append(res)
    return residues


def return_proteins_from_file(pdb_path: Path, use_model: int = 0) -> Structure:
    """
    An advanced structure loading function that returns only the protein residues.
    """
    struct = PDBParser(QUIET=1).get_structure(pdb_path.stem, pdb_path.as_posix())
    proteins = PPBuilder().build_peptides(struct[use_model])
    residues = flatten_polypeptide_to_residue(proteins)
    try:
        return fold_entities(residues, "S")[0]
    except:
        raise ValueError(f"{pdb_path} has failed")

def load_structure(pdb_path:Path)->Structure:
    """
    A simpler structure loading function
    """
    return PDBParser(QUIET=1).get_structure(pdb_path.stem, pdb_path.as_posix())

def check_pdb_title(title:str)->bool:
    """
    Check if PDB title indicates that this structure is a chimera.
    """
    if any(True for x in FILTER_TITLES if x in title.lower()):
        return True
    else:
        return False