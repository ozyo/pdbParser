# See LICENSE for license

from os import getcwd
from pathlib import Path
from enspdb.utils import ParserError
from typing import List
import urllib.request, urllib.error, urllib.parse
import logging

from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.Polypeptide import PPBuilder

FILTER_MOLTYPES = ["DNA", "RNA", "5'", "3'", "FAB", "ScFv"]


def getpdb(
    pdbid: Path, download: bool = False, clean: bool = False, cwd: Path = Path(getcwd())
) -> List[str]:
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
            pdblines = open(f"{pdbid}", "r").readlines()
            if clean:
                (cwd / pdbid).unlink()
            return pdblines
        except urllib.error.HTTPError as err:
            logging.error("FAIL", exc_info=err)
            raise ParserError
    else:
        if (cwd / pdbid).is_file():
            logging.info(f"Opening file {pdbid}")
            with open(cwd / pdbid, "r") as pdbfile:
                pdblines = pdbfile.readlines()
            return pdblines
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
    last = 0
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
    residues = []
    for pep in polypep:
        for res in pep:
            residues.append(res)
    return residues


def return_proteins_from_file(pdb_path: Path, use_model: int = 0) -> Structure:
    struct = PDBParser().get_structure(pdb_path.stem, pdb_path.as_posix())
    proteins = PPBuilder().build_peptides(struct[0])
    residues = flatten_polypeptide_to_residue(proteins)
    return unfold_entities(residues, "S")[0]
