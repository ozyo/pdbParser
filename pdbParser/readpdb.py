# See COPYING for license

# We need to read quite a bit of information from the REMARKS lines. So it is better to keep them in a seperate module.
from os import getcwd
from pathlib import Path
from typing import List
import numpy as np
import urllib.request, urllib.error, urllib.parse
import logging

from numpy.lib.arraysetops import unique

FILTER_MOLTYPES = ["DNA", "RNA", "5'", "3'", "FAB", "ScFv"]


def getpdb(pdbid: Path, download: bool = False, clean: bool = False, cwd: Path = Path(getcwd())) -> List[str]:
    """
    Download or read in PDB files.
    """
    if pdbid.suffix != ".pdb":
        pdbid = pdbid.with_suffix(".pdb")
    if download:
        logging.info(f"Downloading {pdbid}")
        try:
            urllib.request.urlretrieve(f"http://files.rcsb.org/download/{pdbid}", f"{cwd}/{pdbid}")
            pdblines = open(f"{pdbid}", "r").readlines()
            if clean:
                (cwd / pdbid).unlink()
            return pdblines
        except urllib.error.HTTPError as err:
            if err.code == 404:
                raise ValueError("PDB code not found")
            else:
                raise ValueError("Cannot reach PDB database, unknown error")
    else:
        if (cwd / pdbid).is_file():
            logging.info(f"Opening file {pdbid}")
            with open(cwd / pdbid, "r") as pdbfile:
                pdblines = pdbfile.readlines()
            return pdblines
        else:
            raise FileExistsError(f"{pdbid} not found in {cwd.absolute()}")


def check_multimodel(pdblines: List[str]) -> bool:
    """
    Check if there are multiple models in the PDB file.
    """
    count = 0
    for line in pdblines:
        if line.startswith("MODEL"):
            count += 1
    if count > 1:
        return True
    else:
        return False


def get_compound_info(pdblines: List[str]) -> List[str]:
    """
    Filter compound lines
    """
    compndlist = []
    for line in pdblines:
        if "COMPND" in line:
            compndlist.append(line[10:-1].strip())
    return compndlist


def get_molid_index(compndlist: List[str]) -> List[int]:
    """
    Return indexes of molecule ids
    """
    molin = []
    for line in compndlist:
        if line.startswith("MOL_ID:"):
            molin.append(compndlist.index(line))
    return molin


def check_moltype_filter(compndlist: List[str]) -> bool:
    """
    Check if mol can be filtered.
    """
    for line in compndlist:
        if any(s in line for s in FILTER_MOLTYPES):
            return True
    return False


def extract_relevant_molid(molin: List[int], compndlist: List[str]):
    """
    Extract relevant molid, after filtering there must be a single MOLID remaining, otherwise an error is raised.
    """
    if len(molin) > 1:
        passing_mol = []
        for i in range(len(molin) - 1):
            lines = compndlist[molin[i] : molin[i + 1]]
            if not check_moltype_filter(lines):
                passing_mol.append(i)
        # Check the last element
        if not check_moltype_filter(compndlist[len(molin) :]):
            passing_mol.append(i + 1)
        if len(passing_mol) != 1:
            raise ValueError("Current version cannot handle this pdb file, both interaction partners are proteins.")
        if passing_mol[0] == len(molin) - 1:
            return compndlist[molin[passing_mol[0]] :]
        return compndlist[molin[passing_mol[0]] : molin[passing_mol[0] + 1]]
    else:
        return compndlist


def extract_chains(mol_lines: List[str]):
    """
    Extract chain ids that belongs to the molid.
    """
    for line in mol_lines:
        if "CHAIN:" in line:
            return [i.strip().strip(";") for i in line.split(":")[1].split(",")]
    raise ValueError("Cannot extract chain information.")


def find_chains(pdblines: List[str]) -> List[str]:
    """
    Find chains that belongs to the complex from COMPND lines.
    """
    compndlist = get_compound_info(pdblines)
    molin = get_molid_index(compndlist)
    relevant_compounds = extract_relevant_molid(molin, compndlist)
    return extract_chains(relevant_compounds)


# TODO: It's a little dangerous to assume that missing residues start exactly at 7th element.
# This will break for any unconventional PDB file with missing missing residue lines.
def extract_missing_residues(pdblines: List[str], chlist: List[str]) -> np.recarray:
    """
    Extracts missing residues.
    """
    # Where this function is used again?
    remark = []
    for line in pdblines:
        if "REMARK 465" in line:
            remark.append(line)
    r465 = np.genfromtxt(remark[7:], names=["REMARK", "465", "rname", "ch", "rid"], dtype=["U6", int, "U3", "U1", int])
    filt = r465[np.in1d(r465["ch"], chlist)]
    return filt


def read_atom(pdblines: List[str]) -> List[str]:
    """
    Return atom lines. Reading heteroatom lines might present other complications like mistaking calcium atoms for Calpha.
    For this reason, currently we are only reading atom lines.
    """
    atoms = []
    for line in pdblines:
        if line[0:6].strip() in ["ATOM"]:
            atoms.append(line)
    return atoms


def coord(atomlines: List[str]) -> np.recarray:
    """
    Create a coordinate array from atomlines.
    """
    coords = []
    for atom in atomlines:
        atnr = int(atom[6:11].strip())
        atname = atom[12:16].strip()
        altloc = atom[16].strip()
        resname = atom[17:20].strip()
        ch = atom[21]
        resnr = int(atom[22:26])
        icode = atom[26].strip()
        x = float(atom[30:38].strip())
        y = float(atom[38:46].strip())
        z = float(atom[46:54].strip())
        if len(atom) >= 79:
            occu = float(atom[54:61].strip())
            tfact = float(atom[61:66].strip())
            element = atom[76:78].strip()
            charge = atom[78:80].strip()
        else:
            occu = float()
            tfact = float()
            element = ""
            charge = ""
        coords.append((atnr, atname, altloc, resname, ch, resnr, icode, x, y, z, occu, tfact, element, charge))
    coords = np.array(coords, dtype=("i,U4,U4,U4,U4,i,U4,f,f,f,f,f,U4,U4"))
    coords.dtype.names = (
        "atnr",
        "atname",
        "altloc",
        "resname",
        "ch",
        "resnr",
        "icode",
        "x",
        "y",
        "z",
        "occu",
        "tfact",
        "element",
        "charge",
    )
    unique_indices = np.sort(
        np.unique(coords[["atname", "altloc", "resname", "ch", "resnr"]], return_index=True, axis=0)[1]
    )
    return coords[unique_indices]
