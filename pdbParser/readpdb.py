# See COPYING for license

# We need to read quite a bit of information from the REMARKS lines. So it is better to keep them in a seperate module.
import numpy as np
import urllib.request, urllib.error, urllib.parse
import logging

FILTER_MOLTYPES = ["DNA", "RNA", "5'", "3'", "FAB", "ScFv"]


def getpdb(pdbid, download=False):
    """
    Download or read in PDB files.
    """
    if download:
        try:
            urllib.request.urlretrieve("http://files.rcsb.org/download/%s.pdb" % pdbid, pdbid + ".pdb")
            pdblines = open(pdbid + ".pdb", "r").readlines()
            return pdblines
        except urllib.error.HTTPError as err:
            if err.code == 404:
                raise ValueError("PDB code not found")
            else:
                raise ValueError("Cannot reach PDB database, unknown error")
    else:
        logging.info(f"Opening file {pdbid}")
        with open(pdbid, "r") as pdbfile:
            pdblines = pdbfile.readlines()
        return pdblines


def checkmulti(pdb):
    """
    Check if there are multiple models in the PDB file.
    """
    count = 0
    for line in pdb:
        if line.startswith("MODEL"):
            count += 1
    if count > 1:
        return True
    else:
        return False


def get_compound_info(pdb):
    """
    Filter compound lines
    """
    compndlist = []
    for line in pdb:
        if "COMPND" in line:
            compndlist.append(line[10:-1].strip())
    return compndlist


def get_molid_index(compndlist):
    """
    Return indexes of molecule ids
    """
    molin = []
    for line in compndlist:
        if line.startswith("MOL_ID:"):
            molin.append(compndlist.index(line))
    return molin


def check_moltype_filter(lines):
    """
    Check if mol can be filtered.
    """
    for line in lines:
        if any(s in line for s in FILTER_MOLTYPES):
            return True
    return False


def extract_relevant_molid(molin, compndlist):
    """
    Extract relevant molid, after filtering there must be a single MOLID remaining, otherwise an error is raised.
    """
    if len(molin) > 1:
        passing_mol = []
        mols = {1: compndlist[molin[0] : molin[1]], 2: compndlist[molin[1] :]}
        for mol, lines in mols.items():
            if not check_moltype_filter(lines):
                passing_mol.append(mol)
        if len(passing_mol) != 1:
            raise ValueError("Current version cannot handle this pdb file, both interaction partners are proteins.")
        return mols[passing_mol[0]]
    else:
        return compndlist


def extract_chains(mol_lines):
    """
    Extract chain ids that belongs to the molid.
    """
    for line in mol_lines:
        if "CHAIN:" in line:
            return [i.strip().strip(";") for i in line.split(":")[1].split(",")]
    raise ValueError("Cannot extract chain information.")


def findchains(pdb):
    """
    Find chains that belongs to the complex from COMPND lines.
    """
    compndlist = get_compound_info(pdb)
    molin = get_molid_index(compndlist)
    relevant_compounds = extract_relevant_molid(molin, compndlist)
    return extract_chains(relevant_compounds)


def extract_missing_residues(pdb, chlist):
    """
    Extracts missing residues.
    """
    # Where this function is used again?
    remark = []
    for line in pdb:
        if "REMARK 465" in line:
            remark.append(line)
    r465 = np.genfromtxt(
        remark[7:-1], names=["REMARK", "465", "rname", "ch", "rid"], dtype=["U6", int, "U3", "U1", int]
    )
    filt = r465[np.in1d(r465["ch"], chlist)]
    return filt


def read_atom(pdb):
    """
    Return atom lines. Reading heteroatom lines might present other complications like mistaking calcium atoms for Calpha.
    For this reason, currently we are only reading atom lines.
    """
    atoms = []
    for line in pdb:
        if line[0:6].strip() in ["ATOM"]:
            atoms.append(line)
    return atoms


def coord(atomlines):
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
        if len(atom) == 79:
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
    return coords
