from pathlib import Path
import pytest
import numpy as np

from pdbParser.readpdb import (
    FILTER_MOLTYPES,
    getpdb,
    check_multimodel,
    get_compound_info,
    get_molid_index,
    check_moltype_filter,
    extract_relevant_molid,
    extract_chains,
    find_chains,
    extract_missing_residues,
    read_atom,
    coord,
)

DATA_DIR = Path(__file__).resolve().parents[0] / "test_data"

TEST_PDB = open(f"{DATA_DIR}/1BA2.pdb").readlines()

COMPND_SINGLE = [
    "MOL_ID: 1;",
    "MOLECULE: D-RIBOSE-BINDING PROTEIN;",
    "CHAIN: A, B;",
    "ENGINEERED: YES;",
    "MUTATION: YES",
]

COMPND_MULTI_TEMPLATE = [
    "MOL_ID: 1;",
    "MOLECULE: D-RIBOSE-BINDING PROTEIN;",
    "CHAIN: A, B;",
    "ENGINEERED: YES;",
    "MUTATION: YES",
    "MOL_ID: 2;",
    "MOLECULE: %s;",
    "CHAIN: C,D;",
]

MISSING_RESIDUES = [
    "REMARK 465 MISSING RESIDUES",
    "REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE",
    "REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN",
    "REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)",
    "REMARK 465",
    "REMARK 465   M RES C SSSEQI",
    "REMARK 465     ILE A   103",
    "REMARK 465     ASP A   104",
    "REMARK 465     LYS B   105",
    "REMARK 465     HIS B   341",
    "REMARK 465     HIS C   342",
]

LAST_ATOMS_TEST_PDB = [
    "ATOM   4005  O   GLN B 271     -20.348   4.657   4.099  1.00 31.34           O  \n",
    "ATOM   4006  CB  GLN B 271     -20.068   3.067   6.959  1.00 30.59           C  \n",
    "ATOM   4007  CG  GLN B 271     -20.350   1.718   7.641  1.00 30.95           C  \n",
    "ATOM   4008  CD  GLN B 271     -18.966   1.143   7.945  1.00 32.20           C  \n",
    "ATOM   4009  OE1 GLN B 271     -18.243   1.762   8.735  1.00 32.58           O  \n",
    "ATOM   4010  NE2 GLN B 271     -18.646   0.036   7.286  1.00 32.60           N  \n",
    "ATOM   4011  OXT GLN B 271     -21.819   3.067   4.468  1.00 31.48           O  \n",
]

LAST_ATOMS_TEST_PDB_SHORT = [i[:60] for i in LAST_ATOMS_TEST_PDB]

DUPLICATE_TEST_PDB = [
    "ATOM   4005  O   GLN B 271     -20.348   4.657   4.099  1.00 31.34           O  \n",
    "ATOM   4006  CB  GLN B 271     -20.068   3.067   6.959  1.00 30.59           C  \n",
    "ATOM   4007  CG  GLN B 271     -20.350   1.718   7.641  1.00 30.95           C  \n",
    "ATOM   4008  CD  GLN B 271     -18.966   1.143   7.945  1.00 32.20           C  \n",
    "ATOM   4005  O   GLN B 271     -20.348   4.657   4.099  1.00 31.34           O  \n",
    "ATOM   4009  OE1 GLN B 271     -18.243   1.762   8.735  1.00 32.58           O  \n",
    "ATOM   4010  NE2 GLN B 271     -18.646   0.036   7.286  1.00 32.60           N  \n",
    "ATOM   4005  O   GLN B 271     -20.348   4.657   4.099  1.00 31.34           O  \n",
    "ATOM   4011  OXT GLN B 271     -21.819   3.067   4.468  1.00 31.48           O  \n",
    "ATOM   4010  NE2 GLN B 271     -18.646   0.036   7.286  1.00 32.60           N  \n",
    "ATOM   4007  CG  GLN B 271     -20.350   1.718   7.641  1.00 30.95           C  \n",
]

COORDS = np.array(
    [
        (4005, "O", "", "GLN", "B", 271, "", -20.348, 4.657, 4.099, 1.00, 31.34, "O", ""),
        (4006, "CB", "", "GLN", "B", 271, "", -20.068, 3.067, 6.959, 1.00, 30.59, "C", ""),
        (4007, "CG", "", "GLN", "B", 271, "", -20.350, 1.718, 7.641, 1.00, 30.95, "C", ""),
        (4008, "CD", "", "GLN", "B", 271, "", -18.966, 1.143, 7.945, 1.00, 32.20, "C", ""),
        (4009, "OE1", "", "GLN", "B", 271, "", -18.243, 1.762, 8.735, 1.00, 32.58, "O", ""),
        (4010, "NE2", "", "GLN", "B", 271, "", -18.646, 0.036, 7.286, 1.00, 32.60, "N", ""),
        (4011, "OXT", "", "GLN", "B", 271, "", -21.819, 3.067, 4.468, 1.00, 31.48, "O", ""),
    ],
    dtype=("i,U4,U4,U4,U4,i,U4,f,f,f,f,f,U4,U4"),
)

COORDS_SHORT = np.array(
    [
        (4005, "O", "", "GLN", "B", 271, "", -20.348, 4.657, 4.099, float(), float(), "", ""),
        (4006, "CB", "", "GLN", "B", 271, "", -20.068, 3.067, 6.959, float(), float(), "", ""),
        (4007, "CG", "", "GLN", "B", 271, "", -20.350, 1.718, 7.641, float(), float(), "", ""),
        (4008, "CD", "", "GLN", "B", 271, "", -18.966, 1.143, 7.945, float(), float(), "", ""),
        (4009, "OE1", "", "GLN", "B", 271, "", -18.243, 1.762, 8.735, float(), float(), "", ""),
        (4010, "NE2", "", "GLN", "B", 271, "", -18.646, 0.036, 7.286, float(), float(), "", ""),
        (4011, "OXT", "", "GLN", "B", 271, "", -21.819, 3.067, 4.468, float(), float(), "", ""),
    ],
    dtype=("i,U4,U4,U4,U4,i,U4,f,f,f,f,f,U4,U4"),
)

COORDS.dtype.names = (
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

COORDS_SHORT.dtype.names = (
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


def test_getpdb():
    assert getpdb(Path("1BA2"), download=True, clean=True) == TEST_PDB
    assert getpdb(Path("1BA2"), cwd=DATA_DIR) == TEST_PDB


@pytest.mark.parametrize(
    "pdblines,result", [(TEST_PDB, False), (["MODEL", "MODEL"], True), (["MODEL"] + TEST_PDB, False)]
)
def test_check_multimodel(pdblines, result):
    assert check_multimodel(pdblines) == result


def test_get_compound_info():
    assert get_compound_info(TEST_PDB) == COMPND_SINGLE


def test_get_molid_index():
    assert get_molid_index(COMPND_MULTI_TEMPLATE) == [0, 5]


def yield_compnd():
    for mol_type in FILTER_MOLTYPES:
        tmp = []
        for i in COMPND_MULTI_TEMPLATE:
            if "%" in i:
                tmp.append(i % mol_type)
            else:
                tmp.append(i)
        return tmp


@pytest.mark.parametrize("lines,result", [(yield_compnd(), True), (COMPND_SINGLE, False)])
def test_check_moltype_filter(lines, result):
    assert check_moltype_filter(lines) == result


# TODO: Add tests where this function will raise an error.
@pytest.mark.parametrize("lines,mol_index,result", [(yield_compnd(), [0, 5], COMPND_SINGLE)])
def test_extract_relevant_molid(lines, mol_index, result):
    assert extract_relevant_molid(mol_index, lines) == result


def test_extract_chains():
    assert extract_chains(COMPND_SINGLE) == ["A", "B"]


# TODO: Add tests where sub fuctions will raise an error.
def test_find_chains():
    assert find_chains(TEST_PDB) == ["A", "B"]


@pytest.mark.parametrize(
    "lines,chids,result",
    [
        (TEST_PDB, ["A", "B"], []),
        ((TEST_PDB + MISSING_RESIDUES), ["A", "B"], MISSING_RESIDUES[7:-1]),
    ],
)
def test_extract_missing_residues(lines, chids, result):
    result = np.genfromtxt(result, names=["REMARK", "465", "rname", "ch", "rid"], dtype=["U6", int, "U3", "U1", int])
    np.testing.assert_array_equal(extract_missing_residues(lines, chids), result)


def test_read_atom():
    assert read_atom(TEST_PDB)[-7:] == LAST_ATOMS_TEST_PDB


def test_coord():
    np.testing.assert_array_equal(coord(DUPLICATE_TEST_PDB), COORDS, verbose=True)
    np.testing.assert_array_equal(coord(LAST_ATOMS_TEST_PDB), COORDS, verbose=True)
    np.testing.assert_array_equal(coord(LAST_ATOMS_TEST_PDB_SHORT), COORDS_SHORT, verbose=True)
