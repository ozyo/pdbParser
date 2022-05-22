from pathlib import Path
import pytest
import numpy as np

from enspdb.readpdb import (
    FILTER_MOLTYPES,
    getpdb,
    check_moltype,
    extract_missing_residues,
    extract_relevant_chains,
    return_proteins_from_file,
)

DATA_DIR = Path(__file__).resolve().parents[0] / "test_data"

TEST_PDB = DATA_DIR / "1BA2.pdb"

struct = return_proteins_from_file(TEST_PDB)


def test_getpdb():
    assert getpdb(Path("1BA2"), download=True, clean=True) == TEST_PDB
    assert getpdb(Path("1BA2"), cwd=DATA_DIR) == TEST_PDB


def test_extract_missing_residues():
    result = []
    assert extract_missing_residues(struct, ["A", "B"]) == result

def test_extract_relevant_chains():
    pass

def test_check_moltype():
    pass
