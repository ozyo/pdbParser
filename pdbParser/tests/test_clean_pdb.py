import numpy as np

from pdbParser.clean_pdb import (
    get_chain_order,
    order_chainids,
    getca_forchains,
    filter_coordinates,
    clean_altloc,
    clean_icode,
)

COORDS = np.array(
    [
        (4005, "CA", "A", "GLN", "A", 271, "", -20.348, 4.657, 4.099, 1.00, 31.34, "C", ""),
        (4005, "CA", "B", "GLN", "A", 271, "", -20.348, 4.657, 4.099, 1.00, 31.34, "C", ""),
        (4006, "O", "", "GLN", "A", 271, "", -20.068, 3.067, 6.959, 1.00, 30.59, "O", ""),
        (4007, "CA", "", "GLN", "B", 272, "", -20.350, 1.718, 7.641, 1.00, 30.95, "C", ""),
        (4008, "O", "", "GLN", "B", 272, "", -18.966, 1.143, 7.945, 1.00, 32.20, "O", ""),
        (4009, "CA", "", "GLN", "C", 271, "", -18.243, 1.762, 8.735, 1.00, 32.58, "C", ""),
        (4009, "CA", "", "GLN", "C", 271, "A", -18.243, 1.762, 8.735, 1.00, 32.58, "C", ""),
        (4010, "O", "", "GLN", "C", 271, "", -18.646, 0.036, 7.286, 1.00, 32.60, "O", ""),
        (4011, "N", "", "GLN", "C", 271, "", -21.819, 3.067, 4.468, 1.00, 31.48, "N", ""),
    ],
    dtype=("i,U4,U4,U4,U4,i,U4,f,f,f,f,f,U4,U4"),
)

COORDS_CA = np.array(
    [
        (4005, "CA", "A", "GLN", "A", 271, "", -20.348, 4.657, 4.099, 1.00, 31.34, "C", ""),
        (4005, "CA", "B", "GLN", "A", 271, "", -20.348, 4.657, 4.099, 1.00, 31.34, "C", ""),
        (4007, "CA", "", "GLN", "B", 272, "", -20.350, 1.718, 7.641, 1.00, 30.95, "C", ""),
    ],
    dtype=("i,U4,U4,U4,U4,i,U4,f,f,f,f,f,U4,U4"),
)


COORDS_CA_ORDER = np.array(
    [
        (4007, "CA", "", "GLN", "B", 272, "", -20.350, 1.718, 7.641, 1.00, 30.95, "C", ""),
        (4005, "CA", "A", "GLN", "A", 271, "", -20.348, 4.657, 4.099, 1.00, 31.34, "C", ""),
        (4005, "CA", "B", "GLN", "A", 271, "", -20.348, 4.657, 4.099, 1.00, 31.34, "C", ""),
    ],
    dtype=("i,U4,U4,U4,U4,i,U4,f,f,f,f,f,U4,U4"),
)

COORDS_NOALTLOC = np.array(
    [
        (4005, "CA", "", "GLN", "A", 271, "", -20.348, 4.657, 4.099, 1.00, 31.34, "C", ""),
        (4006, "O", "", "GLN", "A", 271, "", -20.068, 3.067, 6.959, 1.00, 30.59, "O", ""),
        (4007, "CA", "", "GLN", "B", 272, "", -20.350, 1.718, 7.641, 1.00, 30.95, "C", ""),
        (4008, "O", "", "GLN", "B", 272, "", -18.966, 1.143, 7.945, 1.00, 32.20, "O", ""),
        (4009, "CA", "", "GLN", "C", 271, "", -18.243, 1.762, 8.735, 1.00, 32.58, "C", ""),
        (4009, "CA", "", "GLN", "C", 271, "A", -18.243, 1.762, 8.735, 1.00, 32.58, "C", ""),
        (4010, "O", "", "GLN", "C", 271, "", -18.646, 0.036, 7.286, 1.00, 32.60, "O", ""),
        (4011, "N", "", "GLN", "C", 271, "", -21.819, 3.067, 4.468, 1.00, 31.48, "N", ""),
    ],
    dtype=("i,U4,U4,U4,U4,i,U4,f,f,f,f,f,U4,U4"),
)

COORDS_NOICODE = np.array(
    [
        (4005, "CA", "A", "GLN", "A", 271, "", -20.348, 4.657, 4.099, 1.00, 31.34, "C", ""),
        (4005, "CA", "B", "GLN", "A", 271, "", -20.348, 4.657, 4.099, 1.00, 31.34, "C", ""),
        (4006, "O", "", "GLN", "A", 271, "", -20.068, 3.067, 6.959, 1.00, 30.59, "O", ""),
        (4007, "CA", "", "GLN", "B", 272, "", -20.350, 1.718, 7.641, 1.00, 30.95, "C", ""),
        (4008, "O", "", "GLN", "B", 272, "", -18.966, 1.143, 7.945, 1.00, 32.20, "O", ""),
        (4009, "CA", "", "GLN", "C", 271, "", -18.243, 1.762, 8.735, 1.00, 32.58, "C", ""),
        (4009, "CA", "", "GLN", "C", 271, "", -18.243, 1.762, 8.735, 1.00, 32.58, "C", ""),
        (4010, "O", "", "GLN", "C", 271, "", -18.646, 0.036, 7.286, 1.00, 32.60, "O", ""),
        (4011, "N", "", "GLN", "C", 271, "", -21.819, 3.067, 4.468, 1.00, 31.48, "N", ""),
    ],
    dtype=("i,U4,U4,U4,U4,i,U4,f,f,f,f,f,U4,U4"),
)

RECARRAY_COLS = [
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
]
COORDS.dtype.names = RECARRAY_COLS
COORDS_CA.dtype.names = RECARRAY_COLS
COORDS_CA_ORDER.dtype.names = RECARRAY_COLS
COORDS_NOALTLOC.dtype.names = RECARRAY_COLS
COORDS_NOICODE.dtype.names = RECARRAY_COLS


def test_get_chain_order():
    get_chain_order(COORDS) == ["A", "B", "C"]


def test_order_chainids():
    assert order_chainids(["C", "B", "A", "D"], ["A", "B", "C", "D"]) == ["C", "B", "A", "D"]
    assert order_chainids(["C", "B", "A", "D", "E", "Z"], ["A", "B", "C", "D"]) == ["C", "B", "A", "D"]


def test_getca_forchains():
    np.testing.assert_array_equal(getca_forchains(COORDS, ["B", "A"], order=False), COORDS_CA)
    np.testing.assert_array_equal(getca_forchains(COORDS, ["B", "A"], order=True), COORDS_CA_ORDER)


def test_filter_coordinates():
    np.testing.assert_array_equal(filter_coordinates(COORDS[:5], "atname", "CA"), COORDS_CA)
    np.testing.assert_array_equal(filter_coordinates(COORDS[:5], "element", "C"), COORDS_CA)


def test_clean_altloc():
    np.testing.assert_array_equal(clean_altloc(COORDS, "B"), COORDS_NOALTLOC)


def test_clean_icode():
    np.testing.assert_array_equal(clean_icode(COORDS), COORDS_NOICODE)
