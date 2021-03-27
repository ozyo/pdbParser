# See COPYING for license

from pdbParser.readpdb import read_atom, coord
from pdbParser.clean_pdb import clean_icode, clean_altloc, getca_forchains, get_chain_order
import logging

TITLE_FILTERS = ["fused", "chimeric", "chimera", "chimaeric"]


def pdb_title(pdb):
    """
    Check if PDB title indicates that this structure is a chimera.
    """
    titles = []
    for line in pdb:
        if "TITLE" in line:
            titles.append(line)
    for i in titles:
        if any(True for x in TITLE_FILTERS if x in i.lower()):
            return True
        else:
            return False


def parse_ca(pdb, chains, altloc="A", remove_icode=True, order_chainids=False):
    """
    Parse the CA coordinates of the chain ids and clean alloc. Optionally remove icode characters.
    The order of the coordinates will be the same as original coordinates, unless order_chainids requested.
    """
    logging.info("Retriving CA coordinates")
    atomlines = read_atom(pdb)
    coords = coord(atomlines)
    coords = clean_altloc(coords, altloc)
    coords = getca_forchains(coords, chains, order_chainids)
    if remove_icode:
        clean_icode(coords)
    return coords


def parse_chlist(pdb):
    """
    Parse chain information of the coordinates.
    """
    logging.info("Parsing the chain information only")
    atomlines = read_atom(pdb)
    coords = coord(atomlines)
    return get_chain_order(coords)
