# See COPYING for license

from pdbParser import readpdb as rp
from pdbParser import clean_pdb as cp
from importlib import reload
reload(rp)
reload(cp)
import logging


def pdb_title(pdb):
    titles = []
    for line in pdb:
        if "TITLE" in line:
            titles.append(line)
    for i in titles:
        if any(True for x in ["fused", "chimeric", "chimera", "chimaeric"] if x in i.lower()):
            return True
        else:
            return False


def parse_ca(pdb, chains, altloc="A", remove_icode=False, order_chainids=False):
    """
    Parse the CA coordinates of the chain ids and clean alloc. Optionally remove icode characters.
    The order of the coordinates will be the same as original coordinates, unless order_chainids requested.
    """
    logging.info("Retriving CA coordinates")
    atomlines = rp.read_atom(pdb)
    coords = rp.coord(atomlines)
    coords = cp.clean_altloc(coords, altloc)
    coords = cp.getca_forchains(coords, chains, order_chainids)
    if remove_icode:
        remove_icode(coords)
    return coords


def parse_chlist(pdb):
    """
    Parse chain information of the coordinates.
    """
    logging.info("Parsing the chain information only")
    atomlines = rp.read_atom(pdb)
    coords = rp.coord(atomlines)
    return cp.get_chain_order(coords)
