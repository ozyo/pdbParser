# See COPYING for license

from typing import List, Union
import numpy as np


def get_chain_order(coord: np.recarray) -> List[str]:
    """
    Return order of the chains in the coordinates.
    """
    _, idx = np.unique(coord["ch"], return_index=True)
    return list(coord["ch"][np.sort(idx)])


def order_chainids(order: List[str], chlist: List[str]) -> List[str]:
    """
    Order chain list.
    """
    return [ch for ch in order if ch in chlist]


def getca_forchains(coord: np.recarray, chlist: List[str], order: bool = False) -> np.recarray:
    """
    Filter CA coordinates that belongs to certain chain ids.
    """
    ca = filter_coordinates(coord, "atname", "CA")
    if order:
        arrays = [ca[np.in1d(ca["ch"], ch)] for ch in chlist]
        ca = arrays[0].__array_wrap__(np.hstack(arrays))
    else:
        # This will always return the order in coordinates.
        ca = ca[np.in1d(ca["ch"], chlist)]
    return ca


def filter_coordinates(coord: np.recarray, column: str, value: Union[str, int, float]) -> np.recarray:
    """
    Filter coordinates based on the value of a column.
    """
    return coord[np.in1d(coord[column], value)]


def clean_altloc(coord: np.recarray, altloc: str):
    """
    Filter altloc lines and clean altloc character.
    """
    filt = coord[(coord["altloc"] == altloc) | (coord["altloc"] == "")]
    filt["altloc"] = ""
    return filt


def clean_icode(coord: np.recarray) -> np.recarray:
    """
    Clean icode character.
    """
    coord["icode"] = ""
    return coord
