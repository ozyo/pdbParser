# See COPYING for license

import numpy as np


def get_chain_order(coord):
    """
    Return order of the chains in the coordinates.
    """
    _, idx = np.unique(coord["ch"], return_index=True)
    return coord["ch"][np.sort(idx)]


def order_chainids(order, chlist):
    """
    Order chain list.
    """
    return [ch for ch in order if ch in chlist]


def getca_forchains(coord, chlist, order=False):
    """
    Filter CA coordinates that belongs to certain chain ids.
    """
    ca = filter_coordinates(coord, "atname", "CA")
    if order:
        ca[np.in1d(ca["ch"], chlist)]
    else:
        ordch = order_chainids(get_chain_order(ca), chlist)
        ca = ca[np.in1d(ca["ch"], ordch)]
    return ca


def filter_coordinates(coord, column, value):
    """
    Filter coordinates based on the value of a column.
    """
    return coord[np.in1d(coord[column], value)]


def clean_altloc(coord, altloc):
    """
    Filter altloc lines and clean altloc character.
    """
    filt = coord[(coord["altloc"] == altloc) | (coord["altloc"] == "")]
    filt["altloc"] = ""
    return filt


def clean_icode(coord):
    """
    Clean icode character.
    """
    coord["icode"] = ""
