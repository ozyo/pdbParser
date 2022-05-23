import os
from pathlib import Path
import numpy as np
from copy import deepcopy

from Bio.PDB.PDBIO import PDBIO
from enspdb.clean_pdb import getca_forchains


def com(xyz, weight, tweight):
    """
    Calculate center-of-mass
    """
    coms = np.sum(xyz, axis=0)
    coms = (coms * weight) / tweight
    return coms


def inertia(xyz, com, weight):
    """
    Calculate intertia vectors
    """
    xyz = deepcopy(xyz) - com
    Ixx = np.sum(weight * (xyz[:, 1] ** 2 + xyz[:, 2] ** 2))
    Iyy = np.sum(weight * (xyz[:, 0] ** 2 + xyz[:, 2] ** 2))
    Izz = np.sum(weight * (xyz[:, 0] ** 2 + xyz[:, 1] ** 2))
    Ixy = np.sum(-(weight * (xyz[:, 0] * xyz[:, 1])))
    Ixz = np.sum(-(weight * (xyz[:, 0] * xyz[:, 2])))
    Iyz = np.sum(-(weight * (xyz[:, 1] * xyz[:, 2])))
    return np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])


def get_rot_mat(a_vec, b_vec):
    """
    Calculate roration matrix
    """
    # where this function is used?
    a_vec = a_vec / np.linalg.norm(a_vec)
    b_vec = b_vec / np.linalg.norm(b_vec)
    cross = np.cross(a_vec, b_vec)
    ab_angle = np.arccos(np.dot(a_vec, b_vec))

    vx = np.array([[0, -cross[2], cross[1]], [cross[2], 0, -cross[0]], [-cross[1], cross[0], 0]])
    R = np.identity(3) * np.cos(ab_angle) + (1 - np.cos(ab_angle)) * np.outer(cross, cross) + np.sin(ab_angle) * vx

    return R


def get_rotation_matrix(i_v, unit=[0, 0, 1]):
    """
    Calculate rotation matrix
    """
    # https://stackoverflow.com/questions/43507491/imprecision-with-rotation-matrix-to-align-a-vector-to-an-axis
    # From http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q38
    if unit is None:
        unit = [1.0, 0.0, 0.0]
    # Normalize vector length
    i_v /= np.linalg.norm(i_v)

    # Get axis
    uvw = np.cross(unit, i_v)

    # compute trig values - no need to go through arccos and back
    rcos = np.dot(i_v, unit)
    rsin = np.linalg.norm(uvw)

    # normalize and unpack axis
    if not np.isclose(rsin, 0):
        uvw /= rsin
    u, v, w = uvw

    # Compute rotation matrix - re-expressed to show structure
    return (
        rcos * np.eye(3)
        + rsin * np.array([[0, -w, v], [w, 0, -u], [-v, u, 0]])
        + (1.0 - rcos) * uvw[:, None] * uvw[None, :]
    )


def align_i_to_axes(xyz, axis="axis3", unit=[0, 0, 1], m=12.0107):
    """
    Align coordinates to axis
    """
    m = np.float(m)
    totm = m * xyz.shape[0]
    coms = com(xyz, m, totm)
    I = inertia(xyz, coms, m)
    e_values, e_vectors = np.linalg.eig(I)  # inertia)
    order = np.argsort(e_values)
    axis3, axis2, axis1 = e_vectors[:, order].transpose()
    axischoice = {"axis1": axis1, "axis2": axis2, "axis3": axis3}
    rotxyz = get_rotation_matrix(axischoice[axis], unit=unit)
    # Rotation first
    xyz = xyz - coms
    for i in range(xyz.shape[0]):
        xyz[i] = np.dot(xyz[i].T, rotxyz.T)
    return xyz


def direction(ca, chains, ids, cwd:Path):
    """
    Determine direction of subunits
    """
    xyz = ca.coord.copy()
    xyz = align_i_to_axes(xyz)
    for ind,atom in enumerate(ca.get_atoms()):
        atom.coord = np.array([round(xyz[ind][0], 3),round(xyz[ind][1], 3),round(xyz[ind][2], 3)])
    writer=PDBIO()
    writer.set_structure(ca)
    writer.save(cwd/f"aln_{ids}")

    coms = []
    for ch in chains:
        for s_chain in xyz.get_chains():
            if s_chain.id == ch:
                chcom = np.mean(ch.coords, 0)
                coms.append(chcom)
    edges = []
    for i in range(len(coms)):
        if i == len(coms) - 1:
            nextc = coms[0]
        else:
            nextc = coms[i + 1]
        c = coms[i]
        edges.append((nextc[0] - c[0]) * (nextc[1] + c[1]))
    mult = 1
    if xyz.coords[0][2] > np.mean(xyz.coords, 0)[2]:
        mult = -1
    else:
        mult = 1
    s = sum(edges) * mult
    if s < 0:
        return "ac"
    else:
        return "c"


def change_direction(chlist):
    """
    Reverse the chain order
    """
    chlist = list(sorted(chlist))
    rev = list(reversed(chlist))
    # To make sure we always keep the lowest letter fixed.
    fixed = chlist[0:1] + rev[0:-1]
    return fixed


def direction_check(uni_list, chlist,cwd:Path):
    """
    Check and determine which complexes to reorder
    """
    directions = {}
    for ids, uni in uni_list.items():
        directions[ids] = direction(uni, chlist[ids], ids, cwd)
    acs = list(directions.values()).count("ac")
    cs = list(directions.values()).count("c")
    reorder = []
    if acs > cs:
        for ids in directions:
            if directions[ids] == "c":
                reorder.append(ids)
    else:
        for ids in directions:
            if directions[ids] == "ac":
                reorder.append(ids)

    neworders = {}
    for ids in reorder:
        neworders[ids] = change_direction(chlist[ids])
    return neworders


def sepchains(ca, chlist):
    """
    Obtain coordinates of individual chains
    """
    sep = {}
    for ch in chlist:
        sep[ch] = ca[0][0][ch].copy()
    return sep


def reorder(id_ch_dict, cwd, tag=None):
    """
    Reorder the chains
    """
    uni_list = {}
    for ids in list(id_ch_dict.keys()):
        if tag:
            name = tag + ids
        else:
            name = ids
        frch = parse_chlist(open(cwd / name, "r"))
        frch = [i for i in frch if i in id_ch_dict[ids]]
        uni_list[ids] = getca_forchains(coord(open(cwd/name).readlines()),frch)
        id_ch_dict[ids] = frch

    neworders = direction_check(uni_list, id_ch_dict,cwd)
    for ids in uni_list:
        # We could simply loop over neworders improve this part.
        if ids in neworders.keys():
            newchains = []
            sep_chains = sepchains(uni_list[ids], id_ch_dict[ids])
            oldchs = id_ch_dict[ids]
            for ind, ch in enumerate(neworders[ids]):
                oldch = oldchs[ind]
                tmp = deepcopy(sep_chains[ch])
                changelock = np.where(tmp["ch"] == ch)
                tmp["ch"][changelock] = oldch
                newchains.append(tmp)
            merged = np.concatenate(newchains)
            os.rename(cwd/f"aln_{ids}",cwd/f"old_{name}")
            writeca(merged, cwd / name)
        else:
            os.rename(cwd/f"aln_{ids}",cwd/name)
    return list(neworders.keys())


def forceorder(id_ch_dict, cwd, alph=True, ordlist=[], forordict={}, tag=None):
    """
    Force an alphabetical or given order on the chains
    """
    if alph:
        for ids in ordlist:
            if tag:
                name = tag + ids
            else:
                name = ids
            newstamp = list(sorted(id_ch_dict[ids]))
            ca = getca_forchains(coord(open(cwd/name).readlines()),newstamp,order=False)
            sep_chains = sepchains(ca, newstamp)
            newchains = []
            for ch in newstamp:
                newchains.append(sep_chains[ch])
            merged = np.concatenate(newchains)
            os.rename(cwd/name,cwd/f"old_{name}")
            writeca(merged, cwd / name)
    else:
        for ids, ords in forordict.items():
            if tag:
                name = tag + ids
            else:
                name = ids
            newstamp = list(sorted(ords))
            ca = getca_forchains(coord(open(cwd/name).readlines()),id_ch_dict[ids],order=False)
            sep_chains = sepchains(ca, id_ch_dict[ids])
            newchains = []
            for ind, ch in enumerate(ords):
                tmp = deepcopy(sep_chains[ch])
                changelock = np.where(tmp["ch"] == ch)
                tmp["ch"][changelock] = newstamp[ind]
                newchains.append(tmp)
            merged = np.concatenate(newchains)
            os.rename(cwd/name,cwd/f"old_{name}")
            writeca(merged, cwd / name)
