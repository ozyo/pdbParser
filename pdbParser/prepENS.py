import logging
import urllib.request, urllib.parse, urllib.error
import os
import numpy as np
from importlib import reload
from pdbParser import parser as pP

reload(pP)
from pdbParser import writepdb as wp
from pdbParser import clean_pdb as cp
from pdbParser import alignment as a

reload(a)
import pandas as pd


class PDBInfo:
    def __init__(self, query, mer, exclude=None):
        self.exclude = exclude
        self.query = query
        self.mer = mer
        self.result, self.refseq = None, None
        self.broken = None
        self.seqfilename = query + "_seq.txt"
        self.residmapfilename = query + "_residmap.txt"
        self.alnfasta = query + "_init.aln.txt"
        self.coremer = None
        self.coreresids = None

    def core_show(self, cwd, positions=[]):
        alndata, _ = a.parse_fasta_aln_multi(cwd + "/" + self.alnfasta)
        self.alndata = alndata
        if len(positions) == 2:
            alndata = self.alndata.iloc[:, positions[0] : positions[-1]]
        elif len(positions) == 0:
            pass
        else:
            logging.error("Please provide a start and an end number in a list")
            return None
        # def here_block(s,se):
        #    return ['background-color: yellow' if v in se else '' for v in s.index]
        def here_gap(val):
            color = "red" if val == "-" else ""
            return "background-color: %s" % color

        styled = alndata.style.applymap(here_gap)
        return styled


class multiPDBInfo:
    def __init__(self, querylist, mer, tag, exclude=None):
        self.exclude = exclude
        self.querylist = querylist
        self.tag = tag
        self.mer = mer
        self.result, self.refseq = self.get_pdbinfo_multi()
        self.broken = None
        self.seqfilename = self.tag + "_seq.txt"
        self.residmapfilename = self.tag + "_residmap.txt"
        self.alnfasta = self.tag + "_init.aln.txt"
        self.coremer = None
        self.coreresids = None
        self.refseqaln = None
        self.structaln = None
        self.core_blocks = dict()

    def core_show(self, seqs=["ref", "struct"], positions=[], blocks=[]):
        if seqs == ["ref", "struct"]:
            filt = pd.concat([self.refseqaln, self.structaln])
        elif seqs == ["struct"]:
            filt = self.structaln
        elif seqs == ["ref"]:
            filt = self.refseqaln
        else:
            logging.error('Please choose from ["ref","struct"] or one of them ["ref"], ["struct"]')
            return None
        if len(blocks) == 0:
            blocks = self.core_blocks.keys()
        else:
            blocks = blocks
        if len(positions) == 2:
            filt = filt.iloc[:, positions[0] : positions[-1]]
        elif len(positions) == 0:
            filt = filt
        else:
            logging.error("Please provide a start and an end number in a list")
            return None
        se = []
        for block in blocks:
            se = se + self.core_blocks[block]

        def here_block(s, se):
            return ["background-color: yellow" if v in se else "" for v in s.index]

        def here_gap(val):
            color = "red" if val == "-" else ""
            return "background-color: %s" % color

        styled = filt.style.apply(here_block, se=se, axis=1).applymap(here_gap)
        return styled


def write_pdb_seq(pdblist, refseq, seq_filename, residmap_filename, cwd, download=False):
    altloc = "A"
    outseq = open(cwd + "/" + seq_filename, "w")
    outresmap = open(cwd + "/" + residmap_filename, "w")
    orderedchinfo = {}
    for query in refseq.keys():
        outseq.write(">refseq_" + query + "\n" + refseq[query] + "\n")

    for query, pdbs in pdblist.items():
        for pdb in pdbs:
            if download:
                urllib.request.urlretrieve("http://files.rcsb.org/download/%s.pdb" % pdb, cwd + "/" + pdb + ".pdb")
            pdblines = open(cwd + "/" + pdb + ".pdb").readlines()
            for mol in range(pdbs[pdb][0]):
                if pP.pdb_title(pdblines) is True:
                    logging.critical("PDB ID %s is a chimera, skipping this file" % pdb)
                    pdbs.pop(pdb)
                    pdblist[query] = pdbs
                    break
                coord = pP.parse_ca(pdblines, pdblist[pdb][1][mol], altloc)
                ordch = pP.parse_chlist(coord)
                orderedchinfo[pdb] = ordch
                wp.writeca(coord, cwd + "/" + pdb + "_" + str(mol + 1) + ".pdb")
                for ch in ordch:
                    ca = cp.getca_forchains(coord, altloc, ch)
                    seq, maps = a.getseq(ca)
                    outseq.write(f"> {pdb}_{mol + 1}.pdb|{ch}|\n")
                    outseq.write(f"{seq}\n")
                    # code, name, nr = list(zip(*maps))
                    _, _, nr = list(zip(*maps))
                    outresmap.write(f">{pdb}_{mol + 1}.pdb|{ch}|\n")
                    outresmap.write(f"{'-'.join([str(i) for i in nr])}\n")
        # os.remove(cwd+'/'+pdb+'.pdb')
    outseq.close()
    outresmap.close()
    return orderedchinfo


def msa(info, cwd, clustalopath, alnf=None, multiseq=False, updates=False, cores=None):
    query = info.query if not multiseq else info.tag
    if multiseq and not alnf:
        logging.error("If using multiple UniProt IDs, profile alignment file is necesseary.")
        return None
    seqfile = info.seqfilename
    resmap = info.residmapfilename
    merinfo = info.result
    totmer = info.mer
    outfile = info.alnfasta
    if not multiseq:
        complete, resids, broken = a.msa_clustal(
            seqfile, resmap, outfile, clustalopath, cwd, merinfo, query, totmer, alnf
        )
    else:
        complete, resids, broken, refaln, structaln, blocks = a.aln_struct_to_core(
            alnf, outfile, seqfile, resmap, cwd, merinfo, query, totmer, clustalopath, updates=updates, cores=cores
        )
        info.refseqaln = refaln
        info.structaln = structaln
        if not updates:
            info.core_blocks = blocks
    info.broken = broken
    info.coremer = complete
    info.coreresids = resids


def getcore(info, cwd, multiseq=False):
    altloc = "A"
    complete = info.coremer
    resids = info.coreresids
    broken = info.broken
    # totmer = info.mer
    filtresid = np.array([])
    if multiseq:
        fulllist = []
        for _, val in info.core_blocks.items():
            fulllist = fulllist + val
        filtresid = resids.iloc[:, fulllist]

    for pdb in complete:
        if pdb in broken:
            try:
                os.rename(cwd + "/" + pdb, cwd + "/" + "broken_" + pdb)
                continue
            except (OSError, IOError):
                continue
        chains = complete[pdb]
        try:
            pdblines = open(cwd + "/" + pdb, "r").readlines()
        except (OSError, IOError):
            logging.warning("File does not exist: " + pdb + " skipped")
            continue
        ca = pP.parse_ca(pdblines, [chains], altloc)
        # We have to reorder the chains since we loop over them.
        chains = cp.get_chain_order(ca)
        newca = None
        for ch in chains:
            if not multiseq:
                nter, cter = [int(i) for i in resids[pdb + "|" + ch + "|"]]
                if newca is None:
                    newca = ca[(ca["ch"] == ch) & (ca["resnr"] >= nter) & (ca["resnr"] <= cter)]
                else:
                    newca = np.concatenate(
                        [newca, ca[(ca["ch"] == ch) & (ca["resnr"] >= nter) & (ca["resnr"] <= cter)]]
                    )
            else:
                tmp = filtresid.loc[f"{pdb}|{ch}|"].to_list()
                if newca is None:
                    newca = ca[(ca["ch"] == ch) & np.in1d(ca["resnr"], tmp)]
                else:
                    newca = np.concatenate([newca, ca[(ca["ch"] == ch) & np.in1d(ca["resnr"], tmp)]])
        wp.writeca(newca, f"{cwd}/correct_{pdb}")
        os.remove(f"{cwd}/{pdb}")
