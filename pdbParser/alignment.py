# See COPYING for license
# Things to do
# obtain files from previous processes and create a sequence alignment between CA atoms
from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio.SeqUtils import seq1 as letter
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
import numpy as np
import pandas as pd
from itertools import groupby
from operator import itemgetter
import logging

# tags: a aligned, s sequence, map residue name-nr map

matrix = substitution_matrices.load("BLOSUM62")


def getseq(ca):
    """
    Extract sequence from
    """
    seq = letter("".join(ca["resname"].tolist()))
    resmap = list(zip(list(seq), ca["resname"].tolist(), ca["resnr"].tolist()))
    return seq, resmap


def align(seq1, seq2):
    # basestring is a left over from python3
    if isinstance(seq1, str):
        alignments = pairwise2.align.globalds(
            seq1,
            seq2,
            matrix,
            -11,
            -11,
            penalize_end_gaps=True,
        )
        return alignments
    elif isinstance(seq1, list):
        alignments = pairwise2.align.globalds(seq1, seq2, matrix, -11, -11, gap_char=["-"], penalize_end_gaps=True)
        return alignments


def findgap(aca):
    # This function finds the first patch without gap. For an MSA we need something that removes the structure if it is broken.
    if aca[0] == "-" and aca[-1] != "-":
        start = next(ind for ind, gap in enumerate(aca) if gap != "-")
        end = len(aca) + 1
    elif aca[-1] == "-" and aca[0] != "-":
        end = aca.index("-") - len(aca)  # -1
        start = 0
    elif aca[0] == "-" and aca[-1] == "-":
        start = next(ind for ind, gap in enumerate(aca) if gap != "-")
        try:
            end = aca[start:-1].index("-") + start - 1
        except ValueError:  # If only the last is - and the ones before are a sequence
            end = aca[start:].index("-") + start - 1
    else:
        start = 0
        end = len(aca) + 1
    return start, end


def findinsertions(aca):
    # Use this function with a reference sequence alignment result to find insertions
    i = 0
    start = i
    end = start + 1
    for i in range(0, len(aca)):
        if aca[i] != "-":
            start = i
            break
    subaca = aca[i:]
    lastelem = len(subaca) - 1
    while lastelem >= 0:
        if subaca[lastelem] != "-":
            end = lastelem + start
            break
        lastelem = lastelem - 1
    insertions = []
    for i in range(start, end + 1):
        if aca[i] == "-":
            insertions.append(i)
    return insertions


def findtergap(aca, insertions):
    i = 0
    start = i
    end = start + 1
    for i in range(0, len(aca)):
        if aca[i] != "-":
            start = i
            break
    subaca = aca[i:]
    lastelem = len(subaca) - 1
    while lastelem >= 0:
        if subaca[lastelem] != "-":
            end = lastelem + start
            break
        lastelem = lastelem - 1
    broken = False
    insertion = False
    for i in range(start, end + 1):
        if i in insertions and aca[i] == "-":
            continue
        elif i in insertions and aca[i] != "-":
            insertion = True
            break
        else:
            if aca[i] == "-":
                broken = True
                break
    if broken is True:
        return None, None
    if insertion is True:
        return None, None
    else:
        return start, end


def getaligned(ca1, ca2):
    sca1, mapca1 = getseq(ca1)
    sca2, mapca2 = getseq(ca2)
    # Grabs the first alignment
    aligned = align(sca1, sca2)[0]
    print(aligned)
    aca1 = list(aligned[0])
    aca2 = list(aligned[1])
    shift1 = findgap(aca1)
    shift2 = findgap(aca2)
    # map the gaps to the opposite structure, this is important!!!
    nter1 = mapca1[shift2[0] : shift2[1]][0]
    cter1 = mapca1[shift2[0] : shift2[1]][-1]
    nter2 = mapca2[shift1[0] : shift1[1]][0]
    cter2 = mapca2[shift1[0] : shift1[1]][-1]
    # since we mapped it to the opposite in the steps above, we can return back to normal
    core1 = ca1[(ca1["resnr"] >= nter1[2]) & (ca1["resnr"] <= cter1[2])]
    core2 = ca2[(ca2["resnr"] >= nter2[2]) & (ca2["resnr"] <= cter2[2])]
    if len(core1) == len(core2):
        logging.info("Run successful can proceed eBDIMS calculation")
        logging.critical("SUCCESS")
        return core1, core2, True
    else:
        logging.critical("There are non-terminal missing residues.")
        logging.critical("Please fix the structures manually to continue.")
        logging.critical("FAIL")
        return core1, core2, False


def multialigned(ca1, ca2, mer):
    whole1 = None
    whole2 = None
    correct = []
    if len(np.unique(ca1["ch"])) == len(np.unique(ca2["ch"])) and len(np.unique(ca1["ch"])) == mer:
        for a, b in zip(np.unique(ca1["ch"]), np.unique(ca2["ch"])):
            cores = getaligned(ca1[ca1["ch"] == a], ca2[ca2["ch"] == b])
            if whole1 is None:
                whole1 = cores[0]
                whole2 = cores[1]
                correct.append(cores[2])
            elif whole1 is not None:
                whole1 = np.append(whole1, cores[0], axis=0)
                whole2 = np.append(whole2, cores[1], axis=0)
                correct.append(cores[2])
    if False in correct:
        logging.critical("There are non-terminal missing residues.")
        logging.critical("Please fix the structures manually to continue.")
        logging.critical("FAIL")
        return whole1, whole2, False
    else:
        return whole1, whole2, True


def findresid(shifts, nter, cter, resmap):
    resmap = open(resmap, "r").read().split(">")
    resid = {}
    for pdb in resmap:
        if len(pdb) > 2:
            ids, nr = pdb.split("\n")[0:2]
            if ids in shifts.keys():
                if len(shifts[ids]) == 0:
                    resid[ids] = [None, None]
                    continue
                noffset = nter - shifts[ids][0]
                coffset = cter - nter + noffset
                nrs = nr.split("-")[noffset : coffset + 1]
                try:
                    resid[ids] = [int(nrs[0]), int(nrs[-1])]
                except IndexError:
                    logging.critical(
                        "PDB ID %s contains too little sequence or the alignment is problematic. Remove this structure and try again."
                        % id
                    )
                    exit()
    return resid


def msa_clustal(infile, resmap, outfile, clustalopath, cwd, merinfo, query, totmer, alnf=None):
    resmap = cwd + "/" + resmap
    if alnf is None:
        clustalomega_cline = ClustalOmegaCommandline(
            infile=cwd + "/" + infile, outfile=cwd + "/" + outfile, verbose=False, auto=True, force=True
        )
        clustalomega_cline()
        msa = AlignIO.read(cwd + "/" + outfile, "fasta")
    else:
        msa = AlignIO.read(cwd + "/" + alnf, "fasta")
    broken = []
    nter = 0
    cter = msa.get_alignment_length()
    shifts = {}
    insertions = []
    for record in msa:
        if record.id == "refseq":
            insertions = findinsertions(str(record.seq))
            break
    for record in msa:
        if record.id == "refseq":
            continue
        aln = str(record.seq)
        start, end = findtergap(aln, insertions)
        if start is None and end is None:
            broken.append(record.id.split("|")[0])
            shifts[record.id] = []
            logging.critical(
                "%s contains insertions or missing residues, skipping this chain and the assembly it belongs to."
                % record.id
            )
            continue
        else:
            shifts[record.id] = [start, end]
            if start is not None and start > nter:
                nter = start
            if end is not None and end < cter:
                cter = end
    core = msa[:, nter : cter + 1]
    resid = findresid(shifts, nter, cter, resmap)
    completemers = {}
    fullids = list(set([key.split("|")[0] for key in resid.keys()]))
    for pdb in fullids:
        pdbid, mer = pdb.split(".")[0].split("_")
        amer = []
        # merinfo[pdbid] #What this line was supposed to be ?
        for ch in merinfo[pdbid][1][int(mer) - 1]:
            if pdb + "|" + ch + "|" in broken:
                continue
            else:
                amer.append(ch)
        completemers[pdb] = amer
    AlignIO.write(core, cwd + "/" + query + "_core.fasta", format="fasta")
    return (completemers, resid, list(set(broken)))


def parse_fasta_aln_multi(alnf):
    """
    Parses an alignment file generated with sequences from multiple organisms and the structures.
    """
    alndata = {}
    alncomment = {}
    ids = ""
    with open(alnf, "r") as aln:
        for line in aln:
            if line.startswith(">"):
                ids = line.split()[0][1:]
                comment = line.split()[1:]
                alndata[ids] = ""
                alncomment[ids] = comment
            else:
                alndata[ids] = alndata[ids] + line.rstrip()
    for key, values in alndata.items():
        alndata[key] = list(values)
    alndata = pd.DataFrame(alndata).T
    return (alndata, alncomment)


def find_core(refaln):
    # https://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
    # https://stackoverflow.com/questions/56704700/select-columns-if-any-of-their-rows-contain-a-certain-string/56704768
    missing = refaln.columns[refaln.select_dtypes(object).applymap(lambda x: "-" in str(x)).any()].tolist()
    keep = refaln.drop(missing, axis=1)
    blocks = []
    for _, g in groupby(enumerate(keep.columns), lambda ix: ix[0] - ix[1]):
        blocks.append(list(map(itemgetter(1), g)))
    filtblocks = {}
    ind = 1
    for b in blocks:
        if len(b) >= 5:
            filtblocks["B" + str(ind)] = b
            ind += 1
    return filtblocks


def find_resid_onetoone(structaln, resmap, blocks):
    """
    Replaces the sequence with the residue numbers and keeps the gaps etc. Also filters the broken structures.
    This is an easier function to use than the above functions. The previous ones are kept since they are working but
    this is much better and easier to follow. So in the future port everything to pandas.
    """
    resids = structaln.copy()
    resmap = open(resmap, "r").read().split(">")
    broken = []
    for entry in resmap:
        if len(entry) > 2:
            ids, nr = entry.split("\n")[0:2]
            nr = nr.split("-")
            gaps = structaln.loc[[ids]][structaln == "-"].dropna(axis=1).columns.to_list()
            keep = (
                structaln.loc[
                    [ids],
                ]
                .drop(gaps, axis=1)
                .columns
            )
            if len(nr) != len(keep):
                logging.error("%s has incosistent sequence and residue numbering. Marked as broken." % ids)
                resids = resids.drop(ids, axis=0)
                broken.append(ids.split("|")[0])
                continue
            else:
                for ind, i in enumerate(keep):
                    resids[i][ids] = nr[ind]

    for block in list(blocks.keys()):
        start = blocks[block][0]
        end = blocks[block][-1] + 1
        filt = resids.iloc[:, start:end]
        missing = 0
        tmpbroken = []
        for ids in filt.index:
            if ids in broken:
                continue
            else:
                gaps = filt.loc[[ids]][filt == "-"].dropna(axis=1).columns.to_list()
                if len(gaps) == len(filt.columns) or len(gaps) >= len(filt.columns) / 2.0:
                    missing += 1
                    tmpbroken.append(ids)
                else:
                    if len(gaps) > 0:
                        logging.warning("%s has missing residues in block %s, marked as broken." % (ids, block))
                        broken.append(ids.split("|")[0])
        if missing == len(filt.index):
            logging.warning("Block %s is not resolved in any of the structures." % block)
            blocks.pop(block)
        else:
            if len(tmpbroken) > 0:
                for tmp in tmpbroken:
                    logging.warning("%s has missing residues in block %s, marked as broken" % (tmp, block))
                    broken.append(tmp.split("|")[0])
    return (resids, list(set(broken)))


def aln_struct_to_core(alnf, outf, seqf, resmap, cwd, merinfo, query, totmer, clustalopath, updates=False, cores=None):
    clustalomega_cline = ClustalOmegaCommandline(
        infile=cwd + "/" + seqf,
        profile1=cwd + "/" + alnf,
        outfile=cwd + "/" + outf,
        verbose=False,
        auto=True,
        force=True,
    )
    clustalomega_cline()
    alndata, _ = parse_fasta_aln_multi(cwd + "/" + outf)
    refaln = alndata.filter(regex="refseq_", axis=0)
    structaln = alndata.filter(regex=".pdb", axis=0)
    if not updates:
        core = find_core(refaln)
    else:
        core = cores
    resid, broken = find_resid_onetoone(structaln, cwd + "/" + resmap, core)
    completemers = {}
    fullids = list(set([key.split("|")[0] for key in resid.index]))
    for pdb in fullids:
        pdbid, mer = pdb.split(".")[0].split("_")
        amer = []
        # merinfo[pdbid] #What this line was supposed to be ?
        for ch in merinfo[pdbid][1][int(mer) - 1]:
            if pdb + "|" + ch + "|" in broken:
                continue
            else:
                amer.append(ch)
        completemers[pdb] = amer
    return (completemers, resid, broken, refaln, structaln, core)
