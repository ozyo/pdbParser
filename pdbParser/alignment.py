# See COPYING for license
from pathlib import Path
from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio.SeqUtils import seq1 as letter
from Bio.Align.Applications import ClustalOmegaCommandline
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


def getaligned(ca1, ca2):
    sca1, mapca1 = getseq(ca1)
    sca2, mapca2 = getseq(ca2)
    # Grabs the first alignment
    aligned = align(sca1, sca2)[0]
    logging.info(str(aligned))
    aca1 = list(aligned[0])
    aca2 = list(aligned[1])
    shift1 = findgap(aca1)
    shift2 = findgap(aca2)
    # There is a case 1CIL-3EWV, where the sequences doesn't match but instead of printing an error
    # code fails with index error. That's because shifts are larger than the mapca.
    # It happens only when 1CIL is the starting structure.
    # To ensure a smoot integration with the server catching all types of index errors here.
    try:
        # map the gaps to the opposite structure, this is important!!!
        nter1 = mapca1[shift2[0] : shift2[1]][0]
        cter1 = mapca1[shift2[0] : shift2[1]][-1]
        nter2 = mapca2[shift1[0] : shift1[1]][0]
        cter2 = mapca2[shift1[0] : shift1[1]][-1]
    except IndexError:
        logging.error("Sequences contain large shifts. Please check the sequence alignment.")
        logging.error(
            "This usually happens when two unrelated proteins are given as start and target and the target structure is shorter than the start."
        )
        return [], [], False
    # since we mapped it to the opposite in the steps above, we can return back to normal
    core1 = ca1[(ca1["resnr"] >= nter1[2]) & (ca1["resnr"] <= cter1[2])]
    core2 = ca2[(ca2["resnr"] >= nter2[2]) & (ca2["resnr"] <= cter2[2])]
    if len(core1) == len(core2):
        return core1, core2, True
    else:
        sch = np.unique(core1["ch"])
        ech = np.unique(core2["ch"])
        logging.error(
            f"There are non-terminal missing residues missing in one of the chains {sch} and {ech} withiin start and target pdb files."
        )
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
    return whole1, whole2, all(correct)


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


def return_gaps(row:pd.DataFrame):
    return row[row=="-"].dropna(axis=1).columns.to_list()


def align_resid_to_seq(resmap,seqaln):
    resmap = open(resmap, "r").read().split(">")
    resids = seqaln.copy()
    broken=[]
    for entry in resmap:
        if len(entry) > 2:
            ids, nr = entry.split("\n")[0:2]
            nr = nr.split("-")
            gaps = return_gaps(seqaln.loc[[ids]])
            keep = (
                seqaln.loc[
                    [ids],
                ]
                .drop(gaps, axis=1)
                .columns
            )
            if len(nr) != len(keep):
                logging.warning("%s has incosistent sequence and residue numbering. Marked as broken." % ids)
                resids = resids.drop(ids, axis=0)
                broken.append(ids.split("|")[0])
                continue
            else:
                for ind, i in enumerate(keep):
                    resids[i][ids] = nr[ind]
    return broken,resids


def find_nter_gap(seqaln):
    pos = None
    for col in seqaln.columns:
        content = list(seqaln[col].unique())
        if "-" in content:
            pos = col
        else:
            return pos + 1

def find_cter_gap(seqaln):
    pos = None
    for col in list(seqaln.columns)[::-1]:
        content = list(seqaln[col].unique())
        if "-" in content:
            pos = col
        else:
            return pos - 1


def clean_empty_blocks(seqaln,blocks):
    block_keys=list(blocks.keys()).copy()
    for block in block_keys:
        filt = seqaln.loc[:, blocks[block]]
        missing = filt.apply(lambda x:1 if set("".join(x)) == {"-"} else 0,axis=1)
        if sum(missing) == len(filt.index):
            logging.info("Block %s is not resolved in any of the structures." % block)
            blocks.pop(block)


def find_resid_onetoone(seqaln, resmap, blocks):
    """
    Replaces the sequence with the residue numbers and keeps the gaps etc. Also filters the broken structures.
    This is an easier function to use than the above functions. The previous ones are kept since they are working but
    this is much better and easier to follow. So in the future port everything to pandas.
    """
    broken, resids = align_resid_to_seq(resmap,seqaln)
    clean_empty_blocks(seqaln,blocks)
    block_keys = list(sorted(blocks.keys()))
    for block in block_keys:
        start = blocks[block][0]
        end = blocks[block][-1]
        nter = find_nter_gap(seqaln.loc[:,start:end])
        cter = find_cter_gap(seqaln.loc[:,start:end])
        to_remove = []
        if nter is not None and  blocks[block][0] < nter:
            start = nter
            for i in blocks[block]:
                if i < nter:
                    to_remove.append(i)
                else:
                    break
        if cter is not None and blocks[block][-1] > cter:
            end = cter
            for i in  blocks[block][::-1]:
                if i > cter:
                    to_remove.append(i)
                else:
                    break
        blocks[block]=[i for i in blocks[block] if i not in to_remove]
        filt = resids.loc[:, range(start,end+1)]
        tmpbroken = []
        for ids in filt.index:
            if ids in broken:
                continue
            else:
                gaps = return_gaps(filt.loc[[ids]])
                if len(gaps) == len(filt.columns) or len(gaps) >= len(filt.columns) / 2.0:
                    tmpbroken.append(ids)
                else:
                    if len(gaps) > 0:
                        logging.warning("%s has missing residues in block %s, marked as broken." % (ids, block))
                        broken.append(ids.split("|")[0])
        if len(tmpbroken) > 0:
            for tmp in tmpbroken:
                logging.warning("%s has missing residues in block %s, marked as broken" % (tmp, block))
                broken.append(tmp.split("|")[0])
    return (resids, list(set(broken)))


def aln_struct_to_core(infile:str, outfile:str, resmap, cwd:Path, merinfo, clustalopath, cores=None, profile=None,alnfile=None):
    if alnfile is None:
        if profile is None:
            clustalomega_cline = ClustalOmegaCommandline( cmd=clustalopath,
                    infile=(cwd/infile).as_posix(), outfile=(cwd/outfile).as_posix(), verbose=False, auto=True, force=True
                )
            clustalomega_cline()
            alndata,_ = parse_fasta_aln_multi((cwd/outfile).as_posix())
        else:
            clustalomega_cline = ClustalOmegaCommandline( cmd=clustalopath,
                infile=(cwd / infile).as_posix(),
                profile1=(cwd / profile).as_posix(),
                outfile=(cwd / outfile).as_posix(),
                verbose=False,
                auto=True,
                force=True,
                )
            clustalomega_cline()
            alndata, _ = parse_fasta_aln_multi(cwd / outfile)
    else:
        alndata,_=parse_fasta_aln_multi(cwd/alnfile)
    refaln = alndata.filter(regex="refseq_", axis=0)
    structaln = alndata.filter(regex=".pdb", axis=0)
    if cores is None:
        core = find_core(refaln)
    else:
        core = cores
    resid, broken = find_resid_onetoone(structaln, cwd / resmap, core)
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
