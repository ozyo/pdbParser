import logging
from enspdb.utils import ParserError
import urllib.request
from Bio import SeqIO

URLBASE = "http://www.uniprot.org/uniprot/"


def get_pdb_chainid_table(query_result):
    pdb_chainid={}
    for line in query_result:
        line=line.decode("utf-8")
        if line.startswith("DR   PDB;"):
            pdb_chainid[line.split(";")[1].strip()]=line.split(";")[4].split("=")[0].strip()
    return pdb_chainid


def filter_ids(pdbids, exclude):
    for ex in exclude:
        try:
            pdbids.remove(ex)
            logging.info(f"Removing PDB ID {ex}")
        except ValueError:
            logging.warning(f"PDB ID {ex} not in the list.")
    return pdbids


def split_mer(mer, chlist):
    """
    Split chain ids into mers.
    """
    if mer == len(chlist):
        return [chlist]
    elif len(chlist) > mer or mer == 1:
        if len(chlist) % mer == 0:
            moldivision = []
            for i in range(0, len(chlist), mer):
                moldivision.append(chlist[i : i + mer])
            return moldivision
        else:
            try:
                raise ValueError(
                    "Chain labels are not equal or multiple of total number of chains given."
                )
            except ValueError as err:
                logging.exception("FAIL", exc_info=err)
                raise ParserError
    else:
        try:
            raise ValueError(f"List of chain ids: {chlist} cannot be split to {mer}.")
        except ValueError as err:
            logging.exception("FAIL", exc_info=err)
            raise ParserError


def mer_info(pdbids, chainids, mer):
    returninfo = {}
    for pdb in pdbids:
        chains = chainids[pdb].split("/")
        mers = split_mer(mer, chains)
        returninfo[pdb] = [len(mers), mers]
    return returninfo


def parse_uniprot_query(query):

    URL = f"https://rest.uniprot.org/uniprotkb/{query}.txt"
    result = urllib.request.urlopen(URL).readlines()
    return result

def parse_refseq(query):
    URL = f"https://rest.uniprot.org/uniprotkb/{query}.fasta"
    result = SeqIO.read(urllib.request.urlretrieve(URL)[0],'fasta')
    return str(result.seq)

def get_query_info(querylist, exclude, mer):
    refseqs = {}
    pdbinfos = {}
    for query in querylist:
        seq, info = get_pdbinfo(query, exclude, mer)
        refseqs[query] = seq
        pdbinfos[query] = info
    return refseqs, pdbinfos
