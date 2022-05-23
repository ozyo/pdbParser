import logging
from enspdb.utils import ParserError
import requests
import urllib.request, urllib.parse, urllib.error

URLBASE = "http://www.uniprot.org/uniprot/"


def get_chainids(query):
    chainids = {
        z.split(";")[1].strip(): z.split(";")[4].split("=")[0].strip()
        for z in [
            i
            for i in urllib.request.urlopen(URLBASE + query + ".txt")
            .read()
            .decode("utf-8")
            .splitlines()
            if i.startswith("DR   PDB;")
        ]
        if z.split()[3] not in ["NMR;", "model;"]
    }
    return chainids


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


def get_pdbinfo(query, exclude, mer):

    idparam = {
        "query": "ID:{}".format(query),
        "format": "tab",
        "columns": "database(PDB),sequence",
    }

    result1 = requests.get(URLBASE, params=idparam).text
    if len(result1) > 0:
        pdbids, refseq = result1.split("\n")[1].split("\t")
        refseq = str(refseq)
        pdbids = ["{}".format(i) for i in pdbids.split(";") if len(i) > 1]
        pdbids = filter_ids(pdbids, exclude)
        chainids = get_chainids(query)
        returninfo = mer_info(pdbids, chainids, mer)
        return refseq, returninfo
    else:
        logging.critical("Cannot retrive the information for query number %s" % (query))
        return (None, None)


def get_query_info(querylist, exclude, mer):
    refseqs = {}
    pdbinfos = {}
    for query in querylist:
        seq, info = get_pdbinfo(query, exclude, mer)
        refseqs[query] = seq
        pdbinfos[query] = info
    return refseqs, pdbinfos
