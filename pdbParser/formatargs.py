# See COPYING for license
import logging
from pdbParser.utils import ParserError


def format_ch(ch_arg):
    """
    Fortmat chain arguments
    """
    chlist = []
    tmplist = ch_arg.split(",")
    for sub in tmplist:
        if len(sub) == 3 and sub[1] == "-":
            sublist = sub.split("-")
            firstch = sublist[0]
            lastch = sublist[1]
            finsublist = [chr(i) for i in range(ord(firstch), ord(lastch) + 1)]
            for ch in finsublist:
                chlist.append(ch)
        elif len(sub) == 1 and isinstance(sub, str):
            chlist.append(sub)
        else:
            try:
                raise ValueError("Please supply the chain ids in these formats: A-E or A-C,E or A,B")
            except ValueError as err:
                logging.error("FAIL", exc_info=err)
                raise ParserError

    return chlist


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
                raise ValueError("Chain labels are not equal or multiple of total number of chains given.")
            except ValueError as err:
                logging.exception("FAIL", exc_info=err)
                raise ParserError
    else:
        try:
            raise ValueError(f"List of chain ids: {chlist} cannot be split to {mer}.")
        except ValueError as err:
            logging.exception("FAIL", exc_info=err)
            raise ParserError
