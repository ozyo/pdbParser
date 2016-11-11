import numpy as np
import logging

def missinginfo(r465,compnd,ca):
    missinginfo={}
    if len(r465) != 0:
        logging.info('Identified missing residues.')
        logging.info('Checking the integrity of the structure')
        for ch in compnd:
            rmis=r465[r465['ch']==ch]['rid'].tolist()
            max=ca[ca['ch']== ch ]['resnr'].tolist()[-1]
            min=ca[ca['ch']== ch ]['resnr'].tolist()[0]
            for i in rmis:
                if i < max and i > min:
                    logging.critical('If there are multiple biological assemblies, I will try to pick the one that is not broken. Else, I will quit')
                    missinginfo[ch]=0
                    continue
            missinginfo[ch]=1
    else:
        logging.info('No missing residues detected.')
        for ch in chain:
            missinginfo[ch]=1
    return missinginfo
