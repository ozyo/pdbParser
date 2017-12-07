
#See COPYING for license 

#We need to read quite a bit of information from the REMARKS lines. So it is better to keep them in a seperate module.
import numpy as np
import logging
import numpy.lib.recfunctions
import os

def remove_field_name(a, name):
    names = list(a.dtype.names)
    if name in names:
        names.remove(name)
    b = a[names]
    return b

class Segsep(object):
    def __init__(self,coord):
        self.hoh=coord[(coord['resname'] == 'HOH') | (coord['resname'] == 'TIP3') | (coord['resname'] == 'SOL')]
        self.rest=coord[(coord['resname'] != 'HOH') & (coord['resname'] != 'TIP3') | (coord['resname'] != 'SOL') ]
        self.chains=set(self.rest['ch'].tolist())
        self.segids=set(self.rest['segid'].tolist())
        self.coorsegs=np.empty_like(coord)
        remove=['element','charge']
        for name in remove:
            self.coorsegs=remove_field_name(self.coorsegs,name)
        #self.coorsegs=numpy.lib.recfunctions.append_fields(self.coorsegs, 'segid', dtypes='S5', usemask=False, asrecarray=True)
        if len(self.chains) == 0:
            self.chains=['A']
            print 'Hey no chain found, I set everything to A'

    def sep_segs(self,coord):
        if len(self.segids) == 0:
            remove=['element','charge']
            for chain in self.chains:
                seg_ch=self.rest[self.rest['ch']==chain]
                seg_less=seg_ch
                for name in remove:
                    seg_less=remove_field_name(seg_less,name)
                seg_id=numpy.lib.recfunctions.append_fields(seg_less, 'segid', ['SEG'+chain]*len(seg_less), dtypes='S5', usemask=False, asrecarray=True)
                self.coorsegs=np.concatenate(self.coorsegs,seg_id)
            seg_less_wat=self.hoh
            if len(seg_less_wat) > 0:
                for name in remove:
                    seg_less_wat=remove_field_name(seg_less_wat,name)
                seg_id_wat=numpy.lib.recfunctions.append_fields(seg_less_wat, 'segid', ['SEGW']*len(seg_less_wat), dtypes='S5', usemask=False, asrecarray=True)
                self.coorsegs=np.concontenate(self.coorsegs,seg_id_wat)
            return self.coorsegs
        else:
            return coord
