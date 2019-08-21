
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
        self.whole=coord
        self.hoh=coord[(coord['resname'] == 'HOH') | (coord['resname'] == 'TIP3') | (coord['resname'] == 'SOL')]
        self.rest=coord[(coord['resname'] != 'HOH') & (coord['resname'] != 'TIP3') & (coord['resname'] != 'SOL') ]
        self.chains=set(self.rest['ch'].tolist())
        try:
            len(self.rest['segid'])
            self.segids=set(self.rest['segid'].tolist())
        except ValueError:
            'No segments found will create a new array with segment field and transfer the coordinates'
        try:
            len(self.rest['segid'])
            self.segids=set(self.rest['segid'].tolist())
        except ValueError:
            'No segments found will create a new array with segment field and transfer the coordinates'
        
        self.coorsegs=np.empty_like(coord)
        remove=['element','charge']
        for name in remove:
            self.coorsegs=remove_field_name(self.coorsegs,name)
        self.coorsegs=numpy.lib.recfunctions.append_fields(self.coorsegs, 'segid',[]*len(self.coorsegs), dtypes='S5', usemask=False, asrecarray=True)
        if len(self.chains) == 0:
            self.chains=['A']
            print 'Hey no chain found, I set everything to A, your segid will SEGA'
        
    def sep_segs(self,coord):
        try:
            self.segids
        except AttributeError:
            remove=['element','charge']
            for chain in self.chains:
                seg_ch=self.rest[self.rest['ch']==chain]
                seg_less=seg_ch
                for name in remove:
                    seg_less=remove_field_name(seg_less,name)
                seg_id=numpy.lib.recfunctions.append_fields(seg_less, 'segid', ['SEG'+chain]*len(seg_less), dtypes='S5', usemask=False, asrecarray=True)
                self.coorsegs=np.concatenate((self.coorsegs,seg_id))
            seg_less_wat=self.hoh
            if len(seg_less_wat) > 0:
                for name in remove:
                    seg_less_wat=remove_field_name(seg_less_wat,name)
                seg_id_wat=numpy.lib.recfunctions.append_fields(seg_less_wat, 'segid', ['SEGW']*len(seg_less_wat), dtypes='S5', usemask=False, asrecarray=True)
                self.coorsegs=np.concatenate((self.coorsegs,seg_id_wat))
            self.coorsegs=self.coorsegs[self.coorsegs['segid']!='N/A']
            return self.coorsegs
        else:
            return coord

    def add_segid(self,coord):
        try:
            self.segids
        except AttributeError:
            remove=['element','charge']
            for chain in self.chains:
                seg_ch=self.whole[self.whole['ch']==chain]
                seg_less=seg_ch
                for name in remove:
                    seg_less=remove_field_name(seg_less,name)
                seg_id=numpy.lib.recfunctions.append_fields(seg_less, 'segid', ['SEG'+chain]*len(seg_less), dtypes='S5', usemask=False, asrecarray=True)
                self.coorsegs=np.concatenate((self.coorsegs,seg_id))
            self.coorsegs=self.coorsegs[self.coorsegs['segid']!='N/A']
            return self.coorsegs
