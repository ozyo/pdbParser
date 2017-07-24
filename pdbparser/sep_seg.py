
#See COPYING for license 

#We need to read quite a bit of information from the REMARKS lines. So it is better to keep them in a seperate module.
import numpy as np
import logging
from writepdb import writecharmm, writecharmm_noicode
import numpy.lib.recfunctions
import os

def remove_field_name(a, name):
    names = list(a.dtype.names)
    if name in names:
        names.remove(name)
    b = a[names]
    return b

def replace(a,old,new):
    for x,y in zip(old,new):
        ind=a['resname'] == x
        a['resname'][ind] = y
    return a

class Segsep(object):
    def __init__(self,coord):
        self.hoh=coord[(coord['resname'] == 'HOH') | (coord['resname'] == 'TIP3') | (coord['resname'] == 'SOL')]
        self.rest=coord[(coord['resname'] != 'HOH') & (coord['resname'] != 'TIP3') | (coord['resname'] != 'SOL') ]
        self.chains=set(self.rest['ch'].tolist())
        print self.chains
        self.segids=set(coord['segid'].tolist())
    def sep_segs(self,coord,cwd,rname):
        #remove=['element','charge']
        remove=['segid']
        for chain in self.chains:
            seg_ch=self.rest[self.rest['ch']==chain]
            seg_less=seg_ch
            for name in remove:
                seg_less=remove_field_name(seg_less,name)
            seg_id=numpy.lib.recfunctions.append_fields(seg_less, 'segid', ['SEG'+chain]*len(seg_less), dtypes='S4', usemask=False, asrecarray=True)
            old=['A','C','G','T']
            new=['ADE','CYT','GUA','THY']
            seg_new=replace(seg_id,old,new)
            #seg_final=remove_field_name(seg_new,"icode")
            writecharmm_noicode(seg_new,cwd+'/seg'+chain.lower()+'.pdb')
        seg_less_wat=self.hoh
        if len(seg_less_wat) > 0 and rname is True:
            for name in remove:
                seg_less_wat=remove_field_name(seg_less_wat,name)
            seg_id_wat=numpy.lib.recfunctions.append_fields(seg_less_wat, 'segid', ['SEGW']*len(seg_less_wat), dtypes='S4', usemask=False, asrecarray=True)
            seg_new_wat=replace(seg_id_wat,['HOH'],['TIP3'])
            #seg_final_wat=remove_field_name(seg_new_wat,"icode")
            writecharmm_noicode(seg_new_wat,cwd+'/segw'+'.pdb')
        else:
            writecharmm_noicode(seg_less_wat,cwd+'/segw'+'.pdb')
    def sep_bysegids(self,coord,cwd):
        for segid in self.segids:
            print segid
            seg=coord[coord['segid']==segid]
            print seg[0]
            writecharmm_noicode(seg,cwd+'/'+segid.lower()+'.pdb')

