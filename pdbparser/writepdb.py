#See COPYING for license 

import numpy as np

def crysol_sanity_check(coord):
    for location, atom in enumerate(coord['atname']):
        if len(atom) == 1:
            coord['atname'][location]=atom+'  '
        elif len(atom) == 2:
            coord['atname'][location]=atom+' '
        # This is a bit dangerous, for crysol I should be removing the drudes for analysis.
        # So remove it for now. Otherwise things are named the same altough they are not!
        #elif len(atom) == 4:
        #    coord['atname'][location]=atom[1:3]
    return coord

def write(coord,file,coortype,bysegid):
    print "Writing coordinates"
    if bysegid is False:
        with open(file,'w') as output:
            if coortype == "charmm":
                np.savetxt(output,coord,fmt='ATOM  %5s %-4s%1s%-4s%1s%5s   %8s%8s%8s%6s%6s      %-5s',delimiter='')
            elif coortype == "crysol":
                #This is needed to keep the atom name the way crysol likes. '  C ' or '  N1' not 'N1  ' or '   C'
                print "Writing crysol format"
                coord=crysol_sanity_check(coord)
                np.savetxt(output,coord,fmt='ATOM  %5s %4s%1s%-4s%1s%5s   %8s%8s%8s%6s%6s     %-5s',delimiter='')
            else:
                coord=crysol_sanity_check(coord)
                np.savetxt(output,coord,fmt='ATOM  %5s %4s%1s%-4s%1s%5s   %8s%8s%8s%6s%6s          %-2s%2s',delimiter='')
            np.savetxt(output,np.array('END').reshape(1,),fmt='%3s',delimiter='s')
    else:
        # Add the END stament to the end of the file
        segids=np.unique(coord['segid'])
        for seg in segids:
            to_write=coord[coord['segid']==seg]
            with open(seg.lower()+'.pdb','w') as output:
                np.savetxt(output,to_write,fmt='ATOM  %5s %-4s%1s%-4s%1s%5s   %8s%8s%8s%6s%6s      %-5s',delimiter='')
                np.savetxt(output,np.array('END').reshape(1,),fmt='%3s',delimiter='s')
