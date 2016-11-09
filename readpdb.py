#We need to read quite a bit of information from the REMARKS lines. So it is better to keep them in a seperate module.
import numpy as np

class pdbread(object):
    def readCompnd(self,pdb):
        compnd=None
        read=False
        for line in pdb:
            if 'COMPND' in line:
                if 'MOL_ID: 1' in line:
                    read=True
                elif 'MOL_ID: 2' in line:
                    read=False
                if read is True and 'CHAIN:' in line:
                    compnd=[i.strip().strip(';') for i in line.split(':')[1].split(',')]
        return compnd

    def readMissing(self,pdb):
        compnd=self.readCompnd(pdb)
        remark=[]
        read=False
        for line in pdb:
            if 'REMARK 465' in line:
                remark.append(line)
        r465=np.genfromtxt(remark[7:-1],names=['REMARK','465','rname','ch','rid'],dtype=['S6',int,'S3','S1',int])
        filt=r465[np.in1d(r465['ch'],compnd)]
        return filt

    def readAtom(self,pdb):
        atoms=[]
        read=False
        for line in pdb:
            if line.startswith('ATOM'):
                atoms.append(line)
        return atoms

    def coord (self,pdb):
        coords=[]
        atoms=self.readAtom(pdb)
        for atom in atoms:
            atnr=int(atom[6:11].strip())
            atname=atom[12:16].strip()
            altloc=atom[16]
            resname=atom[17:20].strip()
            ch=atom[21]
            resnr=atom[22:26]
            icode=atom[26]
            x=float(atom[30:38].strip())
            y=float(atom[38:46].strip())
            z=float(atom[46:54].strip())
            occu=float(atom[54:60].strip())
            tfact=float(atom[60:66].strip())
            element=atom[76:78]
            charge=atom[78:80]
            coords.append([atnr,atname,altloc,resname,ch,resnr,icode,x,y,z,occu,tfact,element,charge])
        coords=np.array(coords)
        return coords

#pdb=open('3RIF.pdb').readlines()
#read=pdbread()
#filt=read.readRemark(pdb)
#print filt
