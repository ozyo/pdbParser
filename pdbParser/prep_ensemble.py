import logging
import os
from pathlib import Path

import numpy as np
import urllib

from pdbParser.clean_pdb import getca_forchains
from pdbParser.parser import parse_ca, pdb_title
from pdbParser.readpdb import getpdb
from pdbParser.writepdb import writeca
from pdbParser.alignment import getseq, msa_clustal, parse_fasta_aln_multi

class PDBInfo():
    def __init__(self, query, mer, exclude=None):
        self.exclude = exclude
        self.query = query
        self.mer = mer
        self.result, self.refseq = self.get_pdbinfo()
        self.broken = None
        self.seqfilename = query + "_seq.txt"
        self.residmapfilename = query + "_residmap.txt"
        self.alnfasta = query + "_init.aln.txt"
        self.coremer = None
        self.coreresids = None
        self.altloc = "A"
        
    def core_show(self, cwd, positions=[]):
        alndata, _ = parse_fasta_aln_multi(cwd + "/" + self.alnfasta)
        self.alndata = alndata
        if len(positions) == 2:
            alndata = self.alndata.iloc[:, positions[0] : positions[-1]]
        elif len(positions) == 0:
            pass
        else:
            logging.error("Please provide a start and an end number in a list")
            return None
        # def here_block(s,se):
        #    return ['background-color: yellow' if v in se else '' for v in s.index]
        def here_gap(val):
            color = "red" if val == "-" else ""
            return "background-color: %s" % color

        styled = alndata.style.applymap(here_gap)
        return styled

    def get_pdbinfo(self):
        URLbase = ('http://www.uniprot.org/uniprot/')

        idparam = {
            'query': 'ID:{}'.format(self.query),
            'format': 'tab',
            'columns': 'database(PDB),sequence'
        }
        idparam=urllib.parse.urlencode(idparam)
        result1 = urllib.request.Request(URLbase+'?'+idparam)
        try:
            result1=urllib.request.urlopen(result1).readlines()[1]
        except:
            logging.error('Cannot retrive the information for query number %s' %(self.query))
            return(None,None)

        if len(result1) > 0:
            pdbids,refseq=result1.split(b'\t')
            refseq=str(refseq)
            pdbids=['{}'.format(i.decode('utf-8')) for i in pdbids.split(b';') if len(i)>1]

            if self.exclude is not None:
                try:
                    for ex in self.exclude:
                        pdbids.remove(ex)
                        logging.info('Removing PDB ID %s' %ex)
                except KeyError:
                    logging.warning('Did not find the PDB ID %s' %ex)

            chainids={}
            for i in urllib.request.urlopen(URLbase+self.query+'.txt').readlines():
                if i.startswith (b'DR   PDB;') and i.split(b';')[3] not in ['NMR;','model;']:
                    chainids[i.split(b';')[1].strip().decode('utf-8')]=i.split(b';')[4].split(b'=')[0].strip().decode('utf-8')
        returninfo={}
        tmpids=[*pdbids]
        for pdb in tmpids:
            count=0
            try:
                chains=chainids[pdb]
            except KeyError:
                logging.warning('PDB ID %s is either an NMR structure or a model. Skipping' %pdb)
                pdbids.remove(pdb)
                continue
            try:
                chains=chains.split('/')
            except AttributeError:
                continue
            else:
                nchain=len(chains)
                if nchain == self.mer:
                    returninfo[pdb]=[count+1,[chains]]
                elif nchain > self.mer:
                    if nchain % self.mer == 0:
                        newchains=[]
                        for chnr in range(0,nchain,self.mer):
                            newchains.append(chains[chnr:chnr+self.mer])
                            count=count+1
                        returninfo[pdb]=[count,newchains]
                    else:
                        logging.critical('Cannot process PDB id %s. It does not contain a complete set' %pdb)
                        chainids.pop(pdb)
                        pdbids.remove(pdb)
                else:
                    logging.critical('Cannot process PDB id %s. It does not contain complete set' %pdb)
                    chainids.pop(pdb)
                    pdbids.remove(pdb)
        if len(returninfo) == 0:
            return(None,None) 
        return(returninfo,refseq)

    def downloadPDB(self, pdb_dir:Path):
        pdb_dir.mkdir(exist_ok=True)
        delete = []
        for pdb in self.result.keys():
            pdb_content=getpdb(Path(pdb),True,pdb_dir)
            if pdb_title(pdb_content) is True:
                try:
                    delete.append(pdb)
                    logging.critical('PDB ID %s contains cannot be processed, possibly a chimera, skipping this file' %pdb)
                    continue
                except KeyError:
                    logging.info('This structure was already removed. I am an example of bad programming. Nothing to worry about')
                    continue
            else:
                for mol in range(0,self.result[pdb][0]):
                    ca = parse_ca(pdb_content, [self.result[pdb][1][mol]], self.altloc)
                    writeca(ca,pdb_dir/f"{pdb}_{mol+1}.pdb")
        [self.result.pop(key) for key in delete]

    def write_chain_sequence(self,cwd,pdb_dir):
        outseq=open(cwd+'/'+self.seqfilename,'w')
        outresmap=open(cwd+'/'+self.residmapfilename,'w')
        outseq.write('>refseq'+'\n'+self.refseq+'\n')
        for pdb in self.result.keys():
            pdb_content = getpdb(pdb,False,cwd=pdb_dir)
            for mol in range(0,self.result[pdb][0]):
                for ch in self.result[pdb][1][mol]:
                    ca=getca_forchains(pdb_content,self.altloc,ch)
                    seq,mapx=getseq(ca)
                    outseq.write('>'+pdb+'_'+str(mol+1)+'.pdb'+'|'+ch+'|'+'\n'+seq+'\n')
                    code,name,nr=zip(*mapx)
                    outresmap.write('>'+pdb+'_'+str(mol+1)+'.pdb'+'|'+ch+'|'+'\n'+'-'.join([str(i) for i in nr])+'\n')
        outseq.close()
        outresmap.close()

    def msa(self,cwd,clustalopath,alnf=None):
        outfile=cwd+'/'+self.alnfasta
        self.coremer,self.coreresids,self.broken=msa_clustal(self.seqfilename,self.residmapfilename,outfile,clustalopath,cwd,self.result,self.query,alnf)
    
    def getcore(self,cwd):
        complete=self.coremer
        resids=self.coreresids
        broken=self.broken
        for pdb in complete:
            if pdb in broken:
                try:
                    os.rename(cwd+'/'+pdb,cwd+'/'+"broken_"+pdb)
                    continue
                except (OSError,IOError):
                    continue
            chains=complete[pdb]
            try:
                pdblines=open(cwd+'/'+pdb,'r').readlines()
            except (OSError,IOError):
                logging.warning('File does not exist: '+pdb+' skipped')
                continue
            ca=getca_forchains(pdblines,[chains],order=False)
            newca=None
            for ch in chains:
                nter,cter=[int(i) for i in resids[pdb+'|'+ch+'|']]
                if newca is None:
                    newca=ca[(ca['ch']==ch)&(ca['resnr']>=nter) & (ca['resnr']<=cter)]
                else:
                    newca=np.concatenate([newca,ca[(ca['ch']==ch)&(ca['resnr']>=nter) & (ca['resnr']<=cter)]])
            writeca(newca,cwd+'/'+'correct_'+pdb)
            os.remove(cwd+'/'+pdb)
