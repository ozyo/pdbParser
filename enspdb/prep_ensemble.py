import logging
import os
from pathlib import Path

import numpy as np
import urllib

from enspdb.clean_pdb import get_chain_order, getca_forchains
from enspdb.parser import parse_ca, pdb_title
from enspdb.readpdb import coord, getpdb
from enspdb.writepdb import writeca
from enspdb import alignment as a # import aln_struct_to_core, getseq, parse_fasta_aln_multi
from importlib import reload
reload(a)

class PDBInfo():
    def __init__(self, query, mer, exclude=None):
        self.exclude = exclude
        self.query = [query] if isinstance(query,str) else query
        self.mer = mer
        self.result, self.refseqs = self.get_pdbinfo()
        self.broken = None
        self.seqfilename = ".".join(self.query) + "_seq.txt"
        self.residmapfilename = ".".join(self.query) + "_residmap.txt"
        self.alnfasta = ".".join(self.query) + "_init.aln.txt"
        self.coremer = None
        self.coreresids = None

        self.altloc = "A"
        self.core_blocks = None
        
    def core_show(self, cwd, positions=[]):
        alndata, _ = a.parse_fasta_aln_multi(cwd/self.alnfasta)
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
        refseqs={}
        returninfo={}
        for query in self.query:
            URLbase = ('http://www.uniprot.org/uniprot/')

            idparam = {
                'query': 'ID:{}'.format(query),
                'format': 'tab',
                'columns': 'database(PDB),sequence'
            }
            idparam=urllib.parse.urlencode(idparam)
            result1 = urllib.request.Request(URLbase+'?'+idparam)
            try:
                result1=urllib.request.urlopen(result1).readlines()[1]
            except:
                logging.warning('Cannot retrive the information for query number %s' %(query))
                continue

            if len(result1) > 0:
                pdbids,refseq=result1.split(b'\t')
                refseqs[query]=refseq.decode('utf-8')
                pdbids=['{}'.format(i.decode('utf-8')) for i in pdbids.split(b';') if len(i)>1]

                if self.exclude is not None:
                    try:
                        for ex in self.exclude:
                            pdbids.remove(ex)
                    except ValueError:
                        pass

                chainids={}
                for i in urllib.request.urlopen(URLbase+query+'.txt').readlines():
                    if i.startswith (b'DR   PDB;') and i.split(b';')[3] not in ['NMR;','model;']:
                        chainids[i.split(b';')[1].strip().decode('utf-8')]=i.split(b';')[4].split(b'=')[0].strip().decode('utf-8')
            tmpids=pdbids.copy()
            for pdb in tmpids:
                count=0
                try:
                    chains=chainids[pdb]
                except KeyError:
                    logging.warning('PDB ID %s is either an NMR structure or a model. Skipping' %pdb)
                    continue
                if "/" in chains:
                    chains=chains.split('/')
                if len(chains) == self.mer:
                    returninfo[pdb]=[count+1,[chains]]
                elif len(chains) > self.mer:
                    if len(chains) % self.mer == 0:
                        newchains=[]
                        for chnr in range(0,len(chains),self.mer):
                            newchains.append(chains[chnr:chnr+self.mer])
                            count=count+1
                        returninfo[pdb]=[count,newchains]
                    else:
                        logging.critical('Cannot process PDB id %s. It does not contain a complete set' %pdb)
                else:
                    logging.critical('Cannot process PDB id %s. It does not contain complete set' %pdb)
        if len(returninfo) == 0:
            return(None,None) 
        return(returninfo,refseqs)

    def downloadPDB(self,cwd:Path):
        pdb_dir=cwd/"rcsb"
        pdb_dir.mkdir(exist_ok=True)
        delete = []
        for pdb in self.result.keys():
            if not (pdb_dir/f"{pdb}.pdb").is_file():
                pdb_content=getpdb(Path(pdb),True,pdb_dir)
                with open(pdb_dir/f"{pdb}.pdb","w") as outf:
                    outf.writelines(pdb_content)
            else:
                pdb_content=open(pdb_dir/f"{pdb}.pdb").readlines()
            if pdb_title(pdb_content) is True:
                delete.append(pdb)
                logging.critical('PDB ID %s cannot be processed, possibly a chimera, skipping this file' %pdb)
                continue
        [self.result.pop(key) for key in delete]

    def process_pdbs(self,cwd:Path,overwrite_pdb:bool=False):
        pdb_dir=cwd/"rcsb"
        clean_pdb_dir=cwd/"clean_pdbs"
        clean_pdb_dir.mkdir(exist_ok=True)
        outseq=open(cwd/self.seqfilename,'w')
        outresmap=open(cwd/self.residmapfilename,'w')
        for query,seq in self.refseqs.items():
            outseq.write(f'>refseq_{query}\n{seq}\n')
        for pdb in self.result.keys():
            pdb_content=open(pdb_dir/f"{pdb}.pdb").readlines()
            for mol in range(0,self.result[pdb][0]):
                if overwrite_pdb or not (clean_pdb_dir/f"{pdb}_{mol+1}.pdb").is_file():
                    ca = parse_ca(pdb_content, [self.result[pdb][1][mol]], self.altloc)
                    writeca(ca,clean_pdb_dir/f"{pdb}_{mol+1}.pdb")
                else:
                    ca = coord(open(clean_pdb_dir/f"{pdb}_{mol+1}.pdb").readlines())
                for ch in self.result[pdb][1][mol]:
                    ca_ch=getca_forchains(ca,[ch])
                    seq,mapx=a.getseq(ca_ch)
                    outseq.write('>'+pdb+'_'+str(mol+1)+'.pdb'+'|'+ch+'|'+'\n'+seq+'\n')
                    _code,_name,nr=zip(*mapx)
                    outresmap.write('>'+pdb+'_'+str(mol+1)+'.pdb'+'|'+ch+'|'+'\n'+'-'.join([str(i) for i in nr])+'\n')
        outseq.close()
        outresmap.close()        

    def msa_clustal(self,cwd,clustalopath,alnf=None,profile=None,cores=None):
        self.coremer, self.coreresids, self.broken, refaln,structaln,self.core_blocks = a.aln_struct_to_core(self.seqfilename,self.alnfasta,self.residmapfilename,cwd,self.result,clustalopath,
        profile=profile,alnfile=alnf,cores=cores)
        return refaln, structaln

    def get_core(self,cwd):
        complete = self.coremer
        resids = self.coreresids
        broken = self.broken
        # totmer = info.mer
        filtresid = np.array([])
        fulllist = []
        for _, val in self.core_blocks.items():
            fulllist = fulllist + val
        filtresid = resids.iloc[:, fulllist]

        for pdb in complete:
            if pdb in broken:
                try:
                    os.rename(cwd /"clean_pdbs"/pdb, cwd /"clean_pdbs"/f"broken_{pdb}")
                    continue
                except (OSError, IOError):
                    continue
            chains = complete[pdb]
            try:
                pdblines = open((cwd /"clean_pdbs"/ pdb), "r").readlines()
            except (OSError, IOError):
                logging.warning("File does not exist: " + pdb + " skipped")
                continue
            ca = getca_forchains(coord(pdblines), [chains])
            # We have to reorder the chains since we loop over them.
            chains = get_chain_order(ca)
            newca = None
            for ch in chains:
                tmp = filtresid.loc[f"{pdb}|{ch}|"].to_list()
                if newca is None:
                    newca = ca[(ca["ch"] == ch) & np.in1d(ca["resnr"], tmp)]
                else:
                    newca = np.concatenate([newca, ca[(ca["ch"] == ch) & np.in1d(ca["resnr"], tmp)]])
            writeca(newca, f"{cwd}/clean_pdbs/correct_{pdb}")
            os.remove(f"{cwd}/clean_pdbs/{pdb}")