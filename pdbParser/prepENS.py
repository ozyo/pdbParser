import requests
import logging
import urllib.request, urllib.parse, urllib.error
import sys,os
import numpy as np
from importlib import reload
from pdbParser import parser as pP
reload(pP)
from pdbParser import writepdb as wp
from pdbParser import clean_pdb as cp
from pdbParser import alignment as a
reload(a)
import pandas as pd

class PDBInfo():
    def __init__(self,query,mer,exclude=None):
        self.exclude=exclude
        self.query=query
        self.mer=mer
        self.result,self.refseq=self.get_pdbinfo()
        self.broken=None
        self.seqfilename=query+'_seq.txt'
        self.residmapfilename=query+'_residmap.txt'
        self.alnfasta=query+'_init.aln.txt'
        self.coremer=None
        self.coreresids=None

    def get_pdbinfo(self):
        URLbase = ('http://www.uniprot.org/uniprot/')

        idparam = {
            'query': 'ID:{}'.format(self.query),
            'format': 'tab',
            'columns': 'database(PDB),sequence'
        }

        result1 = requests.get(URLbase,params=idparam).text
        if len(result1) > 0:
            pdbids,refseq=result1.split('\n')[1].split('\t')
            refseq=str(refseq)
            pdbids=['{}'.format(i) for i in pdbids.split(';') if len(i)>1]
            if self.exclude is not None:
                try:
                    for ex in self.exclude:
                        pdbids.remove(ex)
                        logging.info('Removing PDB ID %s' %ex)
                except ValueError:
                    logging.warning('Did not find the PDB ID %s' %ex)
            chainids={z.split(';')[1].strip():z.split(';')[4].split('=')[0].strip() for z in [i for i in urllib.request.urlopen(URLbase+self.query+'.txt').read().decode('utf-8').splitlines() if i.startswith('DR   PDB;')] if z.split()[3] not in ['NMR;','model;']}
        else:
            logging.critical('Cannot retrive the information for query number %s' %(self.query))
            return(None,None)
        returninfo={}
        for pdb in pdbids:
            try:
                count=0
                try:
                    chains=chainids[pdb].split('/')
                except AttributeError:
                    continue
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
                        logging.critical('Cannot process PDB id %s. It does not contain complete set' %pdb)
                        chainids.pop(pdb)
                        #pdbids.remove(pdb)
                else:
                    logging.critical('Cannot process PDB id %s. It does not contain complete set' %pdb)
                    chainids.pop(pdb)
                    #pdbids.remove(pdb)
            except KeyError:
                logging.warning('PDB ID %s is either an NMR structure or a model. Skipping' %pdb)
                #pdbids.remove(pdb)
        return(returninfo,refseq)

    def core_show(self,cwd,positions=[]):
        alndata,alncomment=a.parse_fasta_aln_multi(cwd+'/'+self.alnfasta)
        self.alndata=alndata
        if len(positions)==2:
            alndata=self.alndata.iloc[:,positions[0]:positions[-1]]
        elif len(positions) == 0:
            pass
        else:
            logging.error('Please provide a start and an end number in a list')
            return None
        #def here_block(s,se):
        #    return ['background-color: yellow' if v in se else '' for v in s.index]
        def here_gap(val):
            color='red' if val == '-' else ''
            return 'background-color: %s' % color    
        styled=alndata.style.applymap(here_gap)
        return(styled)

class multiPDBInfo():
    def __init__(self,querylist,mer,tag,exclude=None):
        self.exclude=exclude
        self.querylist=querylist
        self.tag=tag
        self.mer=mer
        self.result,self.refseq=self.get_pdbinfo_multi()
        self.broken=None
        self.seqfilename=self.tag+'_seq.txt'
        self.residmapfilename=self.tag+'_residmap.txt'
        self.alnfasta=self.tag+'_init.aln.txt'
        self.coremer=None
        self.coreresids=None

    def get_pdbinfo_multi(self):
        URLbase = ('http://www.uniprot.org/uniprot/')
        fullpdblist={}
        fullrefseqs={}
        for query in self.querylist:
            idparam = {
                'query': 'ID:{}'.format(query),
                'format': 'tab',
                'columns': 'database(PDB),sequence'
            }

            result1 = requests.get(URLbase,params=idparam).text

            if len(result1) > 0:
                pdbids,refseq=result1.split('\n')[1].split('\t')
                refseq=str(refseq)
                fullrefseqs[query]=refseq
                pdbids=['{}'.format(i) for i in pdbids.split(';') if len(i)>1]
                if self.exclude is not None:
                    try:
                        for ex in self.exclude:
                            pdbids.remove(ex)
                            logging.info('Removing PDB ID %s' %ex)
                    except ValueError:
                        logging.warning('Did not find the PDB ID %s' %ex)

                chainids={z.split(';')[1].strip():z.split(';')[4].split('=')[0].strip() for z in [i for i in urllib.request.urlopen(URLbase+query+'.txt').read().decode('utf-8').splitlines() if i.startswith('DR   PDB;')] if z.split()[3] not in ['NMR;','model;']}
                print(chainids)
            else:
                logging.critical('Cannot retrive the information for query number %s' %(query))
                return(None,None)
            for pdb in pdbids:
                try:
                    fullpdblist[pdb]
                except KeyError:
                    try:
                        fullpdblist[pdb]=chainids[pdb]
                    except KeyError:
                        print(f"{pdb} does not have chain id skipped")
                else:
                    fullpdblist[pdb]=fullpdblist[pdb]+'/'+chainids[pdb]

        returninfo={}
        for pdb in list(fullpdblist.keys()):
            try:
                count=0
                try:
                    chains=fullpdblist[pdb].split('/')
                except AttributeError:
                    continue
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
                        logging.critical('Cannot process PDB id %s. It does not contain complete set' %pdb)
                        fullpdblist.pop(pdb)
                else:
                    logging.critical('Cannot process PDB id %s. It does not contain complete set' %pdb)
                    fullpdblist.pop(pdb)
            except KeyError:
                logging.warning('PDB ID %s is either an NMR structure or a model. Skipping' %pdb)
                fullpdblist.pop(pdb)
        return(returninfo,fullrefseqs)

    def core_show(self,seqs=['ref','struct'],positions=[],blocks=[]):
        if seqs == ['ref','struct']:
            filt=pd.concat([self.refseqaln,self.structaln])
        elif seqs == ['struct']:
            filt=self.structaln
        elif seqs == ['ref']:
            filt=self.refseqaln
        else:
            logging.error('Please choose from ["ref","struct"] or one of them ["ref"], ["struct"]')
            return None
        if len(blocks)==0:
            blocks=self.core_blocks.keys()
        else:
            blocks=blocks
        if len(positions)==2:
            filt=filt.iloc[:,positions[0]:positions[-1]]
        elif len(positions) == 0:
            filt=filt
        else:
            logging.error('Please provide a start and an end number in a list')
            return None
        se=[]
        for block in blocks:
            se=se+self.core_blocks[block]
        def here_block(s,se):
            return ['background-color: yellow' if v in se else '' for v in s.index]
        def here_gap(val):
            color='red' if val == '-' else ''
            return 'background-color: %s' % color    
        styled=filt.style.apply(here_block,se=se,axis=1).applymap(here_gap)
        return(styled)

    def write_seq_resid(self):
        aanr={}
        
        return(styled)


def downloadPDB(info,cwd,multiseq=False,returnordch=True):
    query=info.query if not multiseq else info.querylist
    pdblist=info.result
    mer=info.mer
    refseq=info.refseq
    altloc="A"
    outseq=open(cwd+'/'+info.seqfilename,'w')
    outresmap=open(cwd+'/'+info.residmapfilename,'w')
    orderedchinfo={pdb:pdblist[pdb][0] for pdb in pdblist}
    orderedchinfo={key:[val,[0]*val] for key,val in orderedchinfo.items()}
    if multiseq:
        for query in refseq.keys():
            outseq.write('>refseq_'+query+'\n'+refseq[query]+'\n')
    else:
        outseq.write('>refseq'+'\n'+refseq+'\n')
    for pdb in list(pdblist.keys()):
        print(pdb)
        try:
            pdblines=open(cwd+'/'+pdb+'.pdb').readlines()
        except:
            try:
                urllib.request.urlretrieve('http://files.rcsb.org/download/%s.pdb' %pdb, cwd+'/'+pdb+'.pdb')
            except:
                print("%s cannot be downloaded as pdb" %pdb)
                continue
            else:
                pdblines=open(cwd+'/'+pdb+'.pdb').readlines()
        #time.sleep(-1)
        for mol in range(pdblist[pdb][0]):
            #pdblines=open(cwd+'/'+pdb+'.pdb').readlines()
            if pP.pdbTitle(pdblines) is True:
                try:
                    pdblist.pop(pdb)
                    logging.critical('PDB ID %s is a chimera, skipping this file' %pdb)
                    continue
                except KeyError:
                    logging.info('This structure was already removed. I am an example of bad programming. Nothing to worry about')
                    continue
            coord,ordch=pP.parse_ca(pdblines,pdb,mer,altloc,chains=pdblist[pdb][1][mol],caonly=True,returnch=True)
            orderedchinfo[pdb][1][mol]=ordch
            wp.writeca(coord,cwd+'/'+pdb+'_'+str(mol+1)+'.pdb')
            for ch in orderedchinfo[pdb][1][mol]:
                ca=cp.getca(coord,altloc,ch)
                seq,maps=a.getseq(ca)
                outseq.write('>'+pdb+'_'+str(mol+1)+'.pdb'+'|'+ch+'|'+'\n'+seq+'\n')
                code,name,nr=list(zip(*maps))
                outresmap.write('>'+pdb+'_'+str(mol+1)+'.pdb'+'|'+ch+'|'+'\n'+'-'.join([str(i) for i in nr])+'\n')
        #os.remove(cwd+'/'+pdb+'.pdb')
    outseq.close()
    outresmap.close()
    return(orderedchinfo)



def msa(info,cwd,clustalopath,alnf=None,multiseq=False,updates=False,cores=None):
    query=info.query if not multiseq else info.tag
    if multiseq and not alnf:
        logging.error('If using multiple UniProt IDs, profile alignment file is necesseary.')
        return(None)
    seqfile=info.seqfilename
    resmap=info.residmapfilename
    merinfo=info.result
    totmer=info.mer
    outfile=info.alnfasta
    if not multiseq:
        complete,resids,broken=a.msa_clustal(seqfile,resmap,outfile,clustalopath,cwd,merinfo,query,totmer,alnf)
    else:
        complete,resids,broken,refaln,structaln,blocks=a.aln_struct_to_core(alnf,outfile,seqfile,resmap,cwd,merinfo,query,totmer,clustalopath,updates=updates,cores=cores)
        info.refseqaln=refaln
        info.structaln=structaln
        if not updates:
            info.core_blocks=blocks
    info.broken=broken
    info.coremer=complete
    info.coreresids=resids
    
def getcore(info,cwd,multiseq=False):
    altloc='A'
    complete=info.coremer
    resids=info.coreresids
    broken=info.broken
    totmer=info.mer
    if multiseq:
        fulllist=[]
        for key,val in info.core_blocks.items():
            fulllist=fulllist+val
        filtresid=resids.iloc[:,fulllist]

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
        ca=pP.parse_ca(pdblines,pdb,totmer,altloc,[chains],caonly=True)
        #We have to reorder the chains since we loop over them.
        _, idx = np.unique(ca['ch'], return_index=True)
        chains=ca['ch'][np.sort(idx)]
        newca=None
        for ch in chains:
            if not multiseq:
                nter,cter=[int(i) for i in resids[pdb+'|'+ch+'|']]
                if newca is None:
                    newca=ca[(ca['ch']==ch)&(ca['resnr']>=nter) & (ca['resnr']<=cter)]
                else:
                    newca=np.concatenate([newca,ca[(ca['ch']==ch)&(ca['resnr']>=nter) & (ca['resnr']<=cter)]])
            else:
                tmp=filtresid.loc[pdb+'|'+ch+'|'].to_list()
                if newca is None:
                    newca=ca[(ca['ch']==ch)&np.in1d(ca['resnr'], tmp)]
                else:
                    newca=np.concatenate([newca,ca[(ca['ch']==ch)&np.in1d(ca['resnr'], tmp)]])
        wp.writeca(newca,cwd+'/'+'correct_'+pdb)
        os.remove(cwd+'/'+pdb)
