import imp
import logging
import os
from pathlib import Path
from urllib import request

import numpy as np

from Bio.PDB.Selection import unfold_entities

from enspdb.alignment import aln_struct_to_core, getseq, parse_fasta_aln_multi
from enspdb import alignment as a
from enspdb.utils import SelectResidues, fold_entities, get_selected, write_pdb
from enspdb import utils as u
from importlib import reload
reload(u)
reload(a)
from enspdb.readpdb import check_pdb_title, getpdb, load_structure
from enspdb.clean_pdb import clean_altloc, getca_forchains
from enspdb import clean_pdb as c
reload(c)
from enspdb.uniprot_pdb_info import get_pdb_chainid_table, parse_refseq, parse_uniprot_query  


class PDBInfo:
    def __init__(self, query, mer, workdir , exclude=None):
        self.cwd = Path(workdir)
        self.cwd.mkdir(exist_ok=True)
        self.exclude = exclude
        self.query = [query] if isinstance(query, str) else query
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

    def get_pdbinfo(self):
        refseqs = {}
        returninfo = {}
        for query in self.query:
            refseq=parse_refseq(query)
            refseqs[query]=refseq
            query_result=parse_uniprot_query(query)
            pdbids=get_pdb_chainid_table(query_result)

            for pdb,chains in pdbids.items():
                count = 0
                if "/" in chains:
                    chains = chains.split("/")
                if len(chains) == self.mer:
                    returninfo[pdb] = [count + 1, [chains]]
                elif len(chains) > self.mer:
                    if len(chains) % self.mer == 0:
                        newchains = []
                        for chnr in range(0, len(chains), self.mer):
                            newchains.append(chains[chnr : chnr + self.mer])
                            count = count + 1
                        returninfo[pdb] = [count, newchains]
                    else:
                        logging.critical(
                            "Cannot process PDB id %s. It does not contain a complete set"
                            % pdb
                        )
                else:
                    logging.critical(
                        "Cannot process PDB id %s. It does not contain complete set"
                        % pdb
                    )
        if len(returninfo) == 0:
            return ({}, {})
        return (returninfo, refseqs)
    @staticmethod
    def downloadPDB(info_class, cwd: Path) -> None:
        pdb_dir = cwd / "rcsb"
        pdb_dir.mkdir(exist_ok=True)
        delete = []
        for pdb in info_class.result.keys():
            if not (pdb_dir / f"{pdb}.pdb").is_file():
                pdb_content = getpdb(Path(pdb), True, cwd=pdb_dir)
            else:
                pdb_content = getpdb(pdb_dir / f"{pdb}.pdb")
            if check_pdb_title(pdb_content.header["name"]) is True:
                delete.append(pdb)
                logging.critical(
                    "PDB ID %s cannot be processed, possibly a chimera, skipping this file"
                    % pdb
                )
                continue
        [info_class.result.pop(key) for key in delete]

    @staticmethod
    def process_pdbs(info_class, cwd: Path, overwrite_pdb: bool = False):
        pdb_dir = cwd / "rcsb"
        clean_pdb_dir = cwd / "clean_pdbs"
        clean_pdb_dir.mkdir(exist_ok=True)
        outseq = open(cwd / info_class.seqfilename, "w")
        outresmap = open(cwd / info_class.residmapfilename, "w")
        for query, seq in info_class.refseqs.items():
            outseq.write(f">refseq_{query}\n{seq}\n")
        for pdb in info_class.result.keys():
            pdb_content = load_structure(pdb_dir / f"{pdb}.pdb")
            for mol in range(0, info_class.result[pdb][0]):
                if (
                    overwrite_pdb
                    or not (clean_pdb_dir / f"{pdb}_{mol+1}.pdb").is_file()
                ):
                    ca = getca_forchains(pdb_content, info_class.result[pdb][1][mol])
                    ca = c.clean_altloc(ca, info_class.altloc)
                    write_pdb(ca, clean_pdb_dir / f"{pdb}_{mol+1}.pdb")
                else:
                    ca = load_structure(clean_pdb_dir / f"{pdb}_{mol+1}.pdb")
                for ch in info_class.result[pdb][1][mol]:
                    ca_ch = getca_forchains(ca, [ch])
                    seq, nr = a.getseq(ca_ch)
                    outseq.write(f">{pdb}_{mol+1}.pdb|{ch}|\n{seq}\n")
                    outresmap.write(
                        f">{pdb}_{mol+1}.pdb|{ch}|\n{'-'.join([str(i) for i in nr])}\n"
                    )
        outseq.close()
        outresmap.close()
    @staticmethod
    def msa_clustal(info_class, cwd, clustalopath, alnf=None, profile=None, cores=None):
        (
            info_class.coremer,
            info_class.coreresids,
            info_class.broken,
            refaln,
            structaln,
            info_class.core_blocks,
        ) = a.aln_struct_to_core(
            info_class.seqfilename,
            info_class.alnfasta,
            info_class.residmapfilename,
            cwd,
            info_class.result,
            clustalopath,
            profile=profile,
            alnfile=alnf,
            cores=cores,
        )
        return refaln, structaln
    @staticmethod
    def get_core(info_class,cwd):
        complete = info_class.coremer
        resids = info_class.coreresids
        broken = info_class.broken
        # totmer = info.mer
        filtresid = np.array([])
        fulllist = []
        for _, val in info_class.core_blocks.items():
            fulllist = fulllist + val
        filtresid = resids.iloc[:, fulllist]

        for pdb in complete:
            if pdb in broken:
                try:
                    os.rename(
                        cwd / "clean_pdbs" / pdb, cwd / "clean_pdbs" / f"broken_{pdb}"
                    )
                    continue
                except (OSError, IOError):
                    continue
            chains = complete[pdb]
            try:
                ca = load_structure(cwd / "clean_pdbs" / pdb)
            except (OSError, IOError):
                logging.warning("File does not exist: " + pdb + " skipped")
                continue
            # ca = getca_forchains(pdblines, [chains]), We are already writing an ensemble
            # No need to reread it
            # We have to reorder the chains since we loop over them.
            chains_resids = {ch.id:list(map(int,filtresid.loc[f"{pdb}|{ch.id}|"])) for ch in ca.get_chains()}
            newca=u.get_selected(ca, u.SelectResidues(chains_resids))
            newca = fold_entities(newca, "S")[0]
            write_pdb(newca, cwd / "clean_pdbs" / f"correct_{pdb}")
            #os.remove(f"{cwd}/clean_pdbs/{pdb}")

    def core_show(self, cwd, positions=[]):
        alndata, _ = parse_fasta_aln_multi(cwd / self.alnfasta)
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
