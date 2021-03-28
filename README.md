# pdbParser for eBDIMS

pdbParser cleans the PDB files and extracts the CA coordinates to be used with eBDIMs. Code can be run on local files, or with PDB codes to download files from the RCSB server. 

## REQUIREMENTS

- numpy
- biopython
- requests
- argparse

## USAGE

To download structures from the RCSB database that has a single chain `A`.

``` pdbParser-cmd --start PDBID --end PDBID --download```

To process multimeric complexes one has to pass the relevant chain ids and how many chains make the complex:

``` pdbParser-cmd --start PDBID --end PDBID --schains A-D --echains F-I --multimeric 2```

Above examples informs the code that chains A to D belongs to the starting complex, and A-B and C-D are making up this complex. It is possible that there are multiple biological assemblies in a PDB file. Within these assembles while A-B contains missing residues C-D could be intact. The code will automatically scan through these assemblies until it finds an intact one that can be used.
Length of chain ids given must be divisible by `multimeric` argument.

### Ensemble Mode

**WARNING** This mode has not been upgraded and tested with python 3.8 yet.

Ensemble mode requires Uniprot IDs to collect information about all the relevant PDB files. This processes is not completely automated and might still require user intervention, thus it is recommended to use this mode from an interpreter and assess the complexes.

``` pdbParser-cmd --prepENS UniProt ID --multimeric Nr.ofChains ```

## OUTPUTS

For the non-ensemble mode two files, start.pdb and target.pdb will be written. These files will contain only the regions that are available in both structures. Terminal mismatches and terminal missing regions will be trimmed. The core region is extracted via sequence alignment that relies on `biopython`.The sequence alignment is printed to the terminal for the user to consult if there any issues. The alignment penalizes the gaps heavily, so any single point mutations is likely to be aligned correctly.

Non-terminal gaps are not fixed and they will result with an error, no file will be written. The user can consult example Modeller scripts on tutoril or any other program to fix the starting structures before running this code.

## ROADMAP

* Ensmble preparation to be upgraded to python 3.8 (supporing heteromeric structures, multiple Uniprot IDs)

* Extract partial structures when insertions/deletions are present. There might be cases where a user wants to proceed with broken structures. For example, the missing region might belong to an intracellular domain of a membrane protein that can be disregarded for quick analysis.

* Read HETATM lines to support non-natural amino acids that contain a canonical backbone. (Check of calcium atoms and clean them in an extra step)

* Add unit tests with multimeric assemblies

* Add a special mode where MOLID type exclusions can be ignored. To allow processing of complexes crystallized with antibodies keywords like `antibody,ScFv` is used to ignore them. However there might be cases where this leads to an empty complex if the MOLID titles are not handled correctly.
