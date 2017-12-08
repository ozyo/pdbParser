pdbParser aims to retrive CA coordinates of the structures deposited at the RCSB data bank to use with eBDIMs

To make use of these scripts you need to install Python 2.7 and the packages listed below:

- numpy
- biopython

To use the standalone version simply run:
``` run.py --start PDBID --end PDBID ```

This will write two files start and end. These files will contain only the regions that are available in both structures. The core region is extracted via sequence alignment.

The issues pdbParser can't tackle at the moment are:
* Multiple biological assemblies (Not even started)
* PDB Headers could be problematic. Currently the program doesn't support mmCIF records. (Disabled at the moment)
* Structures with missing residues are not repaired. (I might write a modeller input file)

Below cases should give correct structures but caution is needed

* Support for multimeric assemblies 
* Insertions/deletions may result with partial structures, testing extreme cases but should be resolved. This is a bit tricky to handle since the gap penalty is set to extreme.
* Mismatches shouldn't create a problem. Each structure's Calpha is extracted.


pdbIO branch is for reading and writing CHARMM segments mainly. It will also allow manipulation of the Bfact field, could be useful for eBDIMS or writing NAMD restraints.
There is one option for writing pdb file in the way CRYSOL likes.