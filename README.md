pdbParser aims to retrive CA coordinates of the structures deposited at the RCSB data bank. 

To make use of these scripts you need to install Python 2.7 and the packages listed below:

- numpy
- biopython

To use the standalone version simply run:
``` run.py --start PDBID --end PDBID ```

Use help to see other features. 

Optionally structures with no PDB headers could be cleaned too, however it is users responsibility to make sure that both structures don't contain any missing regions. Since the sequence is extracted from the supplied files, there is no way for me to know when both structures are missing the same residues. I don't trust the residue numbers; some structures contain a gap in the numbering altough they are intact simply because they were engineered. 

This will write two files start and end. These files will contain only the regions that are available in both structures. The core region is extracted via sequence alignment.

There are several issues when automating this process. The issues pdbParser can tackle at the moment are:
* Mismatches (Double check this feature, might be broken with the new update).
* Multiple biological assembiles are splitted and the first one with a complete set of chains is choosen. 
* Multimeric structures (Assuming the chain id's between two structures match heteromeric structures also works).

The support for below issues are currently being implemented:
* Insertions/deletions will result with partial structures. This is a bit tricky to handle since the gap penalty is set to extreme.
* PDB Headers could be problematic. Currently the program doesn't support mmCIF records. 
* Structures with missing residues will not be handled, will add a modeller script so people can fix the missing regions if they are reasonable for example not more then 3 consecutive residues are missing. 