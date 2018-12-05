pdbParser retrives CA coordinates of the structures deposited at the RCSB data bank to be used with eBDIMs

To make use of these scripts you need to install Python 2.7 and the packages listed below:

- numpy
- biopython
- requests

To use the standalone version simply run (either with local files or with PDB ID):

``` run.py --start PDBID --end PDBID ```

This will write two files start and end. These files will contain only the regions that are available in both structures. The core region is extracted via sequence alignment. There are options to chose the chains, or the alternative location. If your structure is multimeric and any of the files contains multiple assemblies you should provide the number of chains.

There is also a support for ensemble preperation. You can refer to the tutorial in the github page for the steps. I strongly recommended doing this step by step. But if you wish to execute this from the commandline, all you need to provide is the UniProtKB ID and the total number of chains.

``` run.py --prepENS UniProt ID --multimeric Nr.ofChains ```

The issues pdbParser can't tackle at the moment are:
* Heteromeric structures in ensemble preperation.
* Structures with missing residues are not repaired. (An example modeller script can be found in the tutorial folder)

Below cases should give correct structures but caution is needed

* Support for multimeric assemblies 
* Insertions/deletions may result with partial structures. This a bit tricky to handle since the gap penalty is set to extreme.
