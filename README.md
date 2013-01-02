#PDBTOOLS

##General Purpose PDB manipulation utilities

Written and maintained by Suk Namgoong (suk.namgoong@gmail.com)


# Requirements

* [Biopython 1.6.0](https://github.com/biopython/biopython)
* [Superpose (CCP4)](http://www.ccp4.ac.uk/)
* [PyMol](http://www.pymol.org/)

# Files

* superpose.py : Wrapper for CCP4 superpose, which superimpose two PDB using Secondary Structure Matching
* pdbextract.py : Extract (or exclude) specific chains in multiple PDBs
* render.py : Batch rendering of multiple by PyMOL
* pdbquery.py : Using PDB RESTful web service, query PDB using multiple conditions and download them
* pdbblast.py : fetch pdb using BLAST search

