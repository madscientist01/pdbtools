#pdbtools

General Purpose Protein Data Bank (PDB) file manipulation utilities and library
Written and maintained by Suk Namgoong (suk.namgoong@gmail.com)

### Requirements

* [Superpose (CCP4)](http://www.ccp4.ac.uk/)
* [PyMol](http://www.pymol.org/)

### Files

* superpose.py : Wrapper for CCP4 superpose, which superimpose two PDB using Secondary Structure Matching
* pdbextract.py : Extract (or exclude) specific chains in multiple PDBs
* render.py : Batch rendering of multiple by PyMOL
* pdbquery.py : Using PDB RESTful web service, query PDB using multiple conditions and download them
* pdbblast.py : fetch pdb using BLAST search
# fetchPDB.py : fetch pdb by PDB accession codes

### Install

* download

    python setup.py install

### Usage

