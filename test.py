#!/usr/bin/python
#
#  Test for the Superpose.py
#

from superpose import *
import argparse


if __name__ == "__main__":

	sup = Superpose(queryPDB='1DBH.pdb', subjectPDB='1DYN.pdb', superposedPDB='sup.pdb')
	sup.run()
	print sup.RMSD, sup.qscore, sup.aligned
	fasta = sup.FASTAoutput()
	f = open ('1DBH.fasta','w')
	f.writelines(fasta)
	f.close()