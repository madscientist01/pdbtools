#!/usr/bin/python
#
#  Test for the Superpose.py
#

from superpose import *
import argparse


if __name__ == "__main__":

	sup = Superpose(queryPDB='1ATN.pdb', subjectPDB='1K8K.pdb', superposedPDB='sup.pdb')
	sup.run()
	print sup.RMSD, sup.qscore, sup.aligned
	fasta = sup.FASTAoutput()
	f = open ('1ATN.fasta','w')
	f.writelines(fasta)
	f.close()