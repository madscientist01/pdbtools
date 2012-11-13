#!/usr/bin/python
#
# Wrapper of Superpose (CCP4)
#
#

from superpose import *
import argparse


if __name__ == "__main__":

	sup = Superpose(queryPDB='1ATN.pdb', subjectPDB='1K8K.pdb', superposedPDB='sup.pdb')
	sup.run()
	print sup.RMSD, sup.qscore, sup.aligned
	fasta = sup.FASTAoutput()

	for line in fasta:
		print line
	