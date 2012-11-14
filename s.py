#!/usr/bin/python
#
# Wrapper of Superpose (CCP4)
#
#

from superpose import *
import argparse


if __name__ == "__main__":

	sup = Superpose(queryPDB='1ATN.pdb', subjectPDB='1C0F.pdb', superposedPDB='sup.pdb')
	sup.run()
	print sup.RMSD, sup.qscore, sup.aligned
	for line in sup.FASTAoutput():
		print line
	