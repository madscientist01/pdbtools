#!/usr/bin/python
#
#

from pdbextract import *

if __name__ == "__main__":

	[header,chains] = PDBParse('1BTN.pdb', None, None)
	# for chain in chains:
	# 	for atom in chains[chain]:
	# 		print atom

	extractRegion = "A:30-300"
	pdbExtract = PDBExtract(extractregion = extractRegion)
	newFileName = pdbExtract.extractRegions(header,chains,'1BTN.pdb')
	print newFileName
