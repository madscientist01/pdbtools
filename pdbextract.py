#!/usr/bin/python
#
# PDBExtract
#
#

import sys, subprocess, re, os, glob, argparse
import warnings
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.PDBIO import Select




class ChainSelect(Select):
	def __init__(self,targetchain):
		self.targetchain = targetchain

	def accept_chain(self, chain):
#       print dir(residue)
		if chain.get_id() in self.targetchain:
			return 1
		else:
			return 0

def PDBHeader(pdb):
	header = []
	idcapture = re.compile("^(\S+)")
	capturelist = ["HEADER","TITLE","COMPND","SOURCE","KEYWDS","EXPDTA","AUTHOR","REVDAT","JRNL","REMARK","DBREF","SEQRES","MODRES","HET","HETNAM","HETSYN","FORMUL","HELIX","SHEET","LINK","SITE","CRYST1","ORIGX1","ORIGX2","ORIGX3","SCALE1","SCALE2","SCALE3"]
	if os.path.exists(pdb):
		f = open(pdb)
		line = f.readline()
		while line:
			line = f.readline()
			match = idcapture.match(line)
			if match and match.group(1) in capturelist:
				header.append(line)
	return(header)


def pdbParsing(filename, argument):
#
# Parse pdb content and extract using regex
#
	warnings.filterwarnings('ignore')
	parser = PDBParser()
	structure = parser.get_structure('self', filename)
	header = PDBHeader(filename)
	model  = structure[0]
	print "processing {0}".format(filename)
	if argument.extract:
		chainlist = argument.extract.split()
	if argument.exclude:
		chainlist = []
		excludelist = argument.exclude.split()
		for chain in model:
			if not chain.get_id() in excludelist:
				chainlist.append(chain.get_id())
	if argument.split:
		splitlist = ''
		for chain in model:
			splitlist= splitlist+chain.get_id()

	if argument.extract or argument.exclude:
		if len(chainlist)>0:
			w = PDBIO()
			w.set_structure(structure)
			w.save("temp.pdb",select=ChainSelect(chainlist))
			f = open(filename[:len(filename)-4]+"_"+"".join(chainlist)+".pdb",'w')
			f.writelines(header)
			t = open("temp.pdb")
			f.writelines(t.readlines())
			t.close()
			os.remove("temp.pdb")
			f.close()
					
	if argument.split:
		for chain in splitlist:
			w = PDBIO()
			w.set_structure(structure)
			templist=[]
			templist.append(chain)
			w.save("temp.pdb",select=ChainSelect(templist))
			f = open(filename[:len(filename)-4]+"_"+chain+".pdb",'w')
			f.writelines(header)
			t = open("temp.pdb")
			f.writelines(t.readlines())
			t.close()
			os.remove("temp.pdb")
			f.close()
	
	return ()

def main(argument):

	if len(argument.files)>0 :
		pdbList=argument.files
	else:
		pdbList=glob.glob('*.pdb')

		for pdb in pdbList:	
			if os.path.exists(pdb):
				pdbParsing(pdb,argument)				
			else:
				print "{0} is not exist!".format(pdb)
		sys.exit()		


if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--files', nargs='*', dest='files',default=[],
	                    help='Files to process')
	parser.add_argument('-e', '--extract', action='store', dest='extract', 
						help='Set chain to extract')
	parser.add_argument('-x', '--exclude', action='store', dest='exclude',
						help='Set chain to exclude')	
	parser.add_argument('-s', '--split', action='store_true', dest='split', default=False,
						help='Split all of chains to seperate PDB')

	results = parser.parse_args()
	main(results)