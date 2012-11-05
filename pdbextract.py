#!/usr/bin/python
#
# PDBExtract
#
#
import sys, subprocess, re, os, glob, argparse
import warnings

def PDBParse(pdb,argument):
	headers = []
	headerRecordList = ["HEADER","TITLE","COMPND","SOURCE","KEYWDS","EXPDTA","AUTHOR","REVDAT","JRNL","REMARK","DBREF","SEQRES","MODRES","HET","HETNAM","HETSYN","FORMUL","HELIX","SHEET","LINK","SITE","CRYST1","ORIGX1","ORIGX2","ORIGX3","SCALE1","SCALE2","SCALE3"]
	dataRecordList = ["ATOM", "TER", "HETATM"]
	chains = {}
	currentChain = ''
	if os.path.exists(pdb):
		f = open(pdb)
		line = f.readline()
		while line:
			line = f.readline()
			head = line[0:6].strip()
			if head in headerRecordList:
				headers.append(line)
			if head in dataRecordList:
				chainid = line [21:22]
				if not chainid in chains:
					chainContent = []
					currentChain = chainid
					chains[currentChain]=chainContent	
				chains[currentChain].append(line)	
	return(headers,chains)

def pdbProcess(filename, argument):
#
# Parse pdb content and extract using regex
#
	[header,chains] = PDBParse(filename,argument)
	print "processing {0}".format(filename)

	if argument.extract or argument.exclude:
	
		f = open("temp.pdb",'w')
		if argument.header:
			f.writelines(header)
		includedChain=""
		if argument.extract:
			chainlist = argument.extract.split()
			for chain in chainlist:
				if chain in chains:
					f.writelines(chains[chain])
					includedChain=includedChain+chain

		if argument.exclude:
			excludelist = argument.exclude.split()
			for chain in chains:
				if chain in chains and not chain in excludelist:
					f.writelines(chains[chain])
					includedChain=includedChain+chain
		f.close()
		if len(includedChain)>0:
			os.rename("temp.pdb",filename[:len(filename)-4]+"_"+includedChain+".pdb")
		else:
			os.remove("temp.pdb")

	if argument.split:

		for chain in chains:
			f = open(filename[:len(filename)-4]+"_"+chain+".pdb",'w')
			if argument.header:
				f.writelines(header)
			f.writelines(chains[chain])
			f.close()
	return ()

def main(argument):

	if len(argument.files)>0 :
		pdbList=argument.files
	else:
		pdbList=glob.glob('*.pdb')

		for pdb in pdbList:	
			if os.path.exists(pdb):
				pdbProcess(pdb,argument)
				# pdbProcessing(pdb,argument)				
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
	parser.add_argument('-t', '--no_header', action='store_false', dest='header', default=True,
						help='Strip PDB header')


	results = parser.parse_args()
	main(results)