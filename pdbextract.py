#!/usr/bin/python
#
# PDBExtract.py
# General Purpose Protein Data Bank (PDB) file manipulation utility
#
# Written by MadScientist
# http://github.com/madscientist01
# http://madscientist.wordpress.com
#
#

import sys, subprocess, re, os, glob, argparse
import warnings
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

def PDBParse(pdb,filter,hetero):
		headers = []
		headerRecordList = ["HEADER","TITLE","COMPND","SOURCE","KEYWDS","EXPDTA","AUTHOR",\
							"REVDAT","JRNL","REMARK","DBREF","SEQRES","MODRES","HET","HETNAM",\
							"HETSYN","FORMUL","HELIX","SHEET","LINK","SITE","CRYST1","ORIGX1",\
							"ORIGX2","ORIGX3","SCALE1","SCALE2","SCALE3"]
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
						chains[chainid]=chainContent

					append = True
					if hetero and head == "HETATM":
						append = False

					if filter == line[17:20]:
						append = False

					if append :
						chains[chainid].append(line)

		return(headers,chains)

def fastaParsing (header):
#
# Parse SEQRES lines in PDB REMAKR and save it as one letter codes.
#
#
	oneLetter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
	'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
	'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
	'GLY':'G', 'PRO':'P', 'CYS':'C', 'MSE':'M'}
	seqres = {}
	chains = {}
	for line in header:
		if line[0:6] == "SEQRES":
			chainid = line[11:12]
			
			if not chainid in chains:
				chains[chainid] = ''
			chains[chainid] = chains[chainid]+line[19:70].strip()+' '

	for chain in chains:
		aa = chains[chain].strip().split(' ')
		aminoAcids = ''
		for a in aa:
			if a in oneLetter:
				aminoAcids = aminoAcids+oneLetter[a]
			else:
				aminoAcids = aminoAcids + 'X'
		seqres[chain] = aminoAcids
	return(seqres)

def fastaSplit(fasta,width):
#
# split sequence in fasta string into predefined width of string lists
#

	cursor = 0
	end = 0
	buffer = []
	while (cursor <len(fasta)):
		if (len(fasta)-cursor) > width:
			end+=width
		else:
			end=len(fasta)		
		buffer.append(fasta[cursor:end]+"\n")
		cursor+=width
	return(buffer)	

class PDBExtract(object):

	def __init__(self, **kwargs):

		self.files=kwargs['files']
		self.header=kwargs['header']
		self.list = kwargs['list']
		self.fasta=kwargs['fasta']
		self.extract=kwargs['extract']
		self.exclude=kwargs['exclude']
		self.split = kwargs['split']
		self.unique = kwargs['unique']
		self.filter = kwargs['filter']
		self.hetero = kwargs['hetero']

	def pdbProcess(self, filename):
	#
	# Parse pdb content and extract using regex
	#
		[header,chains] = PDBParse(filename, self.filter, self.hetero)
		print "processing {0}".format(filename)	
		newFileName=""

		if self.extract or self.exclude:
			f = open("temp.pdb",'w')
			if self.header:
				f.writelines(header)
			includedChain=""
			if self.extract:
				chainlist = self.extract.split()
				for chain in chainlist:
					if chain in chains:
						f.writelines(chains[chain])
						includedChain=includedChain+chain
			if self.exclude:
				excludelist = self.exclude.split()
				for chain in chains:
					if chain in chains and not chain in excludelist:
						f.writelines(chains[chain])
						includedChain=includedChain+chain
			f.close()
			if len(includedChain)>0:
				newFileName = filename[:len(filename)-4]+"_"+includedChain+".pdb"
				os.rename("temp.pdb",newFileName)
			else:
				os.remove("temp.pdb")

		if self.split:
			for chain in chains:
				newFileName =filename[:len(filename)-4]+"_"+chain+".pdb" 
				f = open(newFileName,'w')
				if self.header:
					f.writelines(header)
				f.writelines(chains[chain])
				f.close()

		if self.fasta:
			sequence = fastaParsing(header)
			newFileName =filename[:len(filename)-4]+".fasta" 
			f = open(newFileName,'w')
			for chain in sorted(sequence.iterkeys()):
				f.write(">{0}\n".format(chain))
				buffer = fastaSplit(sequence[chain],60)
				f.writelines(buffer)
			f.close()

		if self.unique:
			sequence = fastaParsing(header)
			sequences = []
			includeList = {}
			for chain in sequence:
				sequences.append(sequence[chain])
				includeList[chain] = True

			chainlist = sequence.keys()
			#
			# Using pairwise2 module in BioPython, carry out pairwise alignment with all of PDB chains.
			# if two alignment give score higher than length of aa*3, it is considered as redundant chain
			# and it will be removed.
			#
			if len(sequences)>1:
				gapOpen=-10
				gapExtend=-0.5
				matrix = matlist.blosum62
				for i in range(len(sequences)):
					for j in range(i+1, len(sequences)):
						aln = pairwise2.align.globalds(sequences[i],sequences[j],matrix,gapOpen,gapExtend)
						query,subject,score,begin,end = aln[0]
						cutoff = len(sequences[i])*3
						if cutoff < score:
							includeList[chainlist[j]] = False
	#					print chainlist[i],chainlist[j],len(sequences[i]),len(sequences[j]), score
			
			f = open("temp.pdb",'w')
			if self.header:
				f.writelines(header)
			
			includedChain = ''
			for chain in includeList:
				if includeList[chain]:
					f.writelines(chains[chain])
					includedChain=includedChain+chain
			f.close()

			if len(includedChain)>0:
				newFileName = filename[:len(filename)-4]+"_"+includedChain+".pdb"
				os.rename("temp.pdb",newFileName)
			else:
				os.remove("temp.pdb")

		if (self.hetero or self.filter) and (not self.exclude and not self.unique and not self.extract and not self.split):
			newFileName = filename[:len(filename)-4]+"_filtered"+".pdb"
			f = open(newFileName,'w')
			if self.header:
				f.writelines(header)
				for chain in chains:
					f.writelines(chains[chain])
			f.close()


		return (newFileName)

	def go(self):

		if len(self.files)>0 :
			pdbList=self.files
		else:
			pdbList=glob.glob('*.pdb')

		writtenFileList = []
		for pdb in pdbList:	
			if os.path.exists(pdb):
				filename = self.pdbProcess(pdb)
				if len(filename)>0:
					writtenFileList.append(filename+"\n")
			else:
				print "{0} is not exist!".format(pdb)

		if len(writtenFileList)>0:
			f = open (self.list,'w')
			f.writelines(writtenFileList)
			f.close

		return(writtenFileList)
	

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	group = parser.add_mutually_exclusive_group()
	parser.add_argument('-f', '--files', nargs='*', dest='files',default=[],
	                    help='Files to process')
	parser.add_argument('-t', '--no_header', action='store_false', dest='header', default=True,
						help='Strip PDB header')
	parser.add_argument('-l', '--list', action='store', dest='list', default='extractlist.txt',
						help='Saved List')
	parser.add_argument('-a', '--fasta',action='store_true', dest='fasta',default=False,
						help='Extract Fasta')
	parser.add_argument('-r', '--remove_residue', action='store', dest='filter', 
						help='Remove designated type of residues. if you want to remove waters, -r HOH')
	parser.add_argument('-ht', '--remove_heteroatom', action='store_true', dest='hetero', default=False,
						help='Remove heteroatom')
	group.add_argument('-e', '--extract', action='store', dest='extract', 
						help='Set chain to extract')
	group.add_argument('-x', '--exclude', action='store', dest='exclude',
						help='Set chain to exclude')	
	group.add_argument('-s', '--split', action='store_true', dest='split', default=False,
						help='Split all of chains to seperate PDB')
	group.add_argument('-u', '--unique', action='store_true', dest='unique', default=False,
						help='Extract only chains have uniq sequence')
	
	results = parser.parse_args()
	pdbExtract = PDBExtract(files=results.files, header=results.header, list=results.list, fasta=results.fasta,\
							extract=results.extract, exclude=results.exclude, split=results.split,\
							unique=results.unique, filter=results.filter, hetero=results.hetero)
	fileList = pdbExtract.go()