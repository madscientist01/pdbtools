#!/usr/bin/python
#
# PDBExtract.py
# General Purpose Protein Data Bank (PDB) file manipulation utility
#
# Written by MadScientist
# http://github.com/madscientist01
# http://madscientist.wordpress.com
# Require Python 2.6 or recent and BioPython
#

import sys, subprocess, re, os, glob, argparse
from shutil import copy
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

def readFasta(filename) :

	sequence = ''
	sequenceName = ''
	if os.path.exists(filename):
		
		f = open(filename)
		while True:
			line = f.readline().strip()
			if not line: 
				break 
			match = re.match('^>(.*)',line)
			if match:
				sequenceName = match.group(1)
			else:
				sequence = sequence+line			
		return(sequenceName,sequence)	    	
	else:
		return()

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

def filterATOM(chains,extractRegion):
	#
	# chains : PDB ATOM regions extracted by PDBParse
	# extractRegion : Dictionary contains extraction regions (start-end) and chain as key
	#
	buffer = []

	for chain in extractRegion:
		if chain in chains:
			s,e = extractRegion[chain].split("-")
			start=int(s)
			end = int(e)
			# print start,end
			for line in chains[chain]:
				residueNumber = int(line[22:26].strip())
				# print residueNumber,start,end
				if residueNumber >=start and residueNumber <=end:
					buffer.append(line)
	return (buffer)

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

def filterHeader(header,chainlist):
	#
	# split sequence in fasta string into predefined width of string lists
	#
	filteredHeader = []
	for line in header:
		if line[0:6] == "SEQRES":
			chainid = line[11:12]
			if chainid in chainlist:
				filteredHeader.append(line)
		else:
				filteredHeader.append(line)

	return(filteredHeader)	


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
	
		self.files=kwargs.get('files')
		self.fof = kwargs.get('fof')
		self.header=kwargs.get('header')
		self.list = kwargs.get('list')
		self.fasta=kwargs.get('fasta')
		self.extract=kwargs.get('extract')
		self.exclude=kwargs.get('exclude')
		self.split = kwargs.get('split')
		self.unique = kwargs.get('unique')
		self.filter = kwargs.get('filter')
		self.hetero = kwargs.get('hetero')
		self.excludeseq = kwargs.get('excludeseq')
		self.extractseq = kwargs.get('extractseq')
		self.overwrite = kwargs.get('overwrite')
		self.extractregion=kwargs.get('extractregion')
		if self.extractseq or self.excludeseq:
			if self.extractseq:
				filename = self.extractseq
			elif self.excludeseq:
				filename = self.excludeseq
			(self.querySeqFastaName, self.querySeqFasta) = readFasta(filename)
		return

	def newFileName(self,filename,includedChain):
		if self.overwrite:
			newFileName = filename
			copy(filename,filename+".backup")
		else:
			newFileName = filename[:len(filename)-4]+"_"+includedChain+".pdb"
		return(newFileName)

	def saveChains(self, chains,header,includeList,filename):

		buffer = []
		includedChain = ''
		for chain in sorted(includeList.iterkeys()):
			if includeList[chain]:
				buffer = buffer+chains[chain]
				includedChain=includedChain+chain
		if len(includedChain)>0:
			newFileName = self.newFileName(filename,includedChain)
			f=open(newFileName, 'w')
			f.writelines(filterHeader(header,includedChain.split()))
			f.writelines(buffer)
		else:
			newFileName=''

		return (newFileName)
	
	def extractRegions(self,header,chain,filename):
		
		regionDictionary={}
		regions = self.extractregion.split()
		includedChain = ""
		extractRegex = re.compile("(\S):(\d+)-(\d+)")
		for region in regions:
			match = extractRegex.match(region)
			if match:
				start = int(match.group(2))
				end = int(match.group(3))
				if end<start:
					start = match.group(3)
					end = match.group(2)
				regionDictionary[match.group(1)]="{0}-{1}".format(start,end)
				includedChain=includedChain+match.group(1)
		buffer = filterATOM(chain,regionDictionary)

		if len(buffer)>0:

			newFileName = self.newFileName(filename,"".join(regionDictionary))
			f = open(newFileName,'w')
			if self.header:
				filteredHeader = filterHeader(header,includedChain.split())
				f.writelines(filteredHeader)
			f.writelines(buffer)
			f.close
		else:
			newFileName = ""

		return(newFileName)

	def extractExclude(self,header,chains,filename):

		buffer = []
		includedChain=""
		if self.extract:
			chainlist = self.extract.split()
			for chain in chainlist:
				if chain in chains:
					buffer=buffer+chains[chain]
					includedChain=includedChain+chain

		if self.exclude:
			excludelist = self.exclude.split()
			for chain in chains:
				if chain in chains and not chain in excludelist:
					buffer=buffer+chains[chain]
					includedChain=includedChain+chain

		if len(includedChain)>0:
			newFileName = self.newFileName(filename,includedChain)		
			f = open(newFileName,'w')
			if self.header:
				filteredHeader = filterHeader(header,includedChain.split())
				f.writelines(filteredHeader)
			f.writelines(buffer)
			f.close
		else:
			newFileName = ''
		return (newFileName)


	def splitChains(self,header,chains,filename):

		for chain in chains:
			newFileName =filename[:len(filename)-4]+"_"+chain+".pdb" 
			f = open(newFileName,'w')
			if self.header:
				filteredHeader = filterHeader(header,[chain])
				f.writelines(filteredHeader)
			f.writelines(chains[chain])
			f.close()
		return (newFileName)

	def fastaSeq(self,header,filename):

		sequence = fastaParsing(header)
		newFileName =filename[:len(filename)-4]+".fasta" 
		f = open(newFileName,'w')
		for chain in sorted(sequence.iterkeys()):
			f.write(">{0}\n".format(chain))
			buffer = fastaSplit(sequence[chain],60)
			f.writelines(buffer)
		f.close()
		return (newFileName)

	def extractBasedOnSeq (self,header,chains,filename):

		sequence = fastaParsing(header)
		if len(sequence)>1:
			includeList = {}
			gapOpen=-10
			gapExtend=-0.5
			matrix = matlist.blosum62

			for chain in sequence:
				
				aln = pairwise2.align.globalds(self.querySeqFasta,sequence[chain],matrix,gapOpen,gapExtend)
				query,subject,score,begin,end = aln[0]
				cutoff = len(self.querySeqFasta)*4
				if cutoff < score :
					if self.extractseq:
						includeList[chain]=True
					if self.excludeseq:
						includeList[chain]=False	
				else: 
					if self.extractseq:
						includeList[chain]=False
					if self.excludeseq:
						includeList[chain]=True
			
			newFileName = self.saveChains(chains,header,includeList,filename)
		else:
			newFileName=''
		return (newFileName)

	def uniqueChain (self,header,chains,filename):
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
		newFileName = self.saveChains(chains,header,includeList,filename)
		return(newFileName)

	def filtering (self,header,chains,filename):
		
		newFileName = self.newFileName(filename,'filtered')
		f = open(newFileName,'w')
		if self.header:
			f.writelines(header)
			[f.writelines(chains[chain]) for chain in chains]
		f.close(newFileName)


	def pdbProcess(self, filename):

	#
	# Parse pdb content and extract using regex
	#
		[header,chains] = PDBParse(filename, self.filter, self.hetero)
		print "processing {0}".format(filename)	
		newFileName=""

		if self.extract or self.exclude:
			newFileName=self.extractExclude(header,chains,filename)

		if self.extractregion:
			newFileName=self.extractRegions(header,chains,filename)

		if self.split:
			newFileName=self.splitChains(header,chains,filename)

		if self.fasta:
			newFileName=self.fastaSeq(header,filename)

		if self.extractseq or self.excludeseq:
			newFileName=self.extractBasedOnSeq(header,chains,filename)

		if self.unique:		
			newFileName=self.uniqueChain(header,chains,filename)
		
		if (self.hetero or self.filter) and (not self.exclude and not self.unique and not self.extract and not self.split):
			newFileName=self.filtering(header,chains,filename)

		return (newFileName)

	def go(self):

		if self.files:
			pdbList=self.files
			
		elif self.fof and os.path.exists(self.fof):
			f = open (self.fof)
			pdbList=f.readlines()
			f.close()
			pdbList = [x.strip() for x in pdbList]
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
	filegroup = parser.add_mutually_exclusive_group()
	filegroup.add_argument('-f', '--files', nargs='*', dest='files',default=[],
	                    help='Files to process')
	filegroup.add_argument('-T', '--files-from',  dest='fof',
	                    help='Files from file of names')
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
	parser.add_argument('-o', '--overwrite_original', action='store_true', dest='overwrite', default=False,
						help='Remove heteroatom')
	group.add_argument('-e', '--extract', action='store', dest='extract', 
						help='Set chain to extract')
	group.add_argument('-R', '--extract_region', action='store', dest='extractregion', 
						help='Set range of extraction. Example: A:20-100 B10-40')
	group.add_argument('-x', '--exclude', action='store', dest='exclude',
						help='Set chain to exclude')
	group.add_argument('-E', '--extract_sequence', action='store', dest='extractseq', 
						help='Extract PDB chains based on the protein sequences')
	group.add_argument('-X', '--exclude_sequence', action='store', dest='excludeseq',
						help='Exclude PDB chains based on the protein sequences')
	group.add_argument('-s', '--split', action='store_true', dest='split', default=False,
						help='Split all of chains to seperate PDB')
	group.add_argument('-u', '--unique', action='store_true', dest='unique', default=False,
						help='Extract only unique PDB chains ')
	
	results = parser.parse_args()
	pdbExtract = PDBExtract(files=results.files, fof=results.fof, header=results.header, list=results.list,\
							fasta=results.fasta, extract=results.extract, extractregion = results.extractregion,\
							overwrite=results.overwrite, exclude=results.exclude,\
							split=results.split, unique=results.unique, filter=results.filter,\
							hetero=results.hetero, extractseq=results.extractseq, excludeseq=results.excludeseq)
	fileList = pdbExtract.go()