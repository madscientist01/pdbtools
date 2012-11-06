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
					chains[chainid]=chainContent	
				chains[chainid].append(line)	
	return(headers,chains)

def fastaParsing(header):
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


def pdbProcess(filename, argument):
#
# Parse pdb content and extract using regex
#
	[header,chains] = PDBParse(filename,argument)
	print "processing {0}".format(filename)	
	newFileName=""
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
			newFileName = filename[:len(filename)-4]+"_"+includedChain+".pdb"
			os.rename("temp.pdb",newFileName)
		else:
			os.remove("temp.pdb")

	if argument.split:

		for chain in chains:
			newFileName =filename[:len(filename)-4]+"_"+chain+".pdb" 
			f = open(newFileName,'w')
			if argument.header:
				f.writelines(header)
			f.writelines(chains[chain])
			f.close()

	if argument.fasta:
		sequence = fastaParsing(header)
		newFileName =filename[:len(filename)-4]+".fasta" 
		f = open(newFileName,'w')
		for chain in sorted(sequence.iterkeys()):
			f.write(">{0}\n".format(chain))
			buffer = fastaSplit(sequence[chain],60)
			f.writelines(buffer)
		f.close()
	
	return (newFileName)

def main(argument):

	if len(argument.files)>0 :
		pdbList=argument.files
	else:
		pdbList=glob.glob('*.pdb')


	writtenFileList = []
	for pdb in pdbList:	
		if os.path.exists(pdb):

			filename = pdbProcess(pdb,argument)
			if len(filename)>0:
				writtenFileList.append(filename+"\n")
		else:
			print "{0} is not exist!".format(pdb)

	if len(writtenFileList)>0:
		f = open (argument.list,'w')
		f.writelines(writtenFileList)
		f.close

	sys.exit()		


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
	group.add_argument('-e', '--extract', action='store', dest='extract', 
						help='Set chain to extract')
	group.add_argument('-x', '--exclude', action='store', dest='exclude',
						help='Set chain to exclude')	
	group.add_argument('-s', '--split', action='store_true', dest='split', default=False,
						help='Split all of chains to seperate PDB')


	results = parser.parse_args()
	main(results)